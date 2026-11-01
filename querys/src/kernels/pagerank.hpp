#pragma once
// DQ-5: exact synchronous PageRank on a fixed, lossless summary.
// Input storage is a single-copy triangular CSR, not a query-specific index.
// All degree/reciprocal/iteration state is built inside each query call.
// The same implementation handles positive and signed summaries.

#include "../upstream/storage.hpp"
#include <memory>
#include <type_traits>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#define DQ5_X86 1
#else
#define DQ5_X86 0
#endif

namespace dq5 {
using dq3::VI; using dq3::VD; using dq3::Edges; using dq3::Edge;
using dq3::i64; using dq3::require;
enum class Backend { Automatic, Scalar, AVX512 };

// One copy (u,v), u<=v, with the source implicit in the CSR row.
struct HalfCSR {
    VI off, to;
    HalfCSR() = default;
    HalfCSR(int n, const Edges& edges, bool diagonal) {
        require(n>0 && edges.size()<=size_t(INT32_MAX),"half CSR size");
        off.assign(n+1,0);
        Edge prev{-1,-1};
        for(auto [u,v]:edges) {
            require(0<=u && u<n && u<=v && v<n &&
                    (diagonal || u<v),"half CSR endpoint");
            require(prev<Edge{u,v},"half CSR needs sorted unique records");
            ++off[u+1]; prev={u,v};
        }
        std::partial_sum(off.begin(),off.end(),off.begin());
        to.reserve(edges.size());
        for(auto [u,v]:edges) to.push_back(v);
    }



};

// Ordinary summary storage. No stored degree, reciprocal, completeness flag,
// profiles, query graph, or precomputed answers.
struct Summary {
    int n=0,b=0; i64 m=0;
    VI owner, size;
    HalfCSR super, plus, minus;
    Summary()=default;
    Summary(int nn,i64 mm,VI own,const Edges& se,const Edges& cp,const Edges& cm)
        :n(nn),m(mm),owner(std::move(own)) {
        require(n>0 && owner.size()==size_t(n) && m>=0,"summary dimensions");
        for(int a:owner){require(a>=0&&a<n,"owner range");b=std::max(b,a+1);}
        size.assign(b,0); for(int a:owner)++size[a];
        for(int x:size)require(x>0,"empty supernode");
        super=HalfCSR(b,se,true); plus=HalfCSR(n,cp,false);
        if(!cm.empty())minus=HalfCSR(n,cm,false);
    }
    explicit Summary(const dq3::Summary& s):
        Summary(s.n,s.m,s.owner,s.raw_se,s.raw_cp,s.raw_cm) {
        require(!s.packed,"DQ5 preserves input vertex IDs; pass plain summary");
    }
    i64 ns()const{return super.to.size();}
    i64 np()const{return plus.to.size();}


    i64 nq()const{return minus.to.size();}




};

inline bool has_avx512() {
#if DQ5_X86 && (defined(__GNUC__) || defined(__clang__))
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f");
#else
    return false;
#endif
}
inline bool use_vector(Backend mode) {
    if(mode==Backend::Scalar)return false;
    bool ok=has_avx512();
    if(mode==Backend::AVX512)require(ok,"AVX512 requested but not supported");
    return ok;
}

struct Operations {
    uint64_t degree_super_records=0, degree_correction_records=0;
    uint64_t iterations=0, super_records=0, positive_records=0, negative_records=0;
    uint64_t nonempty_positive_rows=0, nonempty_negative_rows=0;
    uint64_t source_writes=0, target_writes=0, vertex_normalizations=0;
    uint64_t mass_accumulations=0, base_writes=0;
    uint64_t vector_chunks=0, scalar_edges=0;
};

template<bool Negative,bool Count=false>
inline void bilateral_scalar(const HalfCSR& rows, int n,
        const double* __restrict__ z,double* __restrict__ next,Operations*o=nullptr) {
    const int* __restrict__ off=rows.off.data();
    const int* __restrict__ to=rows.to.data();
    for(int u=0;u<n;++u) {
        int j=off[u],end=off[u+1];
        if(j==end)continue;
        // Start with the existing accumulated value, not with zero:
        // exactly 2 floating point edge updates, no final extra source addition.
        double acc=next[u];const double outgoing=z[u];
        for(;j<end;++j) {
            int v=to[j];
            if constexpr(Negative){acc-=z[v];next[v]-=outgoing;}
            else {acc+=z[v];next[v]+=outgoing;}
        }
        next[u]=acc;
        if constexpr(Count){
            uint64_t e=end-off[u];
            o->source_writes++;o->target_writes+=e;o->scalar_edges+=e;
            if constexpr(Negative){o->negative_records+=e;o->nonempty_negative_rows++;}
            else{o->positive_records+=e;o->nonempty_positive_rows++;}
        }
    }
}
#if DQ5_X86 && (defined(__GNUC__) || defined(__clang__))
// Eight distinct target IDs within each sorted row: SIMD scatter is conflict-free.
// This is NOT safe across unrelated rows; vectorization never crosses a row.
template<bool Negative,bool Count=false>
__attribute__((target("avx512f")))
inline void bilateral_avx512(const HalfCSR& rows,int n,
        const double* __restrict__ z,double* __restrict__ next,Operations*o=nullptr) {
    const int* __restrict__ off=rows.off.data();
    const int* __restrict__ to=rows.to.data();
    for(int u=0;u<n;++u) {
        int j=off[u],end=off[u+1];if(j==end)continue;
        const double outgoing=z[u];double acc=next[u];
        uint64_t chunks=0;
        if(end-j>=8) {
            __m512d total=_mm512_setzero_pd(),broadcast=_mm512_set1_pd(outgoing);
            for(;j+8<=end;j+=8) {
                __m256i ids=_mm256_loadu_si256((const __m256i*)(to+j));
                __m512d values=_mm512_i32gather_pd(ids,z,8);
                __m512d old=_mm512_i32gather_pd(ids,next,8);
                if constexpr(Negative)old=_mm512_sub_pd(old,broadcast);
                else old=_mm512_add_pd(old,broadcast);
                _mm512_i32scatter_pd(next,ids,old,8);
                total=_mm512_add_pd(total,values);
                ++chunks;
            }
            double r=_mm512_reduce_add_pd(total);
            if constexpr(Negative)acc-=r;else acc+=r;
        }
        int remainder=end-j;
        for(;j<end;++j) {
            int v=to[j];
            if constexpr(Negative){acc-=z[v];next[v]-=outgoing;}
            else {acc+=z[v];next[v]+=outgoing;}
        }
        next[u]=acc;
        if constexpr(Count){
            uint64_t e=end-off[u];o->vector_chunks+=chunks;o->scalar_edges+=remainder;
            o->source_writes++;o->target_writes+=e;
            if constexpr(Negative){o->negative_records+=e;o->nonempty_negative_rows++;}
            else{o->positive_records+=e;o->nonempty_positive_rows++;}
        }
    }
}
#endif
template<bool Negative,bool Count=false>
inline void bilateral(const HalfCSR&r,int n,const double*z,double*x,bool simd,Operations*o) {
#if DQ5_X86 && (defined(__GNUC__) || defined(__clang__))
    if(simd){bilateral_avx512<Negative,Count>(r,n,z,x,o);return;}
#endif
    bilateral_scalar<Negative,Count>(r,n,z,x,o);
}

// Main kernel: Synchronous Half-Edge Register Accumulation (SHERA).
// No fixed dataset-dependent threshold. Eight is the number of double lanes
// in a 512-bit register; small rows use the scalar tail.
template<bool Signed,bool Count=false>
VD pagerank(const Summary&s,int T=20,double alpha=.85,
            Backend backend=Backend::Automatic,Operations* o=nullptr) {
    require(T>=0&&alpha>=0&&alpha<1&&std::isfinite(alpha),"PageRank parameters");
    if constexpr(!Signed)require(s.nq()==0,"positive kernel cannot read negative summary");
    if constexpr(Count)require(o!=nullptr,"counter pointer required");
    const int n=s.n,b=s.b;
    if(T==0)return VD(n,1./n);
    const bool simd=use_vector(backend);
    VI base(b,0),degree(n);
    std::vector<unsigned char> diag(b,0);
    for(int a=0;a<b;++a) {
        for(int j=s.super.off[a];j<s.super.off[a+1];++j) {
            int v=s.super.to[j];
            if(a==v){base[a]+=s.size[a]-1;diag[a]=1;}
            else {base[a]+=s.size[v];base[v]+=s.size[a];}
        }
    }
    for(int u=0;u<n;++u)degree[u]=base[s.owner[u]];
    for(int u=0;u<n;++u)for(int j=s.plus.off[u];j<s.plus.off[u+1];++j) {
        ++degree[u];++degree[s.plus.to[j]];
    }
    if constexpr(Signed)if(s.nq())for(int u=0;u<n;++u)
        for(int j=s.minus.off[u];j<s.minus.off[u+1];++j) {
            --degree[u];--degree[s.minus.to[j]];
        }
    VD inv(n),z(n),next(n),mass(b,0),incoming(b);
    int n0=0;
    for(int u=0;u<n;++u) {
        require(degree[u]>=0&&degree[u]<n,"invalid summary degree");
        if(degree[u])inv[u]=alpha/degree[u];else ++n0;
        z[u]=inv[u]/n;mass[s.owner[u]]+=z[u];
    }
    if constexpr(Count){
        o->degree_super_records+=s.ns();
        o->degree_correction_records+=s.np()+s.nq();
        o->vertex_normalizations+=n;o->mass_accumulations+=n;
    }
    double sink=double(n0)/n;
    for(int t=0;t<T;++t) {
        std::fill(incoming.begin(),incoming.end(),0.);
        for(int a=0;a<b;++a) {
            int j=s.super.off[a],end=s.super.off[a+1];
            if(j==end)continue;
            double acc=incoming[a],out=mass[a];
            if(s.super.to[j]==a){acc+=out;++j;}
            for(;j<end;++j){int v=s.super.to[j];acc+=mass[v];incoming[v]+=out;}
            incoming[a]=acc;
        }
        const double tele=(1.-alpha+alpha*sink)/n;
        for(int u=0;u<n;++u) {
            int a=s.owner[u];
            next[u]=tele+incoming[a]-(diag[a]?z[u]:0.);
        }
        bilateral<false,Count>(s.plus,n,z.data(),next.data(),simd,o);
        if constexpr(Signed)if(s.nq())
            bilateral<true,Count>(s.minus,n,z.data(),next.data(),simd,o);
        if constexpr(Count){o->iterations++;o->super_records+=s.ns();o->base_writes+=n;}
        if(t+1==T)break;
        std::fill(mass.begin(),mass.end(),0.);
        for(int u=0;u<n;++u){z[u]=inv[u]*next[u];mass[s.owner[u]]+=z[u];}
        if constexpr(Count){o->vertex_normalizations+=n;o->mass_accumulations+=n;}
        sink=n0*tele; // undirected graph, uniform initial/teleport vector
    }
    return next;
}
} // namespace dq5
