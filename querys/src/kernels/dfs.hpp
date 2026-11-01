#pragma once
#include "traversal.hpp"
#include <memory>
namespace dq4 {
using namespace dq3;
struct Counts {
    uint64_t vertex_entries=0, parent_returns=0, root_tests=0;
    uint64_t plus_probes=0, minus_probes=0, csr_probes=0;
    uint64_t super_probes=0, member_tests=0, head_advances=0, super_retired=0;
    uint64_t super_tree_edges=0, plus_tree_edges=0, zipper_links=0, zipper_chains=0;
    uint64_t early_unwind_vertices=0;
    uint64_t pr_super_records=0, pr_positive_records=0, pr_negative_records=0;
    uint64_t pr_source_groups=0, pr_target_writes=0, pr_source_writes=0;
    uint64_t pr_vertex_passes=0, pr_owner_scatter=0;
};


template<bool Count=false, bool Early=true>
DFS dfs_spine(const Summary& s, Counts* stats=nullptr,int src=0) {
    require(s.packed && s.nq==0, "dfs_spine needs block-contiguous positive summary");
    DFS out; out.parent.assign(s.n,-2);out.tin.resize(s.n);out.tout.resize(s.n);
    VI head(s.moff.begin(),s.moff.end()-1);
    VI cursor(s.super.off.begin(),s.super.off.end()-1);
    std::unique_ptr<int[]> cp(new int[s.n]); // initialized only on first entry
    int clock=0, remaining=s.n;
    int* const parent=out.parent.data();
    auto enter=[&](int v,int p) {
        parent[v]=p; out.tin[v]=clock++; cp[v]=s.plus.off[v]; --remaining;
        if constexpr(Count) ++stats->vertex_entries;
    };
    auto first=[&](int block)->int {
        int& h=head[block]; const int end=s.moff[block+1];
        while(h<end) {
            if constexpr(Count) ++stats->member_tests;
            if(parent[h]==-2) return h;
            ++h;
            if constexpr(Count) ++stats->head_advances;
        }
        return -1;
    };
    for(int id=src,seen=0;seen<s.n;++seen) {
        if constexpr(Count) ++stats->root_tests;
        int current=s.query_id(id);
        if(++id==s.n)id=0;
        if(parent[current]!=-2) continue;
        enter(current,-1);
        while(current>=0) {
            if constexpr(Early) {
                if(remaining==0) {
                    // All possible neighbors of all active frames are visited.
                    do {
                        out.tout[current]=clock++;
                        if constexpr(Count) {++stats->parent_returns;++stats->early_unwind_vertices;}
                        current=parent[current];
                    } while(current>=0);
                    return restore_dfs(s,std::move(out));
                }
            }
            const int a=s.owner[current];
            int target=-1, v=-1;
            int& edge=cursor[a];
            while(edge<s.super.off[a+1]) {
                target=s.super.to[edge];
                if constexpr(Count) ++stats->super_probes;
                v=first(target);
                if(v>=0) break;
                ++edge;
                if constexpr(Count) ++stats->super_retired;
            }
            if(v>=0) {
                if constexpr(Count) {++stats->super_tree_edges;++stats->zipper_chains;}
                enter(v,current);current=v;
                // No per-frame preferred block; a,b are known until one side
                // has no unvisited member. Each step is a genuine recursive call.
                int next_block=a, back_block=target;
                while(true) {
                    int y=first(next_block);
                    if(y<0) break;
                    if constexpr(Count) {++stats->super_tree_edges;++stats->zipper_links;}
                    enter(y,current);current=y;
                    std::swap(next_block,back_block);
                }
                continue;
            }
            int& j=cp[current];
            while(j<s.plus.off[current+1]) {
                int x=s.plus.to[j++];
                if constexpr(Count) ++stats->plus_probes;
                if(parent[x]==-2) {v=x;break;}
            }
            if(v>=0) {
                if constexpr(Count) ++stats->plus_tree_edges;
                enter(v,current);current=v;
            } else {
                out.tout[current]=clock++;
                if constexpr(Count) ++stats->parent_returns;
                current=parent[current];
            }
        }
    }
    return restore_dfs(s,std::move(out));
}


template<bool Count=false, bool Early=true>
DFS dfs_csr_spine(const Graph&g,Counts* stats=nullptr,int src=0) {
    DFS out;out.parent.assign(g.n,-2);out.tin.resize(g.n);out.tout.resize(g.n);
    std::unique_ptr<int[]> next(new int[g.n]);int clock=0,remaining=g.n;
    auto enter=[&](int v,int p){out.parent[v]=p;out.tin[v]=clock++;next[v]=g.adj.off[v];--remaining;if constexpr(Count)++stats->vertex_entries;};
    for(int root=src,seen=0;seen<g.n;++seen) {
        if constexpr(Count) ++stats->root_tests;
        int candidate=root;
        if(++root==g.n)root=0;
        if(out.parent[candidate]!=-2)continue;
        int current=candidate;enter(current,-1);
        while(current>=0) {
            if constexpr(Early) if(remaining==0) {
                do {out.tout[current]=clock++;if constexpr(Count){++stats->parent_returns;++stats->early_unwind_vertices;}current=out.parent[current];}while(current>=0);
                return out;
            }
            int v=-1;int&j=next[current];
            while(j<g.adj.off[current+1]) {
                int x=g.adj.to[j++];if constexpr(Count)++stats->csr_probes;
                if(out.parent[x]==-2){v=x;break;}
            }
            if(v>=0){enter(v,current);current=v;}
            else {out.tout[current]=clock++;if constexpr(Count)++stats->parent_returns;current=out.parent[current];}
        }
    }
    return out;
}

}
