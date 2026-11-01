#pragma once
#include "traversal.hpp"
#include <limits>
namespace dq15 {
using namespace dq3;
struct ColoringStats {
  uint64_t se=0,plus=0,minus=0,member=0,neg=0,exc=0;
  uint64_t vertices=0,csr=0,earlier=0,neg_reject=0,mark_writes=0;
  uint64_t color_tests=0,word_tests=0,block_updates=0;
  uint64_t exception_adds=0,exception_comparisons=0,exception_new_pairs=0;
  uint64_t max_block_count=0,max_exception_count=0,state_bytes=0;
};


template<bool Count=false>
VI coloring_csr(const Graph&g,ColoringStats&st) {
  VI c(g.n,-1),stamp(128,-1);
  for(int u=0;u<g.n;++u) {
    if constexpr(Count)++st.vertices;
    for(int j=g.adj.off[u];j<g.adj.off[u+1];++j) {
      if constexpr(Count)++st.csr;
      int v=g.adj.to[j];
      if(v<u) {
        int k=c[v];require(k>=0,"CSR coloring uncolored predecessor");
        if(size_t(k)>=stamp.size())stamp.resize(size_t(k)+1,-1);
        stamp[k]=u;if constexpr(Count)++st.earlier;
      }
    }
    int x=0;
    while(size_t(x)<stamp.size() && stamp[x]==u){if constexpr(Count)++st.color_tests;++x;}
    if constexpr(Count)++st.color_tests;
    if(size_t(x)==stamp.size())stamp.push_back(-1);
    c[u]=x;
  }
  return c;
}


template<bool Signed,bool Count=false>
VI coloring_expand(const Summary&s,ColoringStats&st) {
  require(s.packed,"coloring expand requires packed summary");
  if constexpr(!Signed)require(s.nq==0,"positive coloring cannot read C-");
  VI col(s.n,-1),stamp(128,-1),mark;
  if constexpr(Signed)mark.assign(s.n,-1);
  for(int orig=0;orig<s.n;++orig) {
    int u=s.query_id(orig);if constexpr(Count)++st.vertices;
    if constexpr(Signed)for(int j=s.minus.off[u];j<s.minus.off[u+1];++j) {
      mark[s.minus.to[j]]=u;
      if constexpr(Count){++st.minus;++st.mark_writes;}
    }
    auto take=[&](int v){
      if(s.input_id(v)<orig){
        int k=col[v];require(k>=0,"expand coloring uncolored predecessor");
        if(size_t(k)>=stamp.size())stamp.resize(size_t(k)+1,-1);
        stamp[k]=orig;if constexpr(Count)++st.earlier;
      }
    };
    for(int j=s.plus.off[u];j<s.plus.off[u+1];++j) {
      if constexpr(Count)++st.plus;take(s.plus.to[j]);
    }
    int a=s.owner[u];
    for(int j=s.super.off[a];j<s.super.off[a+1];++j) {
      int b=s.super.to[j];if constexpr(Count)++st.se;
      for(int v=s.moff[b];v<s.moff[b+1];++v) {
        if(v==u)continue;
        if constexpr(Count)++st.member;
        if constexpr(Signed){
          if constexpr(Count)++st.neg;
          if(mark[v]==u){if constexpr(Count)++st.neg_reject;continue;}
        }
        take(v);
      }
    }
    int k=0;
    while(size_t(k)<stamp.size()&&stamp[k]==orig){if constexpr(Count)++st.color_tests;++k;}
    if constexpr(Count)++st.color_tests;
    if(size_t(k)==stamp.size())stamp.push_back(-1);
    col[u]=k;
  }
  VI out(s.n);for(int u=0;u<s.n;++u)out[s.input_id(u)]=col[u];return out;
}


template<bool Count=false>
VI coloring_core_block(const Summary&s,ColoringStats&st) {
  require(s.packed&&s.nq==0,"core block coloring positive packed");
  // Allocate bitmap words only for colors actually encountered in a block.
  std::vector<std::vector<uint64_t>> bm(s.b);VI c(s.n,-1);
  if constexpr(Count)st.state_bytes=uint64_t(s.b)*sizeof(bm[0]);
  for(int orig=0;orig<s.n;++orig) {
    int u=s.query_id(orig),b=s.owner[u];auto m=bm[b];
    if constexpr(Count)++st.vertices;
    for(int j=s.plus.off[u];j<s.plus.off[u+1];++j) {
      int v=s.plus.to[j];if constexpr(Count)++st.plus;
      if(s.input_id(v)<orig){
        int k=c[v];require(k>=0,"block coloring uncolored predecessor");
        if(size_t(k>>6)>=m.size())m.resize(size_t(k>>6)+1,0);
        m[k>>6]|=uint64_t(1)<<(k&63);if constexpr(Count)++st.earlier;
      }
    }
    int col=-1;
    for(size_t w=0;w<m.size();++w) {
      if constexpr(Count)++st.word_tests;
      uint64_t z=~m[w];if(z){col=int(w*64+__builtin_ctzll(z));break;}
    }
    if(col<0)col=int(m.size()*64); // first bit of the first unused word
    require(col<s.n,"color index must be smaller than number of vertices");
    c[u]=col;
    for(int j=s.super.off[b];j<s.super.off[b+1];++j) {
      auto &row=bm[s.super.to[j]];
      if(size_t(col>>6)>=row.size())row.resize(size_t(col>>6)+1,0);
      row[col>>6]|=uint64_t(1)<<(col&63);
      if constexpr(Count){++st.se;++st.block_updates;}
    }
  }
  VI out(s.n);for(int u=0;u<s.n;++u)out[s.input_id(u)]=c[u];return out;
}


template<class Counter=uint32_t>
struct ColorExc {
  std::vector<std::pair<int,Counter>> x;
  template<bool Count>void add(int c,ColoringStats&st) {
    if constexpr(Count)++st.exception_adds;
    for(auto&p:x) {
      if constexpr(Count)++st.exception_comparisons;
      if(p.first==c) {
        require(p.second<std::numeric_limits<Counter>::max(),"block correction counter overflow");
        ++p.second;
        if constexpr(Count)st.max_exception_count=std::max(st.max_exception_count,uint64_t(p.second));
        return;
      }
    }
    x.push_back({c,Counter(1)});
    if constexpr(Count){++st.exception_new_pairs;st.max_exception_count=std::max(st.max_exception_count,uint64_t(1));}
  }
  template<bool Count>uint32_t get(int c,ColoringStats&st)const {
    if constexpr(Count)++st.exc;
    for(auto p:x){if constexpr(Count)++st.exception_comparisons;if(p.first==c)return p.second;}
    return 0;
  }
};


template<bool Count=false,class Counter=uint32_t>
VI coloring_mags_block(const Summary&s,ColoringStats&st) {
  require(s.packed,"mags block coloring packed");
  std::vector<std::vector<Counter>>bc(s.b);
  std::vector<ColorExc<Counter>>exc(s.n);VI c(s.n,-1),plus_tag(128,-1);
  if constexpr(Count)st.state_bytes=uint64_t(s.b)*sizeof(bc[0])+uint64_t(s.n)*sizeof(exc[0]);
  for(int orig=0;orig<s.n;++orig) {
    int u=s.query_id(orig),b=s.owner[u];
    if constexpr(Count)++st.vertices;
    for(int j=s.plus.off[u];j<s.plus.off[u+1];++j) {
      int v=s.plus.to[j];if constexpr(Count)++st.plus;
      if(s.input_id(v)<orig){
        int k=c[v];require(k>=0,"Mags block coloring uncolored predecessor");
        if(size_t(k)>=plus_tag.size())plus_tag.resize(size_t(k)+1,-1);
        plus_tag[k]=orig;if constexpr(Count)++st.earlier;
      }
    }
    int col=0;
    for(;;++col) {
      require(col<s.n,"color index must be smaller than number of vertices");
      if constexpr(Count)++st.color_tests;
      if(size_t(col)<plus_tag.size()&&plus_tag[col]==orig)continue;
      uint32_t count=(size_t(col)<bc[b].size()?bc[b][col]:Counter(0));
      if(count<=exc[u].template get<Count>(col,st))break;
    }
    c[u]=col;
    for(int j=s.super.off[b];j<s.super.off[b+1];++j) {
      auto &row=bc[s.super.to[j]];
      if(size_t(col)>=row.size())row.resize(size_t(col)+1,Counter(0));
      require(row[col]<std::numeric_limits<Counter>::max(),"block counter overflow");
      auto &x=row[col];++x;
      if constexpr(Count){++st.se;++st.block_updates;st.max_block_count=std::max(st.max_block_count,uint64_t(x));}
    }
    for(int j=s.minus.off[u];j<s.minus.off[u+1];++j) {
      int v=s.minus.to[j];if constexpr(Count)++st.minus;
      if(s.input_id(v)>orig)exc[v].template add<Count>(col,st);
    }
  }
  VI out(s.n);for(int u=0;u<s.n;++u)out[s.input_id(u)]=c[u];return out;
}

}
