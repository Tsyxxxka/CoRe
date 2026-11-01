#pragma once
#include "traversal.hpp"
namespace dqnext {
using namespace dq3;
struct OpCounts {
 uint64_t vertex_disc=0, queue_pops=0, plus_probes=0, block_fires=0,
          super_probes=0, block_member_tests=0, block_bulk_disc=0;
};


template<bool Count=false>
VI bfs_core_sparse(const Summary&s,int src,OpCounts* c=nullptr){
 require(s.packed&&s.nq==0,"bfs_core_sparse requires packed positive summary");
 VI d(s.n,-1),q(s.n); std::vector<unsigned char> fired(s.b,0),pending(s.b,0),filled(s.b,0); int h=0,t=0;
 auto discover=[&](int u,int val,bool bulk){
  if(d[u]>=0)return;
  d[u]=val;if constexpr(Count){++c->vertex_disc;if(bulk)++c->block_bulk_disc;}
  int a=s.owner[u]; bool residual=s.plus.off[u]!=s.plus.off[u+1];
  if(residual || (!fired[a]&&!pending[a])){q[t++]=u;if(!fired[a])pending[a]=1;}
 };
 discover(s.query_id(src),0,false);
 while(h<t){
  int u=q[h++],a=s.owner[u],val=d[u]+1;if constexpr(Count)++c->queue_pops;
  for(int j=s.plus.off[u];j<s.plus.off[u+1];++j){if constexpr(Count)++c->plus_probes;discover(s.plus.to[j],val,false);}
  if(fired[a])continue; fired[a]=1;pending[a]=0;if constexpr(Count)++c->block_fires;
  for(int j=s.super.off[a];j<s.super.off[a+1];++j){
   int b=s.super.to[j];if constexpr(Count)++c->super_probes;if(filled[b])continue;filled[b]=1;
   // packed members are contiguous; each block is bulk-scanned at most once.
   for(int v=s.moff[b];v<s.moff[b+1];++v){if constexpr(Count)++c->block_member_tests;discover(v,val,true);}
  }
 }
 return s.restore(std::move(d));
}

}

namespace dq15 {
using namespace dq3;
struct BFSHybridCounts { uint64_t queue_pops=0,plus=0,minus=0,super=0,fast=0,slow=0,member=0,disc=0; };
template<bool Count=false>
VI bfs_mags_runtime_hybrid(const Summary& s,int src,BFSHybridCounts* c=nullptr){
  require(s.packed,"runtime hybrid requires packed summary");
  VI d(s.n,-1),q(s.n),markV(s.n,-1),markB(s.b,-1); NextLive live(s.n),rels(s.super.to.size());
  int h=0,t=0; auto visit=[&](int v,int val){if(d[v]>=0)return;d[v]=val;q[t++]=v;live.erase(v);if constexpr(Count)++c->disc;};
  visit(s.query_id(src),0);
  while(h<t){int u=q[h++],a=s.owner[u],val=d[u]+1;if constexpr(Count)++c->queue_pops;
    for(int j=s.minus.off[u];j<s.minus.off[u+1];++j){int x=s.minus.to[j];markV[x]=u;markB[s.owner[x]]=u;if constexpr(Count)++c->minus;}
    for(int j=s.plus.off[u];j<s.plus.off[u+1];++j){if constexpr(Count)++c->plus;int v=s.plus.to[j];if(d[v]<0)visit(v,val);}
    int e=rels.find(s.super.off[a]);
    while(e<s.super.off[a+1]){int b=s.super.to[e];if constexpr(Count)++c->super;int j=live.find(s.moff[b]);
      if(markB[b]!=u){if constexpr(Count)++c->fast;while(j<s.moff[b+1]){int v=j;if constexpr(Count)++c->member;visit(v,val);j=live.find(j+1);}}
      else {if constexpr(Count)++c->slow;while(j<s.moff[b+1]){int v=j;if constexpr(Count)++c->member;if(markV[v]!=u)visit(v,val);j=live.find(j+1);}}
      if(live.find(s.moff[b])>=s.moff[b+1])rels.erase(e); e=rels.find(e+1);
    }
  }
  return s.restore(std::move(d));
}
}
