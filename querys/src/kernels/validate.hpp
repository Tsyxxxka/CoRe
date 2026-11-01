#pragma once
#include "traversal.hpp"
#include <cmath>
namespace dq3 {
inline void verify_summary(const Graph&g,const Summary&s){for(int id=0;id<g.n;id++){int u=s.query_id(id);VI z;s.neighbors(u,[&](int v){z.push_back(s.input_id(v));});std::sort(z.begin(),z.end());if(z.size()!=size_t(g.adj.deg(id))||!std::equal(z.begin(),z.end(),g.adj.to.begin()+g.adj.off[id]))throw std::runtime_error("summary recovery at "+std::to_string(id));}}


inline void verify_same_encoding(const Summary&raw,const Summary&packed){
 require(!raw.packed&&packed.packed,"encoding comparison layouts");
 require(raw.n==packed.n&&raw.b==packed.b&&raw.m==packed.m,"encoding dimensions");
 VI own(raw.n);for(int u=0;u<packed.n;u++)own[packed.orig[u]]=packed.owner[u];require(own==raw.owner,"partition changed");
 auto edges=[&](const CSR&z,bool blocks){Edges out;int n=blocks?packed.b:packed.n;for(int u=0;u<n;u++)for(int j=z.off[u];j<z.off[u+1];j++){int v=z.to[j];if(u>v)continue;if(blocks)out.push_back({u,v});else out.push_back(std::minmax(packed.orig[u],packed.orig[v]));}std::sort(out.begin(),out.end());return out;};
 require(edges(packed.super,true)==raw.raw_se,"superedge encoding changed");require(edges(packed.plus,false)==raw.raw_cp,"positive correction changed");require(edges(packed.minus,false)==raw.raw_cm,"negative correction changed");
}


inline void verify_dfs(const Graph&g,const DFS&r){require(r.parent.size()==size_t(g.n)&&r.tin.size()==size_t(g.n)&&r.tout.size()==size_t(g.n),"DFS length");VI events(2*g.n,INT32_MIN),stack;auto anc=[&](int u,int v){return r.tin[u]<=r.tin[v]&&r.tout[u]>=r.tout[v];};for(int u=0;u<g.n;u++){require(r.tin[u]>=0&&r.tin[u]<r.tout[u]&&r.tout[u]<2*g.n,"DFS clock");require(events[r.tin[u]]==INT32_MIN&&events[r.tout[u]]==INT32_MIN,"DFS duplicate time");events[r.tin[u]]=u;events[r.tout[u]]=~u;int p=r.parent[u];require(p>=-1&&p<g.n,"DFS parent");if(p>=0)require(anc(p,u)&&std::binary_search(g.adj.to.begin()+g.adj.off[u],g.adj.to.begin()+g.adj.off[u+1],p),"DFS nonedge parent");}for(int x:events)if(x>=0){require(r.parent[x]==(stack.empty()?-1:stack.back()),"DFS immediate parent");stack.push_back(x);}else{require(!stack.empty()&&stack.back()==~x,"DFS stack nesting");stack.pop_back();}require(stack.empty(),"DFS open stack");for(auto[u,v]:g.edges)require(anc(u,v)||anc(v,u),"DFS cross edge");}


inline double error(const VD&a,const VD&b){require(a.size()==b.size(),"PR size");double er=0;for(size_t i=0;i<a.size();i++){require(std::isfinite(a[i]),"PR nonfinite");er=std::max(er,std::abs(a[i]-b[i]));}require(er<1e-11,"PR error");return er;}

}
