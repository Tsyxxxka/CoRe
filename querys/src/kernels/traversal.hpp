#pragma once
#include "../upstream/storage.hpp"
namespace dq3 {
struct DFS {VI parent,tin,tout;};


inline DFS restore_dfs(const Summary&s,DFS r){if(!s.packed)return r;DFS o;o.parent.resize(s.n);o.tin.resize(s.n);o.tout.resize(s.n);for(int u=0;u<s.n;u++){int x=s.orig[u];o.parent[x]=r.parent[u]<0?r.parent[u]:s.orig[r.parent[u]];o.tin[x]=r.tin[u];o.tout[x]=r.tout[u];}return o;}


inline VI bfs_csr(const Graph&g,int src){VI d(g.n,-1),q(g.n);int h=0,t=0;d[src]=0;q[t++]=src;while(h<t){int u=q[h++];for(int j=g.adj.off[u];j<g.adj.off[u+1];j++){int v=g.adj.to[j];if(d[v]<0){d[v]=d[u]+1;q[t++]=v;}}}return d;}


struct NextLive {VI next;explicit NextLive(int n):next(n+1){std::iota(next.begin(),next.end(),0);}int find(int x){int r=x;while(next[r]!=r)r=next[r];while(next[x]!=x){int z=next[x];next[x]=r;x=z;}return r;}void erase(int x){next[x]=find(x+1);}};


template<bool Signed,bool WCC=false>VI flood_expand(const Summary&s,int src=0){VI d(s.n,-1),q(s.n),mark;if constexpr(Signed)mark.assign(s.n,-1);int h=0,t=0;for(int it=0;it<(WCC?s.n:1);it++){int root=WCC?s.query_id(it):s.query_id(src);if(d[root]>=0)continue;d[root]=WCC?it:0;q[t++]=root;while(h<t){int u=q[h++],val=WCC?d[u]:d[u]+1;auto visit=[&](int v){if(d[v]<0){d[v]=val;q[t++]=v;}};for(int j=s.plus.off[u];j<s.plus.off[u+1];j++)visit(s.plus.to[j]);if constexpr(Signed)for(int j=s.minus.off[u];j<s.minus.off[u+1];j++)mark[s.minus.to[j]]=u;int a=s.owner[u];for(int j=s.super.off[a];j<s.super.off[a+1];j++){int b=s.super.to[j];for(int k=s.moff[b];k<s.moff[b+1];k++){int v=s.vertex_at(k);if constexpr(Signed){if(mark[v]==u)continue;}visit(v);}}}}return s.restore(std::move(d));}


template<bool Signed>DFS dfs_expand(const Summary&s,int src=0){DFS r;r.parent.assign(s.n,-2);r.tin.resize(s.n);r.tout.resize(s.n);struct F{int u,cp,e,pos;};std::vector<F>st;st.reserve(s.n);int tm=0;auto push=[&](int u,int par){r.parent[u]=par;r.tin[u]=tm++;st.push_back({u,s.plus.off[u],s.super.off[s.owner[u]],-1});};for(int id=src,seen=0;seen<s.n;++seen){int root=s.query_id(id);if(++id==s.n)id=0;if(r.parent[root]!=-2)continue;push(root,-1);while(!st.empty()){auto&f=st.back();int u=f.u,v=-1,a=s.owner[u];while(f.cp<s.plus.off[u+1]){int x=s.plus.to[f.cp++];if(r.parent[x]==-2){v=x;break;}}while(v<0&&f.e<s.super.off[a+1]){int b=s.super.to[f.e];if(f.pos<0)f.pos=s.moff[b];while(f.pos<s.moff[b+1]){int x=s.vertex_at(f.pos++);if(r.parent[x]!=-2)continue;if constexpr(Signed){if(s.negative(u,x))continue;}v=x;break;}if(v<0){++f.e;f.pos=-1;}}if(v>=0)push(v,u);else{r.tout[u]=tm++;st.pop_back();}}}return restore_dfs(s,std::move(r));}

}
