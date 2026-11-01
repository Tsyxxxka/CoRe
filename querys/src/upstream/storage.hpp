#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
namespace dq3 {
using VI=std::vector<int>;using VD=std::vector<double>;using Edge=std::pair<int,int>;using Edges=std::vector<Edge>;using i64=std::int64_t;
using Clock=std::chrono::steady_clock;
inline double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
inline void require(bool b,const char*why){if(!b)throw std::runtime_error(why);}
template<class T> std::vector<T> read_raw(const std::string&path){std::ifstream f(path,std::ios::binary|std::ios::ate);if(!f)throw std::runtime_error("cannot read "+path);auto sz=f.tellg();require(sz>=0&&size_t(sz)%sizeof(T)==0,"invalid raw length");std::vector<T>v(size_t(sz)/sizeof(T));f.seekg(0);f.read((char*)v.data(),sz);require(bool(f),"short read");return v;}
struct CSR {VI off,to;CSR()=default;CSR(int n,const Edges&e){off.assign(n+1,0);for(auto[u,v]:e){require(u>=0&&u<n&&v>=0&&v<n,"CSR endpoint");++off[u+1];if(u!=v)++off[v+1];}for(int u=0;u<n;u++){i64 x=(i64)off[u]+off[u+1];require(x<=INT32_MAX,"32-bit CSR overflow");off[u+1]=x;}to.resize(off.back());auto p=off;for(auto[u,v]:e){to[p[u]++]=v;if(u!=v)to[p[v]++]=u;}for(int u=0;u<n;u++)std::sort(to.begin()+off[u],to.begin()+off[u+1]);}int deg(int u)const{return off[u+1]-off[u];}size_t bytes()const{return 4*(off.size()+to.size());}};
struct Graph {int n=0;i64 m=0;CSR adj;Edges edges;Graph()=default;Graph(int n_,Edges e):n(n_),m(e.size()),adj(n,e),edges(std::move(e)){}explicit Graph(const std::string&p){std::ifstream f(p,std::ios::binary);i64 nn=-1;f.read((char*)&nn,8);f.read((char*)&m,8);require(bool(f)&&nn>0&&nn<=INT32_MAX/2&&m>=0&&m<=INT32_MAX/2,"graph header");n=nn;edges.resize(m);f.read((char*)edges.data(),m*8);require(bool(f),"graph short read");char c;require(!f.read(&c,1),"graph extra bytes");require(std::is_sorted(edges.begin(),edges.end())&&std::adjacent_find(edges.begin(),edges.end())==edges.end(),"graph not sorted unique");for(auto[u,v]:edges)require(u>=0&&u<v&&v<n,"graph not simple");adj=CSR(n,edges);}};
struct Summary {
 int n=0,b=0; i64 m=0,ns=0,np=0,nq=0,nloop=0;bool packed=false;
 VI owner, moff, member, pos, cut, orig, perm;CSR super,plus,minus;Edges raw_se,raw_cp,raw_cm;
 Summary()=default;
 Summary(int nn,i64 mm,const std::string&pre):n(nn),m(mm){owner=read_raw<int>(pre+".partition.bin");auto se=read_raw<Edge>(pre+".se.bin"),cp=read_raw<Edge>(pre+".cp.bin"),cm=read_raw<Edge>(pre+".cm.bin");init(se,cp,cm);}
 Summary(int nn,i64 mm,VI ow,const Edges&se,const Edges&cp,const Edges&cm):n(nn),m(mm),owner(std::move(ow)){init(se,cp,cm);}
 void init(const Edges&se,const Edges&cp,const Edges&cm){require((int)owner.size()==n&&n>0,"owner length");b=0;for(int a:owner){require(a>=0&&a<n,"owner invalid");b=std::max(b,a+1);}moff.assign(b+1,0);for(int a:owner)++moff[a+1];std::partial_sum(moff.begin(),moff.end(),moff.begin());member.resize(n);pos.resize(n);auto p=moff;for(int u=0;u<n;u++){int k=p[owner[u]]++;member[k]=u;pos[u]=k;}for(int a=0;a<b;a++)require(moff[a]<moff[a+1],"empty block");ns=se.size();np=cp.size();nq=cm.size();nloop=0;for(auto[a,c]:se)if(a==c)++nloop;
 auto ck=[](const Edges&e,bool loop,int lim){require(std::is_sorted(e.begin(),e.end())&&std::adjacent_find(e.begin(),e.end())==e.end(),"summary not sorted unique");for(auto[u,v]:e)require(u>=0&&v<lim&&u<=v&&(loop||u!=v),"summary endpoint invalid");};ck(se,true,b);ck(cp,false,n);ck(cm,false,n);raw_se=se;raw_cp=cp;raw_cm=cm;super=CSR(b,se);plus=CSR(n,cp);minus=CSR(n,cm);
 }
 int size(int a)const{return moff[a+1]-moff[a];}
 int vertex_at(int i)const{return packed?i:member[i];}
 bool negative(int u,int v)const{return std::binary_search(minus.to.begin()+minus.off[u],minus.to.begin()+minus.off[u+1],v);}
 template<class F>void neighbors(int u,F f)const {for(int j=plus.off[u];j<plus.off[u+1];j++)f(plus.to[j]);int a=owner[u];for(int j=super.off[a];j<super.off[a+1];j++){int c=super.to[j];for(int k=moff[c];k<moff[c+1];k++){int v=vertex_at(k);if(v!=u&&!negative(u,v))f(v);}}}
 template<class T>std::vector<T> restore(const std::vector<T>&x)const{if(!packed)return x;std::vector<T>out(n);for(int i=0;i<n;i++)out[orig[i]]=x[i];return out;}
 template<class T>std::vector<T> restore(std::vector<T>&&x)const{if(!packed)return std::move(x);std::vector<T>out(n);for(int i=0;i<n;i++)out[orig[i]]=x[i];return out;}
 int input_id(int u)const{return packed?orig[u]:u;}
 int query_id(int u)const{return packed?perm[u]:u;}
 i64 records()const{return ns+np+nq;}i64 native()const{return 2*(ns+np+nq)-nloop;}
 size_t bytes()const{return 4*(owner.size()+moff.size()+member.size()+pos.size()+cut.size()+orig.size()+perm.size())+super.bytes()+plus.bytes()+minus.bytes()+8*(raw_se.size()+raw_cp.size()+raw_cm.size());}
 Summary pack(bool empty_suffix=true)const{
  require(!packed,"already packed");Summary t;t.n=n;t.m=m;t.b=b;t.ns=ns;t.np=np;t.nq=nq;t.nloop=nloop;t.packed=true;t.super=super;t.moff=moff;t.cut.resize(b);t.owner.resize(n);t.orig.reserve(n);t.perm.resize(n);
  for(int a=0;a<b;a++){if(empty_suffix){for(int i=moff[a];i<moff[a+1];i++){int u=member[i];if(plus.deg(u)||minus.deg(u))t.orig.push_back(u);}t.cut[a]=t.orig.size();for(int i=moff[a];i<moff[a+1];i++){int u=member[i];if(!plus.deg(u)&&!minus.deg(u))t.orig.push_back(u);}}else{for(int i=moff[a];i<moff[a+1];i++)t.orig.push_back(member[i]);t.cut[a]=moff[a+1];}}
  for(int i=0;i<n;i++){t.perm[t.orig[i]]=i;t.owner[i]=owner[t.orig[i]];}
  auto cv=[&](const CSR&c){CSR z;z.off.resize(n+1);z.off[0]=0;for(int u=0;u<n;u++)z.off[u+1]=z.off[u]+c.deg(t.orig[u]);z.to.resize(c.to.size());for(int u=0;u<n;u++){int old=t.orig[u],w=z.off[u];for(int j=c.off[old];j<c.off[old+1];j++)z.to[w++]=t.perm[c.to[j]];std::sort(z.to.begin()+z.off[u],z.to.begin()+z.off[u+1]);}return z;};t.plus=cv(plus);t.minus=cv(minus);
  // These are alternative PHYSICAL views of the original undirected records,
  // not graph-derived query certificates. Their bytes are explicitly counted.
  t.raw_se=raw_se;
  auto single=[&](const CSR&z){Edges out;out.reserve(z.to.size()/2);for(int u=0;u<n;u++)for(int j=z.off[u];j<z.off[u+1];j++)if(u<z.to[j])out.push_back({u,z.to[j]});return out;};
  t.raw_cp=single(t.plus);t.raw_cm=single(t.minus);return t;
 }
 // Versioned, lossless native storage; no derived degrees or certificates.
 void save(const std::string&path)const{require(packed,"save requires packed storage");std::ofstream f(path,std::ios::binary);require(bool(f),"save open");const char magic[8]={'D','Q','3','S','B','R','0','2'};f.write(magic,8);i64 h[7]={n,b,m,ns,np,nq,nloop};f.write((char*)h,sizeof(h));auto w=[&](const VI&v){uint64_t k=v.size();f.write((char*)&k,8);f.write((char*)v.data(),k*4);};w(owner);w(moff);w(cut);w(orig);w(super.off);w(super.to);w(plus.off);w(plus.to);w(minus.off);w(minus.to);auto we=[&](const Edges&e){uint64_t k=e.size();f.write((char*)&k,8);f.write((char*)e.data(),k*8);};we(raw_se);we(raw_cp);we(raw_cm);f.close();require(bool(f),"save failed");}
 static Summary load(const std::string&path){std::ifstream f(path,std::ios::binary);require(bool(f),"load native open");char magic[8];f.read(magic,8);require(!std::memcmp(magic,"DQ3SBR02",8),"bad native magic");i64 h[7];f.read((char*)h,sizeof(h));require(bool(f)&&h[0]>0&&h[0]<=INT32_MAX/2&&h[1]>0&&h[1]<=h[0]&&h[2]>=0&&h[2]<=INT32_MAX/2&&h[3]>=0&&h[3]<=INT32_MAX/2&&h[4]>=0&&h[4]<=INT32_MAX/2&&h[5]>=0&&h[5]<=INT32_MAX/2&&h[6]>=0&&h[6]<=h[3],"bad native header");Summary s;s.n=h[0];s.b=h[1];s.m=h[2];s.ns=h[3];s.np=h[4];s.nq=h[5];s.nloop=h[6];s.packed=true;auto r=[&](VI&v){uint64_t k;f.read((char*)&k,8);require(bool(f)&&k<=INT32_MAX,"native vector length");v.resize(k);f.read((char*)v.data(),k*4);require(bool(f),"native short read");};r(s.owner);r(s.moff);r(s.cut);r(s.orig);r(s.super.off);r(s.super.to);r(s.plus.off);r(s.plus.to);r(s.minus.off);r(s.minus.to);auto re=[&](Edges&e){uint64_t k;f.read((char*)&k,8);require(bool(f)&&k<=INT32_MAX/2,"native record length");e.resize(k);f.read((char*)e.data(),k*8);require(bool(f),"native record short read");};re(s.raw_se);re(s.raw_cp);re(s.raw_cm);char c;require(!f.read(&c,1),"native extra bytes");require(s.owner.size()==size_t(s.n)&&s.orig.size()==size_t(s.n)&&s.moff.size()==size_t(s.b+1)&&s.cut.size()==size_t(s.b),"native array dimensions");s.perm.assign(s.n,-1);for(int i=0;i<s.n;i++){require(s.orig[i]>=0&&s.orig[i]<s.n&&s.perm[s.orig[i]]==-1,"bad native permutation");s.perm[s.orig[i]]=i;}auto ck=[&](const CSR&z,int rows,int lim){require(z.off.size()==size_t(rows+1)&&z.off[0]==0&&z.off.back()==(int)z.to.size()&&std::is_sorted(z.off.begin(),z.off.end()),"bad native CSR");for(int x:z.to)require(x>=0&&x<lim,"bad native endpoint");for(int i=0;i<rows;i++){auto l=z.to.begin()+z.off[i],r=z.to.begin()+z.off[i+1];require(std::is_sorted(l,r)&&std::adjacent_find(l,r)==r,"native row not sorted unique");}};ck(s.super,s.b,s.b);ck(s.plus,s.n,s.n);ck(s.minus,s.n,s.n);require(s.moff[0]==0&&s.moff.back()==s.n,"bad native member range");for(int a=0;a<s.b;a++){require(s.moff[a]<s.moff[a+1]&&s.cut[a]>=s.moff[a]&&s.cut[a]<=s.moff[a+1],"native row split");for(int u=s.moff[a];u<s.moff[a+1];u++){require(s.owner[u]==a,"native ownership");if(u>=s.cut[a])require(s.plus.deg(u)==0&&s.minus.deg(u)==0,"nonempty suffix");}}require((i64)s.super.to.size()==2*s.ns-s.nloop&&(i64)s.plus.to.size()==2*s.np&&(i64)s.minus.to.size()==2*s.nq,"native record counts");require((i64)s.raw_se.size()==s.ns&&(i64)s.raw_cp.size()==s.np&&(i64)s.raw_cm.size()==s.nq,"native single-copy record counts");return s;}
};
struct Degrees{VI degree,base;std::vector<unsigned char>loop;explicit Degrees(const Summary&s):degree(s.n),base(s.b),loop(s.b,0){for(int a=0;a<s.b;a++){i64 x=0;for(int j=s.super.off[a];j<s.super.off[a+1];j++){int c=s.super.to[j];x+=s.size(c);if(a==c)loop[a]=1;}x-=loop[a];require(x<=INT32_MAX,"degree base overflow");base[a]=x;}for(int u=0;u<s.n;u++){i64 d=(i64)base[s.owner[u]]+s.plus.deg(u)-s.minus.deg(u);require(d>=0&&d<s.n,"invalid degree");degree[u]=d;}}};
inline uint64_t checksum(const VI&a){uint64_t h=0;for(int x:a)h=(h^uint64_t(x+2))*0x9e3779b185ebca87ULL;return h;}
}
