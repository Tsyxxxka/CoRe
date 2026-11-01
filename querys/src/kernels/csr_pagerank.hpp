#pragma once
#include "../upstream/storage.hpp"
namespace dq3 {
inline double gather4(const double* z,const int* to,int lo,int hi){double a=0,b=0,c=0,d=0;int j=lo;for(;j+4<=hi;j+=4){a+=z[to[j]];b+=z[to[j+1]];c+=z[to[j+2]];d+=z[to[j+3]];}double x=((a+b)+c)+d;for(;j<hi;j++)x+=z[to[j]];return x;}


inline VD pr_csr(const Graph&g,int T=20,double alpha=.85,bool fused=true){
 require(T>=0&&alpha>=0&&alpha<1,"PR parameters");const int n=g.n;VD inv(n),z(n),next(n),out(n);int n0=0;
 for(int u=0;u<n;u++){int d=g.adj.deg(u);if(d)inv[u]=alpha/d;else ++n0;z[u]=inv[u]/n;}
 if(T==0)return VD(n,1./n);
 double sink=double(n0)/n;
 if(!fused){VD r(n,1./n);for(int t=0;t<T;t++){for(int u=0;u<n;u++)z[u]=inv[u]*r[u];double tele=(1-alpha+alpha*sink)/n;for(int u=0;u<n;u++)out[u]=tele+gather4(z.data(),g.adj.to.data(),g.adj.off[u],g.adj.off[u+1]);r.swap(out);sink=n0*tele;}return r;}
 for(int t=0;t<T;t++){double tele=(1-alpha+alpha*sink)/n;bool last=t+1==T;for(int u=0;u<n;u++){double x=tele+gather4(z.data(),g.adj.to.data(),g.adj.off[u],g.adj.off[u+1]);if(last)out[u]=x;else next[u]=inv[u]*x;}z.swap(next);sink=n0*tele;}return out;
}

}
