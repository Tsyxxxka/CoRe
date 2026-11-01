#include "runner.hpp"
#include "kernels/dfs.hpp"
int main(int argc,char**argv){try{
 bench::Runner r("DFS",argc,argv,true);
 using namespace dq3;
 if(r.mode=="CSR")r.add("parent_spine",[&]{return bench::Result{dq4::dfs_csr_spine(*r.graph,nullptr,r.source)};});
 else if(r.mode=="CoRe"){
   r.add("expand",[&]{return bench::Result{dfs_expand<false>(*r.packed,r.source)};});
   r.add("spine",[&]{return bench::Result{dq4::dfs_spine(*r.packed,nullptr,r.source)};});
 }else r.add("expand_frame",[&]{return bench::Result{dfs_expand<true>(*r.packed,r.source)};});
 return r.execute();
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
