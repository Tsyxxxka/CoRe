#include "runner.hpp"
#include "kernels/bfs.hpp"
int main(int argc,char**argv){try{
 bench::Runner r("BFS",argc,argv,true);
 using namespace dq3;
 if(r.mode=="CSR")r.add("csr",[&]{return bench::Result{bfs_csr(*r.graph,r.source)};});
 else if(r.mode=="CoRe"){
   r.add("expand",[&]{return bench::Result{flood_expand<false>(*r.packed,r.source)};});
   r.add("sparse_block",[&]{return bench::Result{dqnext::bfs_core_sparse(*r.packed,r.source)};});
 }else{
   r.add("expand",[&]{return bench::Result{flood_expand<true>(*r.packed,r.source)};});
   r.add("runtime_hybrid",[&]{return bench::Result{dq15::bfs_mags_runtime_hybrid(*r.packed,r.source)};});
 }
 return r.execute();
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
