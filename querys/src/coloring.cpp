#include "runner.hpp"
#include "kernels/coloring.hpp"
int main(int argc,char**argv){try{
 bench::Runner r("Coloring",argc,argv,true);
 using namespace dq3;
 if(r.mode=="CSR")r.add("first_fit",[&]{dq15::ColoringStats s;return bench::Result{dq15::coloring_csr(*r.graph,s)};});
 else if(r.mode=="CoRe"){
   r.add("expand",[&]{dq15::ColoringStats s;return bench::Result{dq15::coloring_expand<false,false>(*r.packed,s)};});
   r.add("block_bit",[&]{dq15::ColoringStats s;return bench::Result{dq15::coloring_core_block<false>(*r.packed,s)};});
 }else{
   r.add("expand",[&]{dq15::ColoringStats s;return bench::Result{dq15::coloring_expand<true,false>(*r.packed,s)};});
   r.add("block_u32",[&]{dq15::ColoringStats s;return bench::Result{dq15::coloring_mags_block<false,uint32_t>(*r.packed,s)};});
 }
 return r.execute();
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
