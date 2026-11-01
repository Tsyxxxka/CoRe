#include "runner.hpp"
#include "kernels/pagerank.hpp"
#include "kernels/csr_pagerank.hpp"
int main(int argc,char**argv){try{
 bench::Runner r("PR",argc,argv,false);
 using namespace dq3;
 if(r.mode=="CSR")r.add("two_pass",[&]{return bench::Result{pr_csr(*r.graph,20,.85,false)};});
 else{
   auto t=Clock::now();auto half=std::make_shared<dq5::Summary>(*r.raw);
   r.layout_seconds+=elapsed(t);
   if(r.mode=="CoRe")r.add("dq5_auto",[half]{return bench::Result{dq5::pagerank<false>(*half,20,.85,dq5::Backend::Automatic)};});
   else r.add("dq5_auto",[half]{return bench::Result{dq5::pagerank<true>(*half,20,.85,dq5::Backend::Automatic)};});
 }
 return r.execute();
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
