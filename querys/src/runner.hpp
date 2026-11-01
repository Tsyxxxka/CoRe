#pragma once
// Benchmark infrastructure only. Kernels are compiled in the four dedicated
// translation units. Loading, verification and timing are never inside a kernel.
#include "kernels/validate.hpp"
#include <variant>
#include <functional>
#include <memory>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <climits>
#include <random>
#include <cstdint>
#include <cerrno>

namespace bench {
using namespace dq3;
using Result=std::variant<VI,VD,DFS>;
struct Job {std::string kernel;std::function<Result()> run;};
inline uint64_t hash_result(const Result&r){
 if(auto x=std::get_if<VI>(&r))return checksum(*x);
 if(auto x=std::get_if<DFS>(&r))return checksum(x->parent)^checksum(x->tin)^checksum(x->tout);
 uint64_t h=0;for(double v:std::get<VD>(r)){uint64_t z;std::memcpy(&z,&v,8);h=(h^z)*0x9e3779b185ebca87ULL;}return h;
}
inline std::pair<int,i64> graph_header(const std::string&p){
 std::ifstream f(p,std::ios::binary);i64 n=-1,m=-1;f.read((char*)&n,8);f.read((char*)&m,8);
 require(bool(f)&&n>0&&n<=INT32_MAX/2&&m>=0&&m<=INT32_MAX/2,"invalid graph header");return {int(n),m};
}
inline int maximum_degree_source(const std::string&p,int n,i64 m){
 std::ifstream f(p,std::ios::binary);f.seekg(16);VI degrees(n,0);
 constexpr size_t B=1<<16;std::vector<Edge>buf(B);
 for(i64 done=0;done<m;){size_t k=size_t(std::min<i64>(B,m-done));f.read((char*)buf.data(),k*sizeof(Edge));require(bool(f),"incomplete graph edge list");
 for(size_t j=0;j<k;++j){auto [u,v]=buf[j];require(u>=0&&u<v&&v<n,"invalid graph edge");++degrees[u];++degrees[v];}done+=k;}
 int source=0;for(int u=1;u<n;++u)if(degrees[u]>degrees[source])source=u;return source;
}
// Stable, per-dataset seeded random source: the basename, graph dimensions,
// and QUERY_SEED give the same source for CSR, Mags and CoRe.
inline uint64_t query_seed(){
 const char* env=std::getenv("QUERY_SEED");
 if(!env||!*env)return 20261101ULL;
 require(env[0]!='-',"QUERY_SEED must be an unsigned 64-bit integer");
 errno=0;char*end=nullptr;
 unsigned long long seed=std::strtoull(env,&end,10);
 require(errno!=ERANGE&&end!=env&&*end=='\0',"invalid QUERY_SEED");
 return uint64_t(seed);
}
inline int random_source(const std::string&path,int n,i64 m,uint64_t seed){
 const std::string filename=std::filesystem::path(path).filename().string();
 uint64_t h=14695981039346656037ULL;
 for(unsigned char c:filename){h^=c;h*=1099511628211ULL;}
 auto mix=[&](uint64_t x){for(int i=0;i<8;++i){h^=(x>>(i*8))&255ULL;h*=1099511628211ULL;}};
 mix(uint64_t(n));mix(uint64_t(m));mix(seed);
 std::mt19937_64 rng(h);
 return std::uniform_int_distribution<int>(0,n-1)(rng);
}
inline void same_result(const Result&a,const Result&b){
 if(a.index()!=b.index())throw std::runtime_error("result variant changed");
 if(auto x=std::get_if<DFS>(&a)){const auto&y=std::get<DFS>(b);require(x->parent==y.parent&&x->tin==y.tin&&x->tout==y.tout,"non-deterministic DFS result");}
 else if(auto x=std::get_if<VD>(&a)){const auto&y=std::get<VD>(b);require(x->size()==y.size(),"PR size mismatch");for(size_t i=0;i<x->size();++i)require(std::abs((*x)[i]-y[i])<1e-12,"non-deterministic PageRank result");}
 else require(std::get<VI>(a)==std::get<VI>(b),"non-deterministic integer result");
}
inline std::string json_escape(const std::string&s){std::string a;for(char c:s){if(c=='"'||c=='\\')a+='\\';if(c=='\n')a+="\\n";else a+=c;}return a;}
struct Runner {
 std::string task,mode,graph_path,prefix,output;
 int repeats=3,warmups=1,source=-1,n=0; i64 m=0;
 uint64_t source_seed=0; std::string source_method="legacy_max_degree";
 double metadata_seconds=0,representation_load_seconds=0,layout_seconds=0,verification_seconds=0;
 std::unique_ptr<Graph> graph;
 std::unique_ptr<Summary> raw,packed;
 std::vector<Job> jobs;
 Runner(const std::string&task_name,int argc,char**argv,bool need_packed):task(task_name){
   if(argc<5||argc>8)throw std::runtime_error("usage: QUERY MODE GRAPH.bin SUMMARY_PREFIX_OR_DIR OUT.jsonl [repeats=3] [warmup=1] [source=-1 (seeded random for BFS/DFS)]");
   mode=argv[1];graph_path=argv[2];prefix=argv[3];output=argv[4];
   repeats=argc>5?std::stoi(argv[5]):3;warmups=argc>6?std::stoi(argv[6]):1;
   int supplied=argc>7?std::stoi(argv[7]):-1;
   require(repeats>=1&&repeats<=100&&warmups>=0&&warmups<=20,"invalid repeat/warmup");
   auto t=Clock::now();std::tie(n,m)=graph_header(graph_path);
   if(supplied>=0){source=supplied;source_method="explicit";}
   else if(task=="BFS"||task=="DFS"){
     source_seed=query_seed();
     source=random_source(graph_path,n,m,source_seed);
     source_method="seeded_random";
   }else source=maximum_degree_source(graph_path,n,m);
   require(source>=0&&source<n,"source outside graph");metadata_seconds=elapsed(t);
   t=Clock::now();
   if(mode=="CSR")graph=std::make_unique<Graph>(graph_path);
   else if(mode=="CoRe"||mode=="Mags")raw=std::make_unique<Summary>(n,m,prefix);
   else throw std::runtime_error("MODE must be CSR, CoRe, or Mags");
   representation_load_seconds=elapsed(t);
   t=Clock::now();
   if(raw){if(mode=="CoRe")require(raw->nq==0,"CoRe summary contains C-");if(need_packed)packed=std::make_unique<Summary>(raw->pack(true));}
   layout_seconds=elapsed(t);
 }
 void add(std::string name,std::function<Result()>f){jobs.push_back({std::move(name),std::move(f)});}
 void verify(){
   if(!raw||(n>1000&&!(std::getenv("VERIFY")&&std::string(std::getenv("VERIFY"))=="1")))return;
   auto t=Clock::now();Graph original(graph_path);verify_summary(original,*raw);
   if(packed){verify_summary(original,*packed);verify_same_encoding(*raw,*packed);}
   verification_seconds=elapsed(t);
 }
 // Each algorithm is its own independently compiled executable; result JSONL
 // keys, timing boundaries and checksum logic follow benchmark_four.cpp.
 int execute(){
   verify();auto out_dir=std::filesystem::path(output).parent_path();if(!out_dir.empty())std::filesystem::create_directories(out_dir);
   std::ofstream os(output);require(bool(os),"cannot create output");os<<std::setprecision(17);
   for(const Job&job:jobs){
     try{
       Result reference=job.run();for(int j=0;j<warmups;++j)same_result(reference,job.run());
       std::vector<double>durations;durations.reserve(repeats);
       uint64_t checksum_value=hash_result(reference);
       for(int j=0;j<repeats;++j){auto start=Clock::now();Result r=job.run();double secs=elapsed(start);same_result(reference,r);durations.push_back(secs);}
       auto ordered=durations;std::sort(ordered.begin(),ordered.end());double median=(ordered[(repeats-1)/2]+ordered[repeats/2])/2;
       double sum=std::accumulate(durations.begin(),durations.end(),0.0)/repeats;
       os<<"{\"status\":\"ok\",\"representation\":\""<<mode<<"\",\"task\":\""<<task<<"\",\"kernel\":\""<<job.kernel<<"\",\"n\":"<<n<<",\"m\":"<<m<<",\"source\":"<<source<<",\"source_method\":\""<<source_method<<"\",\"source_seed\":"<<source_seed<<",\"pr_iters\":20,\"warmups\":"<<warmups<<",\"repeats\":"<<repeats<<",\"read_seconds\":"<<metadata_seconds+representation_load_seconds<<",\"graph_metadata_seconds\":"<<metadata_seconds<<",\"representation_load_seconds\":"<<representation_load_seconds<<",\"layout_seconds\":"<<layout_seconds<<",\"verification_seconds\":"<<verification_seconds<<",\"query_seconds_mean\":"<<sum<<",\"query_seconds_median\":"<<median<<",\"query_seconds_runs\":[";
       for(int j=0;j<repeats;++j)os<<(j?",":"")<<durations[j];
       os<<"],\"checksum\":"<<checksum_value<<"}\n";os.flush();
       std::cerr<<mode<<' '<<task<<' '<<job.kernel<<" median="<<median<<" sec\n";
     }catch(const std::exception&e){
       os<<"{\"status\":\"failed\",\"representation\":\""<<mode<<"\",\"task\":\""<<task<<"\",\"kernel\":\""<<job.kernel<<"\",\"message\":\""<<json_escape(e.what())<<"\"}\n";
       std::cerr<<"ERROR "<<task<<'/'<<job.kernel<<": "<<e.what()<<'\n';
     }
   }
   return 0;
 }
};
} // namespace bench
