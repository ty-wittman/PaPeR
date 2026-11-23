#ifndef PAP_FOUR_GENE_IDENTITY_H
#define PAP_FOUR_GENE_IDENTITY_H
// Exact structural gene-origin enumeration. Admission belongs to the caller.
#include "relationship_topology.h"
#include <functional>
#include <cstdint>
namespace pap_internal::four_gene {
enum Reason { endpoint_count=1, latent_observation=2, ancestral_exterior=4,
 joint_exterior=8, invalid_dag=16, nonancestral=32, factor_ownership=64,
 not_ready=128, noncanonical=256 };
struct Certificate {int mask=0;std::vector<int> axes,order;};
inline Certificate certify(const ::relationship_topology::TopologyModel& core,
 const std::vector<unsigned char>& typed,const std::vector<unsigned char>& observed,
 const std::vector<unsigned char>& exterior,bool joint,bool owned,size_t endpoints=2){
 Certificate out;int n=core.nodes.size();if(!core.ready())out.mask|=not_ready;
 if(!owned)out.mask|=factor_ownership;if(joint)out.mask|=joint_exterior;
 if(typed.size()!=size_t(n)||observed.size()!=size_t(n)||exterior.size()!=size_t(n)){out.mask|=factor_ownership;return out;}
 for(int v=0;v<n;++v){if(typed[v])out.axes.push_back(v);else {if(observed[v])out.mask|=latent_observation;if(exterior[v])out.mask|=ancestral_exterior;}}
 if(out.axes.size()!=endpoints)out.mask|=endpoint_count;
 for(int v=0;v<n;++v)if(core.nodes[v].has_birth)for(int p:core.nodes[v].parents)if(p<0||p>=n)out.mask|=invalid_dag;
 if(out.mask&invalid_dag)return out;
 std::vector<unsigned char> done(n);while(out.order.size()<size_t(n)){bool advanced=false;for(int v=0;v<n;++v)if(!done[v]&&(!core.nodes[v].has_birth||(done[core.nodes[v].parents[0]]&&done[core.nodes[v].parents[1]]))){done[v]=1;out.order.push_back(v);advanced=true;}if(!advanced){out.mask|=invalid_dag;return out;}}
 std::vector<unsigned char> ancestor(n);for(int v:out.axes)ancestor[v]=1;
 for(auto at=out.order.rbegin();at!=out.order.rend();++at)if(ancestor[*at]&&core.nodes[*at].has_birth)for(int p:core.nodes[*at].parents)ancestor[p]=1;
 for(int v=0;v<n;++v)if(!typed[v]&&!ancestor[v])out.mask|=nonancestral;
 return out;
}
using Partition=std::array<int,4>;
struct Identity {std::map<Partition,std::uint64_t> counts;std::uint64_t paths=0;};
inline Identity enumerate(const ::relationship_topology::TopologyModel& core,const Certificate& cert){
 if(cert.mask||cert.axes.size()!=2)throw std::runtime_error("uncertified four-gene topology");
 std::vector<std::array<int,2>> genes(core.nodes.size());int founder=0,births=0;
 for(int v:cert.order)if(!core.nodes[v].has_birth)genes[v]={{founder++,founder++}};else ++births;
 if(births>10)throw std::runtime_error("structural segregation budget exceeded");
 Identity out;std::function<void(size_t)> visit=[&](size_t at){if(at==cert.order.size()){
  Partition origins{{genes[cert.axes[0]][0],genes[cert.axes[0]][1],genes[cert.axes[1]][0],genes[cert.axes[1]][1]}},state{};
  std::map<int,int> labels;for(int i=0;i<4;++i){auto [a,fresh]=labels.emplace(origins[i],labels.size());state[i]=a->second;}++out.counts[state];++out.paths;return;
 };int v=cert.order[at];const auto& node=core.nodes[v];if(!node.has_birth){visit(at+1);return;}
 for(int a=0;a<2;++a)for(int b=0;b<2;++b){genes[v]={{genes[node.parents[0]][a],genes[node.parents[1]][b]}};visit(at+1);}};visit(0);return out;
}
inline std::vector<double> table(const Identity& identity,const std::vector<double>& p){
 int k=p.size(),G=k*(k+1)/2;std::vector<int> index(k*k);int g=0;for(int a=0;a<k;++a)for(int b=a;b<k;++b)index[a*k+b]=index[b*k+a]=g++;
 std::vector<double> out(G*G);for(const auto& [state,count]:identity.counts){int groups=*std::max_element(state.begin(),state.end())+1;std::array<int,4> values{};
 std::function<void(int,double)> fill=[&](int at,double w){if(at==groups){int a=index[values[state[0]]*k+values[state[1]]],b=index[values[state[2]]*k+values[state[3]]];out[a+G*b]+=w;return;}for(int x=0;x<k;++x){values[at]=x;fill(at+1,w*p[x]);}};fill(0,double(count)/identity.paths);}
 return out;
}
}
#endif
