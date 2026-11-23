#include "pap_removed_path.h"
#ifndef PAP_MESSAGE_GRAPH_H
#define PAP_MESSAGE_GRAPH_H
#include <chrono>
// Numerical provenance: stage2_p7_collateral_messages.h at 80ae1813,
// dependencies/compute_availability/get/compute (to/from sum-product), and
// scale retention. No Stage-2 orchestration, publication or cache is imported.
#include "relationship_topology.h"
#include "log_products.h"
#include "evidence_ledger.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace pap_internal::pedigree_message {
constexpr double neginf=-std::numeric_limits<double>::infinity();
inline double add(double a,double b){if(a==neginf)return b;if(b==neginf)return a;double hi=std::max(a,b);return hi+std::log1p(std::exp(std::min(a,b)-hi));}
struct Refusal:std::runtime_error {using std::runtime_error::runtime_error;};
struct Locus {
 int k=0,G=0;
 std::vector<double>prior;
 std::vector<std::vector<double>>emission;
 std::vector<std::array<int,2>>alleles;
 std::vector<int>index;
 int genotype(int a,int b)const{return index[a*k+b];}
 template<class F>void support(int a,int b,F&& fn)const {
  for(int x:alleles[a])for(int y:alleles[b])fn(genotype(x,y),-std::log(4.0));
 }
};
struct Factor {
 std::vector<int>scope;
 std::vector<double>value;
 double scale=0;
 void normalize(){double maximum=neginf;for(double v:value){if(std::isnan(v))throw Refusal("unavailable factor value");maximum=std::max(maximum,v);}if(maximum==neginf)return;for(double& v:value)v-=maximum;scale+=maximum;}
};
enum class Kind {phi,prior,birth};
struct Atom {Kind kind;int child=-1,dam=-1,sire=-1;std::string key;};
inline std::string identity(const relationship_topology::TopologyNode& n){return n.typed_id?"T"+std::to_string(n.typed_id):"L"+std::to_string(n.label.size())+":"+n.label;}
inline std::string encoded(const std::string& x){return std::to_string(x.size())+":"+x;}
inline std::string birth_key(const relationship_topology::TopologyModel& m,int c,int d,int s){return "Birth:"+encoded(identity(m.nodes[c]))+encoded(identity(m.nodes[d]))+encoded(identity(m.nodes[s]));}
inline std::vector<Atom> structural_atoms(const relationship_topology::TopologyModel& m){
 std::vector<Atom> out;
 for(size_t i=0;i<m.nodes.size();++i){const auto& n=m.nodes[i];
  if(n.observed&&!n.neutral_phi)out.push_back({Kind::phi,int(i),-1,-1,{}});
  if(!n.has_birth)out.push_back({Kind::prior,int(i),-1,-1,{}});
  else out.push_back({Kind::birth,int(i),n.parents[0],n.parents[1],{}});
 }return out;
}
inline std::vector<Atom> atoms(const relationship_topology::TopologyModel& m){
 auto out=structural_atoms(m);for(auto& a:out)
  a.key=a.kind==Kind::birth?birth_key(m,a.child,a.dam,a.sire):
    (a.kind==Kind::phi?"Phi:":"Prior:")+identity(m.nodes[a.child]);
 return out;
}
inline std::vector<int> scope(const Atom& a){std::vector<int>v{a.child};if(a.kind==Kind::birth){v.push_back(a.dam);v.push_back(a.sire);}std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());return v;}
inline size_t cells(int G,size_t width,size_t limit=2000000U){size_t n=1;for(size_t i=0;i<width;++i){if(n>limit/size_t(G))throw Refusal("largest intermediate exceeds "+std::to_string(limit)+" cells");n*=G;}return n;}
inline size_t offset(const std::vector<int>& vars,const std::vector<int>& assignment,int G){size_t index=0,stride=1;for(int v:vars){index+=assignment[v]*stride;stride*=G;}return index;}
inline Factor materialize(const Atom& a,const relationship_topology::TopologyModel& model,const Locus& l){
 Factor f;f.scope=scope(a);f.value.assign(cells(l.G,f.scope.size()),neginf);
 if(a.kind==Kind::phi){int id=model.nodes[a.child].typed_id;if(id<=0||id>=int(l.emission.size()))throw Refusal("unknown emission identity");f.value=l.emission[id];}
 else if(a.kind==Kind::prior)f.value=l.prior;
 else {std::vector<int>g(model.nodes.size());for(int d=0;d<l.G;++d)for(int s=0;s<l.G;++s){if(a.dam==a.sire&&d!=s)continue;g[a.dam]=d;g[a.sire]=s;l.support(d,s,[&](int c,double weight){g[a.child]=c;auto i=offset(f.scope,g,l.G);f.value[i]=add(f.value[i],weight);});}}
 f.normalize();return f;
}
inline Factor eliminate(std::vector<Factor>f,std::set<int>remaining,int G,size_t work_limit=50000000U,size_t cell_limit=2000000U){
 double scalar=0;size_t work=0;std::vector<int>assignment;
 int largest=-1;for(const auto& x:f)for(int v:x.scope)largest=std::max(largest,v);assignment.resize(largest+1);
 while(!remaining.empty()){
  int best=-1;std::vector<int>joined;size_t width=SIZE_MAX;
  for(int v:remaining){std::set<int>u;for(const auto& x:f)if(std::find(x.scope.begin(),x.scope.end(),v)!=x.scope.end())u.insert(x.scope.begin(),x.scope.end());if(u.size()<width){width=u.size();best=v;joined.assign(u.begin(),u.end());}}
  std::vector<Factor>others,touch;for(auto& x:f)(std::find(x.scope.begin(),x.scope.end(),best)==x.scope.end()?others:touch).push_back(std::move(x));
  if(touch.empty()){scalar+=std::log(double(G));remaining.erase(best);f=std::move(others);continue;}
  size_t count=cells(G,joined.size(),cell_limit);size_t charge=count*(touch.size()+1);if(charge>work_limit-work)throw Refusal("elimination work exceeds "+std::to_string(work_limit)+" units");work+=charge;
  Factor result;result.scope=joined;result.scope.erase(std::find(result.scope.begin(),result.scope.end(),best));result.value.assign(cells(G,result.scope.size(),cell_limit),neginf);for(const auto& x:touch)result.scale+=x.scale;
  // Precompute positional strides; no graph/scope search in the cell loop.
  std::vector<std::vector<size_t>>strides(touch.size(),std::vector<size_t>(joined.size()));
  for(size_t t=0;t<touch.size();++t){size_t stride=1;for(int v:touch[t].scope){strides[t][std::find(joined.begin(),joined.end(),v)-joined.begin()]=stride;stride*=G;}}
  size_t removed=std::find(joined.begin(),joined.end(),best)-joined.begin();std::vector<int>digits(joined.size());
  for(size_t i=0;i<count;++i){size_t q=i;for(size_t j=0;j<digits.size();++j){digits[j]=q%G;q/=G;}double value=0;
   for(size_t t=0;t<touch.size();++t){size_t index=0;for(size_t j=0;j<digits.size();++j)index+=digits[j]*strides[t][j];value+=touch[t].value[index];}
   size_t index=0,stride=1;for(size_t j=0;j<digits.size();++j)if(j!=removed){index+=digits[j]*stride;stride*=G;}result.value[index]=add(result.value[index],value);
  }
  result.normalize();others.push_back(std::move(result));f=std::move(others);remaining.erase(best);
 }
 std::set<int>u;for(const auto& x:f)u.insert(x.scope.begin(),x.scope.end());Factor out;out.scope.assign(u.begin(),u.end());out.value.assign(cells(G,out.scope.size(),cell_limit),scalar);for(const auto& x:f)out.scale+=x.scale;
 for(size_t i=0;i<out.value.size();++i){size_t q=i;for(int v:out.scope){assignment[v]=q%G;q/=G;}for(const auto& x:f)out.value[i]+=x.value[offset(x.scope,assignment,G)];}
 out.normalize();return out;
}
struct WorkBudget {size_t used=0,limit=50000000;void charge(size_t n){if(n>limit-used)throw Refusal("bounded family exterior work cap");used+=n;}};
inline thread_local WorkBudget* active_work_budget=nullptr;
struct WorkScope {WorkBudget* previous;explicit WorkScope(WorkBudget* p):previous(active_work_budget){active_work_budget=p;}~WorkScope(){active_work_budget=previous;}};
inline void charge_work(size_t n){if(active_work_budget)active_work_budget->charge(n);}
struct SecondsSpan {double* total;std::chrono::steady_clock::time_point start;explicit SecondsSpan(double* t):total(t){if(total)start=std::chrono::steady_clock::now();}~SecondsSpan(){if(total)*total+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}};
struct Cached {std::vector<double>value;double scale=0;std::vector<std::string>axes;};
class Store {
 std::map<std::string,Cached>values_;size_t bytes_=0;
 std::map<std::string,uint32_t> symbols_;std::vector<std::string> names_;size_t symbol_bytes_=0;
 public:double key_seconds=0,binding_seconds=0,arithmetic_seconds=0;int key_depth=0;bool profile_arithmetic=false;bool compact_keys=false,gamete_factorized=false,directed_cache=false,family_products=false;size_t family_product_builds=0,family_product_queries=0;size_t computed=0,reused=0;size_t byte_limit=64U*1024U*1024U;
 uint32_t symbol(const std::string& x){auto at=symbols_.find(x);if(at!=symbols_.end())return at->second;if(compact_keys&&(names_.size()>=8192||x.size()+96>1048576-symbol_bytes_))throw Refusal("bounded message symbol capacity");uint32_t id=names_.size();names_.push_back(x);try{symbols_.emplace(x,id);}catch(...){names_.pop_back();throw;}symbol_bytes_+=x.size()+96;return id;}
 const std::string& name(uint32_t id)const{return names_.at(id);}
 const Cached*find(const std::string& key){auto it=values_.find(key);if(it==values_.end())return nullptr;++reused;return &it->second;}
 void save(const std::string& key,Cached value){++computed;size_t bytes=key.size()+sizeof(Cached)+value.value.size()*sizeof(double);for(const auto& x:value.axes)bytes+=x.size();if(bytes>byte_limit-bytes_)return;if(values_.emplace(key,std::move(value)).second)bytes_+=bytes;}
 // The same dependency identity used by directed forest messages. This is
 // a lookup into the existing Store, not a second numerical cache. Scales and
 // structural zeros are returned unchanged; absent factors stay excluded.
 std::string unary_key(std::vector<uint32_t> deps,uint32_t axis,const std::string& binding,const char* kind){ pap_internal::removed_path("unary_key"); }
 const Cached* unary(const std::vector<uint32_t>& deps,uint32_t axis,const std::string& binding){for(auto kind:{"from","to","cavity"})if(auto hit=find(unary_key(deps,axis,binding,kind)))return hit;return nullptr;}
 size_t bytes()const{return bytes_+symbol_bytes_;}
};
struct ArithmeticSpan {Store& owner;bool running;std::chrono::steady_clock::time_point start;explicit ArithmeticSpan(Store& x):owner(x),running(x.profile_arithmetic){if(running)start=std::chrono::steady_clock::now();}void stop(){if(running){owner.arithmetic_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();running=false;}}~ArithmeticSpan(){stop();}};
struct Family {int dam=-1,sire=-1;std::vector<int>children,members;std::vector<uint32_t> birth_ids;};
// Immutable family incidence plan. It carries no locus values and can be
// shared by all loci of one exact factor/cavity context.
struct TreeStructure {
 struct Group {int dam,sire;std::vector<int> children,members;std::vector<size_t> births;};
 std::vector<Atom> atoms;std::vector<Group> families;std::vector<std::vector<int>> incident;
 std::set<int> variables;bool tree=false;size_t nodes=0;
 TreeStructure(const relationship_topology::TopologyModel& m,const std::vector<Atom>& source):atoms(source),incident(m.nodes.size()),nodes(m.nodes.size()){
  std::map<std::pair<int,int>,int> index;
  for(size_t i=0;i<atoms.size();++i){const auto& a=atoms[i];for(int v:scope(a))variables.insert(v);if(a.kind!=Kind::birth)continue;if(a.dam==a.sire)return;
   auto [at,fresh]=index.emplace(std::make_pair(a.dam,a.sire),families.size());if(fresh)families.push_back({a.dam,a.sire,{},{a.dam,a.sire},{}});auto& f=families[at->second];f.children.push_back(a.child);f.members.push_back(a.child);f.births.push_back(i);
  }
  size_t edges=0;for(size_t f=0;f<families.size();++f){auto& members=families[f].members;std::sort(members.begin(),members.end());members.erase(std::unique(members.begin(),members.end()),members.end());edges+=members.size();for(int v:members)incident[v].push_back(f);}
  tree=!variables.empty()&&edges+1==variables.size()+families.size();
  if(tree){std::set<int> seen{*variables.begin()};std::deque<int> queue{*variables.begin()};while(!queue.empty()){int v=queue.front();queue.pop_front();if(v<int(nodes)){for(int f:incident[v])if(seen.insert(nodes+f).second)queue.push_back(nodes+f);}else for(int node:families[v-nodes].members)if(seen.insert(node).second)queue.push_back(node);}tree=seen.size()==variables.size()+families.size();}
 }
 size_t bytes()const{size_t n=sizeof(*this)+incident.size()*sizeof(std::vector<int>)+atoms.size()*sizeof(Atom);for(const auto& a:atoms)n+=a.key.size();for(const auto& v:incident)n+=v.size()*sizeof(int);for(const auto& f:families)n+=sizeof(f)+(f.children.size()+f.members.size())*sizeof(int)+f.births.size()*sizeof(size_t);return n+variables.size()*48;}
};
#include "indexed_messages.h"
// A directed sum-product view of one EXACT factor set. Core-owned Phi/prior
// and births are absent, rather than masked after a message was constructed.
class Tree {
 IndexedForestMessages* indexed_=nullptr;
 const relationship_topology::TopologyModel& model_;const Locus& l_;Store& cache_;std::string binding_;
 std::vector<std::vector<double>>local_;
 std::vector<std::vector<uint32_t>>local_keys_;
 std::vector<Family>families_;
 std::vector<std::vector<int>>incident_;
 std::set<int>variables_;
 std::map<std::pair<int,int>,Factor>directed_;
 std::vector<std::unique_ptr<pap_internal::log_products::Index>> child_products_;std::vector<std::vector<size_t>> pending_children_;std::map<int,size_t> child_slots_;
 std::map<std::pair<int,int>,std::vector<uint32_t>>dependencies_;
 bool tree_=false;
 std::vector<uint32_t>keys(int node,int excluded){SecondsSpan timing(cache_.profile_arithmetic&&cache_.key_depth++==0?&cache_.key_seconds:nullptr);struct Depth{int& d;bool on;~Depth(){if(on)--d;}} depth{cache_.key_depth,cache_.profile_arithmetic};auto at=dependencies_.find({node,excluded});if(at!=dependencies_.end()){charge_work(at->second.size());return at->second;}
  auto out=local_keys_[node];for(int f:incident_[node])if(f!=excluded){out.insert(out.end(),families_[f].birth_ids.begin(),families_[f].birth_ids.end());for(int v:families_[f].members)if(v!=node){auto more=keys(v,f);out.insert(out.end(),more.begin(),more.end());}}
  charge_work(out.size()*(1+size_t(std::ceil(std::log2(double(out.size()+1))))));std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());dependencies_[{node,excluded}]=out;return out;
 }
 std::string key(const std::vector<uint32_t>& factors,const std::vector<int>& axes,const char* kind){SecondsSpan timing(cache_.profile_arithmetic?&cache_.key_seconds:nullptr);charge_work(factors.size()+axes.size());std::string s=binding_+kind;
  if(cache_.compact_keys){auto append=[&](uint32_t value){s.append(reinterpret_cast<const char*>(&value),sizeof(value));};append(factors.size());for(auto id:factors)append(id);append(axes.size());for(int v:axes)append(cache_.symbol("Axis:"+identity(model_.nodes[v])));}
  else{std::vector<std::string> names;for(auto id:factors)names.push_back(cache_.name(id));std::sort(names.begin(),names.end());for(const auto& x:names)s+=std::to_string(x.size())+":"+x;for(int v:axes){auto x=identity(model_.nodes[v]);s+="|"+std::to_string(x.size())+":"+x;}}
  return s;}

 Factor local_cavity(int node,const std::set<int>& excluded){charge_work(size_t(l_.G)*(incident_[node].size()+1));Factor out;out.scope={node};out.value=local_[node];for(int f:incident_[node])if(!excluded.count(f)){auto m=from(f,node);ArithmeticSpan math(cache_);for(int g=0;g<l_.G;++g)out.value[g]+=m.value[g];out.scale+=m.scale;}out.normalize();return out;}
 Factor to(int node,int excluded){
  if(indexed_){size_t slot=indexed_->slot(excluded,node);if(indexed_->to_ready[slot]){++indexed_->reused;const auto& hit=indexed_->to[slot];return {{node},hit.value,hit.scale};}
   auto out=local_cavity(node,{excluded});indexed_->to[slot]={out.value,out.scale,{}};indexed_->to_ready[slot]=1;++indexed_->computed;return out;}
  auto at=directed_.find({node,excluded});if(at!=directed_.end())return at->second;auto k=key(keys(node,excluded),{node},"to");
  if(auto hit=cache_.find(k)){Factor out;out.scope={node};out.value=hit->value;out.scale=hit->scale;directed_[{node,excluded}]=out;return out;}
  auto out=local_cavity(node,{excluded});cache_.save(k,{out.value,out.scale,{identity(model_.nodes[node])}});directed_[{node,excluded}]=out;return out;
 }
 // Donor compute(): sibling likelihoods condition on one shared parent pair;
 // the other parent's cavity weights are integrated exactly, with scale kept.
 Factor pair_kernel(int f,int first,int second){ pap_internal::removed_path("pair_kernel"); }
 Factor from(int f,int node){const auto& F=families_[f];std::string saved_key;
  if(indexed_){size_t slot=indexed_->slot(f,node);if(indexed_->from_ready[slot]){++indexed_->reused;const auto& hit=indexed_->from[slot];return {{node},hit.value,hit.scale};}}
  if(!indexed_&&cache_.directed_cache){auto deps=F.birth_ids;for(int v:F.members)if(v!=node){auto more=keys(v,f);deps.insert(deps.end(),more.begin(),more.end());}std::sort(deps.begin(),deps.end());deps.erase(std::unique(deps.begin(),deps.end()),deps.end());saved_key=key(deps,{node},"from");if(auto hit=cache_.find(saved_key)){Factor out;out.scope={node};out.value=hit->value;out.scale=hit->scale;return out;}}
int partner=-1;for(int v:F.members)if(v!=node){partner=v;break;}if(partner<0)throw Refusal("degenerate family incidence");
  charge_work(size_t(l_.G)*l_.G);auto pair=pair_kernel(f,node,partner);auto m=to(partner,f);ArithmeticSpan sum_math(cache_);Factor out;out.scope={node};out.value.assign(l_.G,neginf);out.scale=pair.scale+m.scale;
  for(int a=0;a<l_.G;++a)for(int b=0;b<l_.G;++b)out.value[a]=add(out.value[a],pair.value[a+size_t(l_.G)*b]+m.value[b]);out.normalize();sum_math.stop();
  if(indexed_){size_t slot=indexed_->slot(f,node);indexed_->from[slot]={out.value,out.scale,{}};indexed_->from_ready[slot]=1;++indexed_->computed;}
  else if(cache_.directed_cache)cache_.save(saved_key,{out.value,out.scale,{identity(model_.nodes[node])}});return out;
 }
 public:
 Tree(const relationship_topology::TopologyModel& m,const std::vector<Atom>& atoms,const Locus& l,Store& cache,std::string binding,
      const std::vector<Factor>& unary={}):model_(m),l_(l),cache_(cache),binding_(std::move(binding)),local_(m.nodes.size(),std::vector<double>(l.G,0)),local_keys_(m.nodes.size()),incident_(m.nodes.size()){
  SecondsSpan binding_time(cache_.profile_arithmetic?&cache_.binding_seconds:nullptr);charge_work(size_t(l.G)*m.nodes.size()+8*atoms.size());
  std::map<std::pair<int,int>,int>family_index;
  for(const auto& a:atoms){for(int v:scope(a))variables_.insert(v);
   if(a.kind!=Kind::birth){auto factor=materialize(a,m,l);for(int g=0;g<l.G;++g)local_[a.child][g]+=factor.value[g]+factor.scale;local_keys_[a.child].push_back(cache_.symbol(a.key));continue;}
   if(a.dam==a.sire)return;
   auto key=std::make_pair(a.dam,a.sire);auto at=family_index.find(key);int f;
   if(at==family_index.end()){f=families_.size();family_index[key]=f;families_.push_back({a.dam,a.sire,{},{a.dam,a.sire}});}else f=at->second;
   families_[f].children.push_back(a.child);families_[f].members.push_back(a.child);families_[f].birth_ids.push_back(cache_.symbol(birth_key(m,a.child,a.dam,a.sire)));
  }
  for(const auto& factor:unary){if(factor.scope.size()!=1)return;int v=factor.scope[0];variables_.insert(v);for(int g=0;g<l.G;++g)local_[v][g]+=factor.value[g]+factor.scale;}
  size_t edges=0;for(size_t f=0;f<families_.size();++f){auto& members=families_[f].members;std::sort(members.begin(),members.end());members.erase(std::unique(members.begin(),members.end()),members.end());edges+=members.size();for(int v:members)incident_[v].push_back(f);}
  // A connected exterior piece is a tree iff its incidence edges = vertices-1.
  tree_=edges+1==variables_.size()+families_.size();
  if(tree_&&!variables_.empty()){int N=m.nodes.size();std::set<int>seen{*variables_.begin()};std::deque<int>todo{*variables_.begin()};while(!todo.empty()){int v=todo.front();todo.pop_front();if(v<N){for(int f:incident_[v])if(seen.insert(N+f).second)todo.push_back(N+f);}else for(int node:families_[v-N].members)if(seen.insert(node).second)todo.push_back(node);}tree_=seen.size()==variables_.size()+families_.size();}
 }
 Tree(const relationship_topology::TopologyModel& m,const TreeStructure& plan,const Locus& l,Store& cache,std::string binding,const std::vector<Factor>& unary={},IndexedForestMessages* indexed=nullptr)
  :indexed_(indexed),model_(m),l_(l),cache_(cache),binding_(std::move(binding)),local_(m.nodes.size(),std::vector<double>(l.G,0)),local_keys_(m.nodes.size()),incident_(plan.incident),variables_(plan.variables),tree_(plan.tree){
  SecondsSpan binding_time(cache_.profile_arithmetic?&cache_.binding_seconds:nullptr);charge_work(size_t(l.G)*m.nodes.size()+8*plan.atoms.size());
  for(const auto& a:plan.atoms)if(a.kind!=Kind::birth){auto factor=materialize(a,m,l);for(int g=0;g<l.G;++g)local_[a.child][g]+=factor.value[g]+factor.scale;if(!indexed_)local_keys_[a.child].push_back(cache_.symbol(a.key));}
  for(const auto& f:plan.families){families_.push_back({f.dam,f.sire,f.children,f.members,{}});if(!indexed_)for(size_t index:f.births){const auto& a=plan.atoms[index];families_.back().birth_ids.push_back(cache_.symbol(birth_key(m,a.child,a.dam,a.sire)));}}
  for(const auto& factor:unary){if(factor.scope.size()!=1){tree_=false;return;}int v=factor.scope[0];if(!variables_.count(v)){tree_=false;return;}for(int g=0;g<l.G;++g)local_[v][g]+=factor.value[g]+factor.scale;}
 }
 // Publish both directions of every family incidence using products of the
 // remaining incoming messages. This never divides by a posterior.
 void sweep(){if(!tree_)throw Refusal("directed forest component is cyclic");for(int f=0;f<int(families_.size());++f)for(int v:families_[f].members){to(v,f);from(f,v);}}
 bool available()const{return tree_;}
 Factor marginal(std::vector<int>axes){if(!tree_)throw Refusal("directed dependency closure is cyclic");
  if(axes.size()<=1){int root=axes.empty()?*variables_.begin():axes[0];auto out=local_cavity(root,{});
   if(axes.empty()){double value=neginf;for(double v:out.value)value=add(value,v);out.scope.clear();out.value={value};out.normalize();}return out;}
  if(axes.size()!=2)throw Refusal("higher attachment uses exact elimination");
  int a=axes[0],b=axes[1],N=model_.nodes.size();std::vector<int>parent(N+families_.size(),-1);std::deque<int>todo{a};parent[a]=a;
  while(!todo.empty()&&parent[b]<0){int v=todo.front();todo.pop_front();if(v<N){for(int f:incident_[v])if(parent[N+f]<0){parent[N+f]=v;todo.push_back(N+f);}}else for(int node:families_[v-N].members)if(parent[node]<0){parent[node]=v;todo.push_back(node);}}
  if(parent[b]<0)throw Refusal("joint axes are disconnected");std::vector<int>path;for(int v=b;;v=parent[v]){path.push_back(v);if(v==a)break;}std::reverse(path.begin(),path.end());
  auto start=local_cavity(a,{path[1]-N});auto matrix=pair_kernel(path[1]-N,a,path[2]);matrix.scale+=start.scale;
  for(int x=0;x<l_.G;++x)for(int y=0;y<l_.G;++y)matrix.value[x+size_t(l_.G)*y]+=start.value[x];matrix.normalize();
  for(size_t i=2;i+2<path.size();i+=2){charge_work(size_t(l_.G)*l_.G*l_.G);int node=path[i];auto local=local_cavity(node,{path[i-1]-N,path[i+1]-N});auto next=pair_kernel(path[i+1]-N,node,path[i+2]);ArithmeticSpan path_math(cache_);Factor result;result.scope={a,path[i+2]};result.scale=matrix.scale+local.scale+next.scale;result.value.assign(size_t(l_.G)*l_.G,neginf);
   for(int x=0;x<l_.G;++x)for(int y=0;y<l_.G;++y)for(int g=0;g<l_.G;++g)result.value[x+size_t(l_.G)*y]=add(result.value[x+size_t(l_.G)*y],matrix.value[x+size_t(l_.G)*g]+local.value[g]+next.value[g+size_t(l_.G)*y]);result.normalize();matrix=std::move(result);
  }
  auto end=local_cavity(b,{path[path.size()-2]-N});matrix.scale+=end.scale;for(int x=0;x<l_.G;++x)for(int y=0;y<l_.G;++y)matrix.value[x+size_t(l_.G)*y]+=end.value[y];matrix.scope={a,b};matrix.normalize();return matrix;
 }
};
struct Result {int state=3;double log_partition=std::numeric_limits<double>::quiet_NaN();std::string reason;};
struct PiecePlan {std::vector<Atom>atoms;std::vector<int>inside,axes;std::string key;mutable IndexedForestProjection indexed;};
struct Plan {std::vector<Atom>core;std::vector<PiecePlan>pieces;};
// One validated hypothesis, one ordered factor inventory. Fixed scopes retain
// the reference sorted/unique variable order, including repeated-parent births.
struct AtomScope {
 std::array<int,3> nodes{};unsigned char width=0;
 explicit AtomScope(const Atom& a){nodes={a.child,a.dam,a.sire};width=a.kind==Kind::birth?3:1;std::sort(nodes.begin(),nodes.begin()+width);width=std::unique(nodes.begin(),nodes.begin()+width)-nodes.begin();}
 explicit AtomScope(const std::vector<int>& v){width=v.size();std::copy(v.begin(),v.end(),nodes.begin());}
 auto begin()const{return nodes.begin();}auto end()const{return nodes.begin()+width;}size_t size()const{return width;}
 int operator[](size_t i)const{return nodes[i];}
};
struct FactorInventory {
 struct Record {AtomScope scope;int core_mask;};
 std::vector<Atom> atoms;std::vector<Record> records;std::vector<unsigned char> inside;std::vector<int> lookup;
 FactorInventory(const relationship_topology::TopologyModel& m,const std::set<int>& core):atoms(structural_atoms(m)),inside(m.nodes.size()),lookup(3*m.nodes.size(),-1){
  for(int v:core)inside.at(v)=1;records.reserve(atoms.size());
  for(const auto& a:atoms){AtomScope vars(a);int mask=0;for(size_t j=0;j<vars.size();++j)if(inside.at(vars[j]))mask|=1<<j;
   lookup.at(size_t(a.kind)*inside.size()+a.child)=records.size();records.push_back({vars,mask});}
 }
  bool owns(const Plan& plan)const {
  std::vector<int> counts(atoms.size());bool good=true;
  auto consume=[&](const Atom& a,bool core){if(a.child<0||size_t(a.child)>=inside.size()||int(a.kind)<0||int(a.kind)>2){good=false;return;}int i=lookup[size_t(a.kind)*inside.size()+a.child];if(i<0){good=false;return;}const auto& original=atoms[i];good&=a.dam==original.dam&&a.sire==original.sire;if(core)good&=records[i].core_mask==(1<<records[i].scope.size())-1;++counts[i];};
  for(const auto& a:plan.core)consume(a,true);for(const auto& p:plan.pieces)for(const auto& a:p.atoms)consume(a,false);
  for(int count:counts)good&=count==1;return good;
 }
};
inline Plan compile(const relationship_topology::TopologyModel& m,const std::set<int>& explicit_core,bool named=true,const FactorInventory* prepared=nullptr){
 if(named)prepared=nullptr;
 Plan out;std::vector<Atom> owned;if(!prepared)owned=named?atoms(m):structural_atoms(m);const auto& all=prepared?prepared->atoms:owned;auto member=[&](int v){return prepared?bool(prepared->inside.at(v)):explicit_core.count(v)>0;};std::vector<Atom>outside;std::map<int,std::set<int>>adj;
 for(size_t ai=0;ai<all.size();++ai){const auto& a=all[ai];auto vars=prepared?prepared->records[ai].scope:AtomScope(scope(a));bool external=false;for(int v:vars)external|=!member(v);if(!external){out.core.push_back(a);continue;}outside.push_back(a);for(int v:vars)if(!member(v))for(int w:vars)if(!member(w))adj[v].insert(w);}
 std::set<int>seen;
 for(const auto& item:adj)if(!seen.count(item.first)){std::set<int>inside{item.first};std::deque<int>todo{item.first};seen.insert(item.first);while(!todo.empty()){int v=todo.front();todo.pop_front();for(int w:adj[v])if(seen.insert(w).second){inside.insert(w);todo.push_back(w);}}
  PiecePlan piece;piece.inside.assign(inside.begin(),inside.end());std::set<int>axes;std::vector<std::string>keys;
  for(const auto& a:outside){auto vars=prepared?prepared->records.at(prepared->lookup.at(size_t(a.kind)*m.nodes.size()+a.child)).scope:AtomScope(scope(a));bool owns=false;for(int v:vars)owns|=inside.count(v);if(!owns)continue;piece.atoms.push_back(a);if(named)keys.push_back(a.key);for(int v:vars)if(member(v))axes.insert(v);}
  piece.axes.assign(axes.begin(),axes.end());if(named)std::sort(piece.axes.begin(),piece.axes.end(),[&](int a,int b){return identity(m.nodes[a])<identity(m.nodes[b]);});std::sort(keys.begin(),keys.end());for(const auto& key:keys)piece.key+=encoded(key);if(named){piece.key+="|axes:";for(int v:piece.axes)piece.key+=encoded(identity(m.nodes[v]));}out.pieces.push_back(std::move(piece));
 }
 // A literal factor cannot occur in both a potential and the explicit core.
 if(named){std::map<std::string,int>counts;for(const auto& a:out.core)++counts[a.key];for(const auto& p:out.pieces)for(const auto& a:p.atoms)++counts[a.key];
 if(counts.size()!=all.size())throw Refusal("factor identity alias or duplicate ownership");for(const auto& row:counts)if(row.second!=1)throw Refusal("factor counted more than once");}
 return out;
}
class Engine {
 public:explicit Engine(std::vector<Locus> data):loci(std::move(data)){}
 const std::vector<Locus>loci;Store store;size_t tree_pieces=0,joint_pieces=0,absorbed_pieces=0;
 private:
 Factor piece(const relationship_topology::TopologyModel& m,const PiecePlan& p,size_t locus){const auto& l=loci.at(locus);auto key=std::to_string(locus)+"|piece:"+p.key;
  if(auto cached=store.find(key)){Factor out;out.scope=p.axes;out.value=cached->value;out.scale=cached->scale;return out;}
  Factor value;Tree messages(m,p.atoms,l,store,std::to_string(locus)+"|");
  if(messages.available()&&p.axes.size()<=2){value=messages.marginal(p.axes);++tree_pieces;if(p.axes.size()==2)++joint_pieces;}
  else {std::vector<Factor>f;for(const auto& a:p.atoms)f.push_back(materialize(a,m,l));value=eliminate(std::move(f),std::set<int>(p.inside.begin(),p.inside.end()),l.G);++absorbed_pieces;}
  // Generic elimination orders its axes numerically. Rebind to the signature's
  // biological order before caching; model-local node numbers never key values.
  if(value.scope!=p.axes){Factor ordered;ordered.scope=p.axes;ordered.value.resize(value.value.size());ordered.scale=value.scale;std::vector<int>assignment(m.nodes.size());for(size_t i=0;i<ordered.value.size();++i){size_t q=i;for(int v:p.axes){assignment[v]=q%l.G;q/=l.G;}ordered.value[i]=value.value[offset(value.scope,assignment,l.G)];}value=std::move(ordered);}
  std::vector<std::string>names;for(int v:p.axes)names.push_back(identity(m.nodes[v]));store.save(key,{value.value,value.scale,names});return value;
 }
 // Exact leaf-child contraction, from the donor's child_factor. This avoids
 // materializing G^3 tables when a child has only local genotype evidence.
 Factor core_value(const relationship_topology::TopologyModel& m,std::vector<Atom>specs,
                   std::vector<Factor>extras,size_t locus,const std::string& token){ pap_internal::removed_path("core_value"); }
 public:
 Result score(const relationship_topology::TopologyModel& m,const std::set<int>& core,size_t locus,bool reference=false,size_t reference_work_limit=50000000U){
  if(!m.ready())return {m.status==relationship_topology::TopologyStatus::Inapplicable?0:3,m.status==relationship_topology::TopologyStatus::Inapplicable?neginf:std::numeric_limits<double>::quiet_NaN(),m.reason};
  try {const auto& l=loci.at(locus);if(l.k<2||l.k>9)throw Refusal("explicit k outside 2..9 refusal");Factor value;
   if(reference){std::vector<Factor>f;std::set<int>vars;for(const auto& a:atoms(m)){f.push_back(materialize(a,m,l));for(int v:scope(a))vars.insert(v);}value=eliminate(std::move(f),vars,l.G,reference_work_limit);}
   else {if(!core.count(m.focal))throw Refusal("focal genotype must remain explicit");auto plan=compile(m,core);std::vector<Factor>potentials;std::string token;for(const auto& p:plan.pieces){potentials.push_back(piece(m,p,locus));token+=encoded(p.key);}for(const auto& a:plan.core)token+=encoded(a.key);value=core_value(m,std::move(plan.core),std::move(potentials),locus,token);}
   if(!value.scope.empty()||value.value.size()!=1)throw Refusal("partition retains an unexpected axis");double z=value.value[0]+value.scale;
   if(z==neginf)return {2,neginf,"structural zero"};if(!std::isfinite(z))throw Refusal("unavailable partition");return {1,z,"exact declared factor model"};
  }catch(const Refusal& e){return {3,std::numeric_limits<double>::quiet_NaN(),e.what()};}
 }
};
}
#endif
