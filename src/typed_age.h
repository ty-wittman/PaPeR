#ifndef PAP_TYPED_AGE_H
#define PAP_TYPED_AGE_H

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <deque>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace pap_internal { namespace typed_age {
using Year = std::int64_t;
constexpr Year infinity = Year(1) << 55;
struct Bounds { std::vector<int> lo, hi; };
// Main-thread visit lifetime, also used by the pedigree likelihood catalog. Numerical workers
// consume already-built immutable models, never this mutable typed state.
inline thread_local const Bounds* active = nullptr;
struct Scope {
  const Bounds* previous;
  explicit Scope(const Bounds* bounds):previous(active) { active=bounds; }
  ~Scope() { active=previous; }
};
struct Patch {
  bool feasible=true;
  std::vector<int> nodes;
  std::vector<Year> lo, hi;
  std::size_t relaxations=0;
};

// Difference constraints on typed individuals only. At the fixed point the
// lower (upper) vector itself is feasible and attains every marginal minimum
// (maximum). Thus these are exact interval projections, not heuristic ages.
class State {
 public:
  std::vector<std::array<int,2>> parents; // zero-based; -1 means absent
  std::vector<std::vector<int>> children;
  std::vector<Year> original_lo, original_hi, lo, hi, last;
  std::array<Year,2> minimum, maximum;
  bool feasible=true;

  State(const std::vector<int>& lower,const std::vector<int>& upper,
        const std::vector<int>& year_last,const std::array<double,2>& min_age,
        const std::array<double,2>& max_age,
        const std::vector<std::array<int,2>>& links):parents(links) {
    const int n=parents.size();
    if(lower.size()!=parents.size()||upper.size()!=parents.size()||year_last.size()!=parents.size())
      throw std::invalid_argument("typed age input sizes differ");
    children.resize(n); original_lo.resize(n);original_hi.resize(n);last.resize(n);
    for(int side=0;side<2;++side) {
      if(!std::isfinite(min_age[side])||std::isnan(max_age[side])||min_age[side]>max_age[side]||
         std::abs(min_age[side])>double(infinity/4)||
         (std::isfinite(max_age[side])&&std::abs(max_age[side])>double(infinity/4)))
        throw std::invalid_argument("nonrepresentable typed age limits");
      minimum[side]=Year(std::ceil(min_age[side]));
      maximum[side]=std::isfinite(max_age[side])?Year(std::floor(max_age[side])):infinity;
    }
    for(int v=0;v<n;++v) {
      original_lo[v]=lower[v]!=INT_MIN?lower[v]:-infinity;
      original_hi[v]=upper[v]!=INT_MIN?upper[v]:infinity;
      last[v]=year_last[v]!=INT_MIN?year_last[v]:infinity;
      if(original_lo[v]>original_hi[v])feasible=false;
      for(int p:parents[v])if(p>=0) {
        if(p>=n)throw std::invalid_argument("typed age parent out of range");
        children[p].push_back(v);
      }
    }
    lo=original_lo;hi=original_hi;
    Patch p;p.nodes.resize(n);std::iota(p.nodes.begin(),p.nodes.end(),0);
    solve(p,-1,{{-1,-1}},true);
    feasible=feasible&&p.feasible;if(feasible)publish(p);
  }

  bool exact(int v) const { return original_lo[v]==original_hi[v]; }
  Patch propose(int child,std::array<int,2> replacement) const {
    if(child<0||std::size_t(child)>=parents.size())throw std::invalid_argument("typed age child out of range");
    for(int p:replacement)if(p < -1||p>=int(parents.size()))throw std::invalid_argument("typed age parent out of range");
    Patch patch;
    if(!feasible){patch.feasible=false;return patch;}
    if(parents[child]==replacement)return patch;
    std::vector<unsigned char> seen(parents.size());
    auto add=[&](int v){if(v>=0&&!seen[v]){seen[v]=1;patch.nodes.push_back(v);}};
    add(child);for(int p:parents[child])add(p);for(int p:replacement)add(p);
    // Union of old/new affected components, cut at exact-year boundaries.
    // The endpoints are seeds even when exact; no other branch through an
    // exact endpoint can change when this edge is edited.
    for(std::size_t j=0;j<patch.nodes.size();++j) {
      int v=patch.nodes[j];if(exact(v))continue;
      for(int p:parents[v])add(p);
      for(int c:children[v])add(c);
      if(v==child)for(int p:replacement)add(p);
      if(v==replacement[0]||v==replacement[1])add(child);
    }
    bool removed=false;
    for(int s=0;s<2;++s)removed|=parents[child][s]>=0&&parents[child][s]!=replacement[s];
    solve(patch,child,replacement,removed);
    return patch;
  }

  void commit(int child,std::array<int,2> replacement,const Patch& p) {
    if(!p.feasible)throw std::logic_error("cannot commit infeasible typed ages");
    for(int old:parents[child])if(old>=0) {
      auto& row=children[old];row.erase(std::remove(row.begin(),row.end(),child),row.end());
    }
    parents[child]=replacement;
    for(int parent:replacement)if(parent>=0)children[parent].push_back(child);
    publish(p);
  }
  Bounds bounds(const Patch* patch=nullptr) const {
    Bounds out;out.lo.resize(lo.size());out.hi.resize(hi.size());
    auto integer=[](Year x){return x<=INT_MIN||x>INT_MAX?INT_MIN:int(x);};
    for(std::size_t v=0;v<lo.size();++v){out.lo[v]=integer(lo[v]);out.hi[v]=integer(hi[v]);}
    if(patch)for(std::size_t j=0;j<patch->nodes.size();++j){int v=patch->nodes[j];out.lo[v]=integer(patch->lo[j]);out.hi[v]=integer(patch->hi[j]);}
    return out;
  }
 private:
  void publish(const Patch& p){for(std::size_t j=0;j<p.nodes.size();++j){lo[p.nodes[j]]=p.lo[j];hi[p.nodes[j]]=p.hi[j];}}
  void solve(Patch& out,int child,std::array<int,2> replacement,bool reset) const {
    const int m=out.nodes.size();std::vector<int> index(parents.size(),-1);
    out.lo.resize(m);out.hi.resize(m);
    for(int j=0;j<m;++j){int v=out.nodes[j];index[v]=j;out.lo[j]=reset?original_lo[v]:lo[v];out.hi[j]=reset?original_hi[v]:hi[v];}
    auto links=[&](int v){return v==child?replacement:parents[v];};
    struct Edge {int to;Year weight;};
    std::vector<std::vector<Edge>> edges(m);
    for(int j=0;j<m;++j) {
      int v=out.nodes[j];auto row=links(v);
      for(int s=0;s<2;++s)if(row[s]>=0) {
        int p=row[s],i=index[p];
        out.hi[j]=std::min(out.hi[j],last[p]);
        if(i<0)continue; // unchanged edge beyond an exact boundary
        edges[i].push_back({j,minimum[s]});
        if(maximum[s]!=infinity)edges[j].push_back({i,-maximum[s]});
      }
      if(out.lo[j]>out.hi[j]){out.feasible=false;return;}
    }
    // Max-plus lower closure; reversed edges give min-plus upper closure.
    auto propagate=[&](bool lower) {
      auto& values=lower?out.lo:out.hi;
      std::vector<std::vector<Edge>> reverse;
      if(!lower){reverse.resize(m);for(int i=0;i<m;++i)for(const auto& e:edges[i])reverse[e.to].push_back({i,-e.weight});}
      const auto& graph=lower?edges:reverse;
      std::deque<int> todo;std::vector<unsigned char> queued(m,0);std::vector<int> depth(m);
      auto enqueue=[&](int j){if(j>=0&&!queued[j]){queued[j]=1;todo.push_back(j);}};
      if(reset||child<0)for(int j=0;j<m;++j)enqueue(j);
      else {enqueue(index[child]);for(int p:replacement)if(p>=0)enqueue(index[p]);}
      while(!todo.empty()) {
        int i=todo.front();todo.pop_front();queued[i]=0;
        if(values[i]==(lower?-infinity:infinity))continue;
        for(const auto& e:graph[i]) {
          Year value=values[i]+e.weight;
          if(lower?value<=values[e.to]:value>=values[e.to])continue;
          values[e.to]=value;++out.relaxations;depth[e.to]=depth[i]+1;
          if(out.lo[e.to]>out.hi[e.to]||depth[e.to]>=m){out.feasible=false;return;}
          if(!queued[e.to]){queued[e.to]=1;todo.push_back(e.to);}
        }
      }
    };
    propagate(true);if(out.feasible)propagate(false);
    // An entirely undated component has no finite origin. Check feasibility
    // separately with zero-source distances, without inventing absolute dates.
    if(out.feasible&&m&&std::any_of(out.lo.begin(),out.lo.end(),[](Year x){return x==-infinity;})) {
      std::vector<Year> distance(m,0);
      for(int pass=0;pass<m;++pass){bool changed=false;for(int i=0;i<m;++i)for(const auto& e:edges[i])if(distance[e.to]<distance[i]+e.weight){distance[e.to]=distance[i]+e.weight;changed=true;}
        if(!changed)break;
        if(pass==m-1)out.feasible=false;}
    }
  }
};
}} // namespace
#endif
