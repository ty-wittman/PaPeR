#ifndef PAP_RELATIONSHIP_TOPOLOGY_H
#define PAP_RELATIONSHIP_TOPOLOGY_H

// Symbolic ordinary hypotheses only. No likelihood, R API or inference state.
// The caller owns a common typed observation/retained-edge boundary and must
// provide proposal-aware, tested-edge-masked parent slots before calling this.
#include "score_types.h"

#include "topology_cache.h"
#include "relative_age_cache.h"
#include <chrono>
#include <queue>
#include <functional>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace relationship_topology {

namespace topology_detail {

struct Failure : public std::runtime_error {
  TopologyStatus status;
  Failure(TopologyStatus s, const std::string& reason)
      : std::runtime_error(reason), status(s) {}
};

inline std::string integer_string(int value) {
  std::ostringstream out;
  out << value;
  return out.str();
}

inline const char* side_name(int side) { return side == 0 ? "dam" : "sire"; }

inline int life_value(const std::vector<int>& values, int typed_id) {
  if (typed_id <= 0 || values.empty()) return INT_MIN;
  const std::size_t index = static_cast<std::size_t>(typed_id - 1);
  if (index >= values.size())
    throw Failure(TopologyStatus::Unsupported,
                  "Life-history field does not cover original typed ID");
  return values[index];
}

inline int life_value(const std::vector<int>& owned,
                      const std::vector<int>* borrowed, int typed_id) {
  return life_value(borrowed != nullptr ? *borrowed : owned, typed_id);
}

inline int life_sex_value(const Life& life, int typed_id) {
  // Validate source coverage even when a local inference supplies the value.
  const int value = life_value(life.sex, life.borrowed_sex, typed_id);
  if (typed_id <= 0 || life.sex_overrides.empty()) return value;
  const auto found = std::lower_bound(
    life.sex_overrides.begin(), life.sex_overrides.end(), typed_id,
    [](const std::pair<int, int>& item, int id) { return item.first < id; });
  return found != life.sex_overrides.end() && found->first == typed_id ?
    found->second : value;
}

inline void set_life_sex_override(Life& life, int typed_id, int sex) {
  if (typed_id <= 0) return;
  if (life.sex_overrides.empty() || life.sex_overrides.back().first < typed_id) {
    life.sex_overrides.emplace_back(typed_id, sex);
    return;
  }
  const auto found = std::lower_bound(
    life.sex_overrides.begin(), life.sex_overrides.end(), typed_id,
    [](const std::pair<int, int>& item, int id) { return item.first < id; });
  if (found != life.sex_overrides.end() && found->first == typed_id)
    found->second = sex;
  else life.sex_overrides.insert(found, std::make_pair(typed_id, sex));
}

// Exact role multiplicities after substituting one or two typed birth rows.
// The caller supplies original global counts and the same masked rows used by
// its topology. Missing typed parent slots use zero; no model is evaluated.
inline int role_after_birth_row_substitution(
    std::array<int, 2> counts, int typed_id,
    const std::array<std::array<int, 2>, 2>& original_rows,
    const std::array<std::array<int, 2>, 2>& masked_rows,
    int changed_rows) {
  for (int row = 0; row < changed_rows; ++row)
    for (int side = 0; side < 2; ++side)
      counts[static_cast<std::size_t>(side)] +=
        (masked_rows[static_cast<std::size_t>(row)][static_cast<std::size_t>(side)] == typed_id ? 1 : 0) -
        (original_rows[static_cast<std::size_t>(row)][static_cast<std::size_t>(side)] == typed_id ? 1 : 0);
  return (counts[0] > 0 ? 1 : 0) | (counts[1] > 0 ? 2 : 0);
}

inline Life materialize_life(const Life& life) {
  Life result;
  result.sex = life.borrowed_sex != nullptr ? *life.borrowed_sex : life.sex;
  result.birth_lo = life.borrowed_birth_lo != nullptr ?
    *life.borrowed_birth_lo : life.birth_lo;
  result.birth_hi = life.borrowed_birth_hi != nullptr ?
    *life.borrowed_birth_hi : life.birth_hi;
  result.year_last = life.borrowed_year_last != nullptr ?
    *life.borrowed_year_last : life.year_last;
  result.sex_overrides = life.sex_overrides;
  result.min_age = life.min_age;
  result.max_age = life.max_age;
  return result;
}

class Builder {
 public:
  explicit Builder(const Input& input) : input_(input) {
    validate_input();
    nodes_ = input_.nodes;
    if (input_.optional_source_birth_budget >= 0) required_birth_.resize(nodes_.size(), false);
    representative_.resize(nodes_.size());
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
      representative_[i] = static_cast<int>(i);
      nodes_[i].source_input_index = static_cast<int>(i);
      if (nodes_[i].label.empty())
        nodes_[i].label = nodes_[i].typed_id > 0 ?
          "T" + integer_string(nodes_[i].typed_id) :
          "L_input_" + integer_string(static_cast<int>(i));
      if (nodes_[i].parents[0] >= 0 || nodes_[i].parents[1] >= 0)
        nodes_[i].has_birth = true;
    }
    // Finish all declared births before merging repeated typed roles.
    const std::size_t original_size = nodes_.size();
    for (std::size_t i = 0; i < original_size; ++i)
      if (nodes_[i].has_birth) ensure_birth(static_cast<int>(i));
    ensure_birth(input_.focal);
    ensure_birth(input_.candidate);
    if (!required_birth_.empty()) required_birth_[static_cast<std::size_t>(find(input_.focal))] = true;
    std::map<int, int> typed_ids;
    for (std::size_t i = 0; i < original_size; ++i) {
      if (nodes_[i].typed_id <= 0) continue;
      const int id = nodes_[i].typed_id;
      const std::map<int, int>::iterator found = typed_ids.find(id);
      if (found == typed_ids.end()) typed_ids[id] = static_cast<int>(i);
      else unify(found->second, static_cast<int>(i));
    }
    if (find(input_.focal) == find(input_.candidate))
      inapplicable("Focal and candidate are the same canonical typed identity");
  }

  bool bound_to(const Input& input)const{return &input_==&input;}

  // Narrow symbolic construction seam for separately declared compound DAGs.
  // These wrappers retain the existing slot completion, typed unification,
  // resource and joint applicability rules; ordinary apply() is unchanged.
  int parent_slot(int node, int side) {
    if (node < 0 || static_cast<std::size_t>(node) >= nodes_.size() ||
        (side != 1 && side != 2))
      throw Failure(TopologyStatus::Unsupported,
                    "Compound parent-slot request outside builder boundary");
    return slots(node)[static_cast<std::size_t>(side - 1)];
  }

  int unify_nodes(int first, int second) {
    if (first < 0 || second < 0 ||
        static_cast<std::size_t>(first) >= nodes_.size() ||
        static_cast<std::size_t>(second) >= nodes_.size())
      throw Failure(TopologyStatus::Unsupported,
                    "Compound unification request outside builder boundary");
    return unify(first, second);
  }

  void require_distinct_nodes(int first, int second, const std::string& reason) {
    if (first < 0 || second < 0 ||
        static_cast<std::size_t>(first) >= nodes_.size() ||
        static_cast<std::size_t>(second) >= nodes_.size())
      throw Failure(TopologyStatus::Unsupported,
                    "Compound distinct-role request outside builder boundary");
    require_distinct(first, second, reason);
  }

  void apply(const TopologyModel& model) {

    const int focal = find(input_.focal);
    const int candidate = find(input_.candidate);
    if (model.relationship == "PO") {
      if (input_.proposed_parent_side != 1 && input_.proposed_parent_side != 2)
        inapplicable("Proposal PO requires a tested dam or sire side");
      nodes_[static_cast<std::size_t>(focal)].parents[
        static_cast<std::size_t>(input_.proposed_parent_side - 1)] = candidate;
      provenance_.push_back("Supply tested direct PO edge");
    } else if (model.relationship == "U") {
      provenance_.push_back("No added relationship; retain common masked birth context");
    } else if (model.relationship == "FS") {
      const std::array<int, 2> first = slots(focal);
      const std::array<int, 2> second = slots(candidate);
      unify(first[0], second[0]);
      unify(first[1], second[1]);
    } else if (model.relationship == "HS") {
      const int shared = model.shared_side - 1;
      const std::array<int, 2> first = slots(focal);
      const std::array<int, 2> second = slots(candidate);
      unify(first[static_cast<std::size_t>(shared)],
            second[static_cast<std::size_t>(shared)]);
      require_distinct(first[static_cast<std::size_t>(1 - shared)],
                       second[static_cast<std::size_t>(1 - shared)],
                       "HS opposite parent identities coincide");
    } else {
      const int descendant = model.focal_is_descendant ? focal : candidate;
      const int relative = model.focal_is_descendant ? candidate : focal;
      const int bridge = slots(descendant)[static_cast<std::size_t>(model.bridge_side - 1)];
      require_distinct(bridge, relative, "Bridge/intermediate is the relative itself");
      const std::array<int, 2> bridge_parents = slots(bridge);
      if (model.relationship == "GP") {
        unify(bridge_parents[static_cast<std::size_t>(model.shared_side - 1)], relative);
      } else {
        const std::array<int, 2> relative_parents = slots(relative);
        if (model.relationship == "FA") {
          unify(bridge_parents[0], relative_parents[0]);
          unify(bridge_parents[1], relative_parents[1]);
        } else if (model.relationship == "HA") {
          const int shared = model.shared_side - 1;
          unify(bridge_parents[static_cast<std::size_t>(shared)],
                relative_parents[static_cast<std::size_t>(shared)]);
          require_distinct(bridge_parents[static_cast<std::size_t>(1 - shared)],
                           relative_parents[static_cast<std::size_t>(1 - shared)],
                           "HA opposite grandparent identities coincide");
        } else {
          throw Failure(TopologyStatus::Unsupported, "Unknown ordinary relationship");
        }
      }
    }
  }

  // Exact current boundary identity is checked even after eager invalidation:
  // newly introduced relatives cannot be certified from old dependencies.
  static std::string bound_fingerprint(const Input& input){
    using pap_internal::topology_templates::append;using pap_internal::topology_templates::append_string;std::string key;
    append(key,input.focal);append(key,input.candidate);append(key,input.proposed_parent_side);append(key,input.max_nodes);append(key,input.max_typed_nodes);append(key,input.optional_source_birth_budget);
    append(key,input.nodes.size());for(const auto& n:input.nodes){append(key,n.typed_id);append(key,n.observed);append(key,n.neutral_phi);append(key,n.has_birth);append(key,n.parents);append(key,n.source_input_index);append_string(key,n.label);}
    append(key,input.retained_source_birth_ids.size());for(int id:input.retained_source_birth_ids)append(key,id);
    append(key,input.structural_source_parents.size());for(const auto& p:input.structural_source_parents)append(key,p);return key;
  }
  static std::string bound_key(const Input& input,const TopologyModel& model,const std::string& tag){
    using pap_internal::topology_templates::append;using pap_internal::topology_templates::append_string;std::string key;
    if(input.focal<0||input.candidate<0||size_t(input.focal)>=input.nodes.size()||size_t(input.candidate)>=input.nodes.size())return key;
    append(key,input.nodes[input.focal].typed_id);append(key,input.nodes[input.candidate].typed_id);append(key,input.proposed_parent_side);
    append_string(key,model.name);append_string(key,model.relationship);append(key,model.model_order);append(key,model.orientation_index);append(key,model.focal_is_descendant);append(key,model.bridge_side);append(key,model.shared_side);append_string(key,tag);return key;
  }
  static bool ancestry_certificate(const TopologyModel& model){
    if(model.topological_order.size()!=model.nodes.size())return false;std::vector<int> rank(model.nodes.size(),-1);
    for(size_t i=0;i<model.topological_order.size();++i){int v=model.topological_order[i];if(v<0||size_t(v)>=rank.size()||rank[v]>=0)return false;rank[v]=i;}
    for(size_t v=0;v<model.nodes.size();++v)if(model.nodes[v].has_birth)for(int p:model.nodes[v].parents)if(p<0||size_t(p)>=rank.size()||rank[p]>=rank[v])return false;
    return true;
  }
  static void validate_life_current(const TopologyModel& model,const Life& life){

    bool done=false;try{done=verify_relative_life(model,life);}catch(const std::bad_alloc&){}catch(const std::length_error&){}
    if(done)return;auto* cache=pap_internal::topology_life_reuse::active;
    auto original=[&]{pap_internal::topology_life_reuse::Scope uncached(nullptr);verify_life(model,life);};
    if(!cache)verify_life(model,life);else if(!cache->limit){++cache->refusals;original();}else try{verify_life_components(model,life);}catch(const std::bad_alloc&){++cache->refusals;original();}catch(const std::length_error&){++cache->refusals;original();}
  }
  static std::string transfer_input(const Input& input){
    auto key=bound_fingerprint(input);using pap_internal::topology_templates::append;
    append(key,bool(input.life));if(input.life){const auto& l=*input.life;append(key,l.min_age);append(key,l.max_age);
      for(const auto& n:input.nodes)if(n.typed_id>0){append(key,life_sex_value(l,n.typed_id));append(key,life_value(l.birth_lo,l.borrowed_birth_lo,n.typed_id));append(key,life_value(l.birth_hi,l.borrowed_birth_hi,n.typed_id));append(key,life_value(l.year_last,l.borrowed_year_last,n.typed_id));}}
    return key;
  }
  static bool transferred(const Input& input,TopologyModel& model,const std::string& tag){
    auto* c=FinishedModelTransfer::active;if(!c||!c->reader)return false;++c->requests;
    try{auto started=std::chrono::steady_clock::now();auto key=bound_key(input,model,tag);auto at=c->entries.find(key);
      if(at==c->entries.end()){++c->missing;return false;}
      auto source=transfer_input(input);c->key_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
      if(at->second.input!=source){++c->mismatches;return false;}
      started=std::chrono::steady_clock::now();model=at->second.model;c->copy_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();++c->hits;return true;
    }catch(const std::bad_alloc&){++c->refusals;return false;}catch(const std::length_error&){++c->refusals;return false;}
  }
  static bool try_bound(const Input& input,TopologyModel& model,const std::string& tag={}){
    if(transferred(input,model,tag))return true;
    auto* cache=TopologyTemplateCache::active;if(!cache||!cache->bound_enabled)return false;++cache->bound_requests;
    try{auto started=std::chrono::steady_clock::now();auto key=bound_key(input,model,tag);auto at=cache->bound.find(key);if(at==cache->bound.end())return false;
      auto fingerprint=bound_fingerprint(input);if(fingerprint!=at->second.fingerprint){cache->bytes-=at->second.charge;cache->bound.erase(at);++cache->closure_misses;return false;}
      cache->lookup_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();

      auto rebound=at->second.model;
      for(const auto& birth:rebound.structural_births)if(birth.cut){rebound.nodes[birth.child].has_birth=true;rebound.nodes[birth.child].parents=birth.parents;}
      {if(!ancestry_certificate(rebound))return false;}
      if(input.life)validate_life_current(rebound,*input.life);
      for(const auto& birth:rebound.structural_births)if(birth.cut){rebound.nodes[birth.child].has_birth=false;rebound.nodes[birth.child].parents={{-1,-1}};}
      model=std::move(rebound);++cache->bound_hits;return true;
    }catch(const Failure&){return false;}catch(const std::bad_alloc&){++cache->refusals;return false;}catch(const std::length_error&){++cache->refusals;return false;}
  }
  static void remember_bound(const Input& input,const TopologyModel& model,const std::string& tag={}){
    if(auto* c=FinishedModelTransfer::active;c&&!c->reader){++c->offers;
      try{auto started=std::chrono::steady_clock::now();auto key=bound_key(input,model,tag);if(!key.empty()){
        auto fingerprint=transfer_input(input);c->key_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        // One current engine result per biological model identity; exact input
        // equality is still required at consumption. Refusal uses the builder.
        size_t charge=sizeof(FinishedModelTransfer::Entry)+key.capacity()+fingerprint.capacity()+model.nodes.size()*sizeof(TopologyNode)+model.structural_births.size()*sizeof(StructuralBirth)+model.original_to_canonical.size()*sizeof(int)+model.topological_order.size()*sizeof(int)+1024;
        if(c->bytes+charge<=16U*1024U*1024U){c->entries.insert_or_assign(std::move(key),FinishedModelTransfer::Entry{std::move(fingerprint),model});c->bytes+=charge;c->peak_bytes=std::max(c->peak_bytes,c->bytes);}else ++c->refusals;
      }}catch(const std::bad_alloc&){++c->refusals;}catch(const std::length_error&){++c->refusals;}
    }
    auto* cache=TopologyTemplateCache::active;if(!cache||!cache->bound_enabled||!model.ready())return;
    try{auto key=bound_key(input,model,tag);if(key.empty()||cache->bound.count(key))return;typename TopologyTemplateCache::Bound entry;entry.model=model;entry.fingerprint=bound_fingerprint(input);
      for(const auto& n:input.nodes)if(n.typed_id>0)entry.dependencies.push_back(n.typed_id);std::sort(entry.dependencies.begin(),entry.dependencies.end());entry.dependencies.erase(std::unique(entry.dependencies.begin(),entry.dependencies.end()),entry.dependencies.end());
      size_t charge=sizeof(entry)+key.capacity()+entry.fingerprint.capacity()+entry.dependencies.capacity()*sizeof(int)+model.nodes.capacity()*sizeof(TopologyNode)+model.structural_births.capacity()*sizeof(StructuralBirth)+(model.original_to_canonical.capacity()+model.topological_order.capacity())*sizeof(int)+256;
      for(const auto& n:model.nodes)charge+=n.label.capacity();for(const auto& p:model.provenance)charge+=sizeof(std::string)+p.capacity();entry.charge=charge;
      if(cache->room(charge)){cache->bound.emplace(std::move(key),std::move(entry));cache->bytes+=charge;}
    }catch(const std::bad_alloc&){++cache->refusals;}catch(const std::length_error&){++cache->refusals;}
  }

  void finish(TopologyModel& model) {



    std::vector<int> roots;auto* templates=TopologyTemplateCache::active;
    std::string geometry_key;std::shared_ptr<const TopologyTemplateCache::Geometry> geometry;
    if(templates&&templates->shape_enabled&&model.nodes.empty())try{++templates->shape_requests;geometry_key=shape_key();auto found=templates->shapes.find(geometry_key);if(found!=templates->shapes.end()){geometry=found->second;++templates->shape_hits;}}catch(const std::bad_alloc&){++templates->refusals;geometry_key.clear();}catch(const std::length_error&){++templates->refusals;geometry_key.clear();}
    if(geometry)try{
      roots=geometry->roots;for(size_t i=0;i<roots.size();++i){auto node=nodes_[roots[i]];node.parents=geometry->parents[i];model.nodes.push_back(std::move(node));}
      model.focal=geometry->focal;model.candidate=geometry->candidate;model.original_to_canonical=geometry->original;model.topological_order=geometry->order;

      if(!ancestry_certificate(model)){++templates->refusals;geometry.reset();geometry_key.clear();roots.clear();model.nodes.clear();model.original_to_canonical.clear();model.topological_order.clear();}
    }catch(const std::bad_alloc&){++templates->refusals;geometry.reset();geometry_key.clear();roots.clear();model.nodes.clear();model.original_to_canonical.clear();model.topological_order.clear();}
    catch(const std::length_error&){++templates->refusals;geometry.reset();geometry_key.clear();roots.clear();model.nodes.clear();model.original_to_canonical.clear();model.topological_order.clear();}
    if(!geometry){
    // Keep every original typed identity, including missing observations.
    // Only latent roots disconnected from all retained typed nodes are pruned.
    const bool compact_projection = pap_internal::topology_life_reuse::active != nullptr;
    std::set<int> retained;
    std::vector<unsigned char> retained_flags(compact_projection ? nodes_.size() : 0);
    std::vector<int> pending;
    {
    for (std::size_t i = 0; i < input_.nodes.size(); ++i) {
      const int root = find(static_cast<int>(i));
      if (nodes_[static_cast<std::size_t>(root)].typed_id > 0 ||
          nodes_[static_cast<std::size_t>(root)].observed)
        pending.push_back(root);
    }
    while (!pending.empty()) {
      const int root = find(pending.back());
      pending.pop_back();
      if (compact_projection) {
        if (retained_flags[root]) continue;
        retained_flags[root] = 1;
      } else if (!retained.insert(root).second) continue;
      const TopologyNode& node = nodes_[static_cast<std::size_t>(root)];
      if (node.has_birth) {
        pending.push_back(find(node.parents[0]));
        pending.push_back(find(node.parents[1]));
      }
    }
    if (compact_projection) {
      for (size_t v = 0; v < retained_flags.size(); ++v) if (retained_flags[v]) roots.push_back(v);
    } else roots.assign(retained.begin(), retained.end());
    }
    std::sort(roots.begin(), roots.end(), RootOrder(nodes_));
    std::map<int, int> compact;
    std::vector<int> compact_index(compact_projection ? nodes_.size() : 0, -1);
    for (std::size_t i = 0; i < roots.size(); ++i) {
      if (compact_projection) compact_index[roots[i]] = i;
      else compact[roots[i]] = static_cast<int>(i);
    }
    auto canonical = [&](int root) { if (compact_projection && compact_index.at(root) >= 0) return compact_index[root]; return compact.at(root); };
    for (std::size_t i = 0; i < roots.size(); ++i) {
      TopologyNode node = nodes_[static_cast<std::size_t>(roots[i])];
      if (node.has_birth) {
        node.parents[0] = canonical(find(node.parents[0]));
        node.parents[1] = canonical(find(node.parents[1]));
      } else node.parents = std::array<int, 2>{{-1, -1}};
      model.nodes.push_back(node);
    }
    model.focal = canonical(find(input_.focal));
    model.candidate = canonical(find(input_.candidate));
    model.original_to_canonical.resize(input_.nodes.size(), -1);
    for (std::size_t i = 0; i < input_.nodes.size(); ++i) {
      if (compact_projection) model.original_to_canonical[i] = compact_index[find(static_cast<int>(i))];
      else {
        const std::map<int, int>::const_iterator found = compact.find(find(static_cast<int>(i)));
        if (found != compact.end()) model.original_to_canonical[i] = found->second;
      }
    }
    {verify_topological_order(model);}
    if(templates&&!geometry_key.empty())try{
      auto value=std::make_shared<TopologyTemplateCache::Geometry>();value->roots=roots;value->original=model.original_to_canonical;value->order=model.topological_order;value->focal=model.focal;value->candidate=model.candidate;for(const auto& node:model.nodes)value->parents.push_back(node.parents);
      size_t charge=sizeof(*value)+geometry_key.capacity()+(value->roots.capacity()+value->original.capacity()+value->order.capacity())*sizeof(int)+value->parents.capacity()*sizeof(value->parents[0])+128;
      if(templates->room(charge)){templates->shapes.emplace(std::move(geometry_key),std::move(value));templates->bytes+=charge;++templates->shape_builds;}
    }catch(const std::bad_alloc&){++templates->refusals;}catch(const std::length_error&){++templates->refusals;}
    }
    if(input_.life)validate_life_current(model,*input_.life);
    if (input_.optional_source_birth_budget >= 0) {
      // The complete local motif was validated above. Its structural identities
      // remain recorded; only explicitly optional numerical birth factors cut.
      for (std::size_t i = 0; i < roots.size(); ++i) {
        auto& node = model.nodes[i];
        if (!node.has_birth) continue;
        StructuralBirth row;
        row.child = static_cast<int>(i); row.parents = node.parents;
        row.required_by_motif = required_birth_[static_cast<std::size_t>(find(roots[i]))];
        row.protected_source = std::binary_search(input_.retained_source_birth_ids.begin(),
          input_.retained_source_birth_ids.end(), node.typed_id);
        row.cut = !row.required_by_motif && !row.protected_source;
        model.structural_births.push_back(row);
        if (row.cut) {
          node.has_birth = false; node.parents = std::array<int, 2>{{-1, -1}};
        }
      }
      provenance_.push_back("Optional source births cut at common declared B; cut child receives population prior");
    }
    model.provenance = provenance_;
    model.provenance.push_back("One canonical node and observation owner per retained typed ID");
    model.provenance.push_back("One birth or founder prior per node; missing slots are distinct unless motif-unified");
    model.status = TopologyStatus::Ready;
    model.reason = "Canonical ordinary topology and joint applicability valid; likelihood unscored";
  }

 private:
  std::string shape_key(){
    using pap_internal::topology_templates::append;std::string key;
    append(key,input_.nodes.size());append(key,input_.focal);append(key,input_.candidate);append(key,nodes_.size());
    std::vector<int> order(nodes_.size()),rank(nodes_.size());for(size_t v=0;v<nodes_.size();++v)order[v]=v;std::sort(order.begin(),order.end(),RootOrder(nodes_));for(size_t i=0;i<order.size();++i)rank[order[i]]=i;
    for(size_t v=0;v<nodes_.size();++v){const auto& n=nodes_[v];append(key,find(v));append(key,rank[v]);bool typed=n.typed_id>0;append(key,typed);append(key,n.observed);append(key,n.neutral_phi);append(key,n.has_birth);append(key,n.parents);}
    return key;
  }

  const Input& input_;
  std::vector<TopologyNode> nodes_;
  std::vector<int> representative_;
  std::vector<bool> required_birth_;
  std::vector<std::string> provenance_;

  struct RootOrder {
    const std::vector<TopologyNode>& nodes;
    explicit RootOrder(const std::vector<TopologyNode>& n) : nodes(n) {}
    bool operator()(int a, int b) const {
      const TopologyNode& first = nodes[static_cast<std::size_t>(a)];
      const TopologyNode& second = nodes[static_cast<std::size_t>(b)];
      if ((first.typed_id > 0) != (second.typed_id > 0)) return first.typed_id > 0;
      if (first.typed_id > 0 && first.typed_id != second.typed_id)
        return first.typed_id < second.typed_id;
      if (first.label != second.label) return first.label < second.label;
      return a < b;
    }
  };

  static void inapplicable(const std::string& reason) {
    throw Failure(TopologyStatus::Inapplicable, reason);
  }

  void validate_input() const {
    if (input_.nodes.empty() || input_.max_nodes == 0U || input_.max_typed_nodes == 0U)
      throw Failure(TopologyStatus::Unsupported, "Empty boundary or zero topology resource cap");
    if (input_.nodes.size() > input_.max_nodes)
      throw Failure(TopologyStatus::Unsupported, "Input exceeds topology node cap");
    if (input_.focal < 0 || input_.candidate < 0 ||
        static_cast<std::size_t>(input_.focal) >= input_.nodes.size() ||
        static_cast<std::size_t>(input_.candidate) >= input_.nodes.size())
      throw Failure(TopologyStatus::Unsupported, "Endpoint index outside common boundary");
    if (input_.nodes[static_cast<std::size_t>(input_.focal)].typed_id <= 0 ||
        input_.nodes[static_cast<std::size_t>(input_.candidate)].typed_id <= 0)
      throw Failure(TopologyStatus::Unsupported, "Ordinary endpoints must be typed identities");
    if (input_.proposed_parent_side < 0 || input_.proposed_parent_side > 2)
      throw Failure(TopologyStatus::Unsupported, "Invalid tested parental side");
    if (input_.optional_source_birth_budget < -1 ||
        !std::is_sorted(input_.retained_source_birth_ids.begin(), input_.retained_source_birth_ids.end()))
      throw Failure(TopologyStatus::Unsupported, "Invalid common source birth budget/ledger");
    std::set<int> typed;
    for (std::size_t i = 0; i < input_.nodes.size(); ++i) {
      const TopologyNode& node = input_.nodes[i];
      if (node.typed_id < 0 || (node.typed_id == 0 && node.observed))
        throw Failure(TopologyStatus::Unsupported, "Invalid typed/latent observation identity");
      if (node.typed_id > 0) typed.insert(node.typed_id);
      for (int side = 0; side < 2; ++side) {
        const int parent = node.parents[static_cast<std::size_t>(side)];
        if (parent < -1 || (parent >= 0 && static_cast<std::size_t>(parent) >= input_.nodes.size()))
          throw Failure(TopologyStatus::Unsupported, "Parent index outside common boundary");
      }
    }
    if (typed.size() > input_.max_typed_nodes)
      throw Failure(TopologyStatus::Unsupported, "Input exceeds typed-node cap");
  }

  int find(int node) {
    const int parent = representative_[static_cast<std::size_t>(node)];
    if (parent != node)
      representative_[static_cast<std::size_t>(node)] = find(parent);
    return representative_[static_cast<std::size_t>(node)];
  }

  int new_missing(int owner, int side) {
    if (nodes_.size() >= input_.max_nodes)
      throw Failure(TopologyStatus::Unsupported, "Hypothesis exceeds typed-plus-latent node cap");
    TopologyNode latent;
    latent.label = "M(" + nodes_[static_cast<std::size_t>(find(owner))].label +
      "," + side_name(side) + ")";
    if (find(owner) == find(input_.focal) && input_.proposed_parent_side == side + 1)
      latent.label = "Z(" + nodes_[static_cast<std::size_t>(find(owner))].label +
        "," + side_name(side) + ")";
    const int index = static_cast<int>(nodes_.size());
    nodes_.push_back(latent);
    representative_.push_back(index);
    if (input_.optional_source_birth_budget >= 0) required_birth_.push_back(false);
    return index;
  }

  void ensure_birth(int node) {
    const int root = find(node);
    const bool previous_birth = nodes_[static_cast<std::size_t>(root)].has_birth;
    nodes_[static_cast<std::size_t>(root)].has_birth = true;
    for (int side = 0; side < 2; ++side) {
      // No reference into nodes_ may survive new_missing(), which can reallocate.
      if (nodes_[static_cast<std::size_t>(root)].parents[static_cast<std::size_t>(side)] < 0) {
        const int missing = new_missing(root, side);
        nodes_[static_cast<std::size_t>(root)].parents[static_cast<std::size_t>(side)] = missing;
      }
    }
    if (!previous_birth)
      provenance_.push_back("Materialize open birth for " + nodes_[static_cast<std::size_t>(root)].label);
  }

  std::array<int, 2> slots(int node) {
    const int root = find(node);
    ensure_birth(root);
    if (!required_birth_.empty()) required_birth_[static_cast<std::size_t>(root)] = true;
    const std::array<int, 2> pair = nodes_[static_cast<std::size_t>(root)].parents;
    return std::array<int, 2>{{find(pair[0]), find(pair[1])}};
  }

  void require_distinct(int first, int second, const std::string& reason) {
    if (find(first) == find(second)) inapplicable(reason);
  }

  int unify(int first, int second) {
    int a = find(first);
    int b = find(second);
    if (a == b) return a;
    const int first_id = nodes_[static_cast<std::size_t>(a)].typed_id;
    const int second_id = nodes_[static_cast<std::size_t>(b)].typed_id;
    if (first_id > 0 && second_id > 0 && first_id != second_id)
      inapplicable("Conflicting retained typed identities " + integer_string(first_id) +
                   " and " + integer_string(second_id));
    if ((second_id > 0 && first_id == 0) ||
        ((first_id > 0) == (second_id > 0) && RootOrder(nodes_)(b, a)))
      std::swap(a, b);
    const bool a_birth = nodes_[static_cast<std::size_t>(a)].has_birth;
    const bool b_birth = nodes_[static_cast<std::size_t>(b)].has_birth;
    const std::array<int, 2> a_parents = nodes_[static_cast<std::size_t>(a)].parents;
    const std::array<int, 2> b_parents = nodes_[static_cast<std::size_t>(b)].parents;
    const bool a_observed = nodes_[static_cast<std::size_t>(a)].observed;
    const bool b_observed = nodes_[static_cast<std::size_t>(b)].observed;
    // A real observed owner dominates a unit owner. An unobserved latent's
    // default policy must never overwrite its observed typed representative.
    nodes_[static_cast<std::size_t>(a)].neutral_phi = a_observed ?
      (nodes_[static_cast<std::size_t>(a)].neutral_phi &&
       (!b_observed || nodes_[static_cast<std::size_t>(b)].neutral_phi)) :
      (b_observed && nodes_[static_cast<std::size_t>(b)].neutral_phi);
    nodes_[static_cast<std::size_t>(a)].observed = a_observed || b_observed;
    if (!required_birth_.empty()) required_birth_[static_cast<std::size_t>(a)] =
      required_birth_[static_cast<std::size_t>(a)] || required_birth_[static_cast<std::size_t>(b)];
    representative_[static_cast<std::size_t>(b)] = a;
    if (!a_birth && b_birth) {
      nodes_[static_cast<std::size_t>(a)].has_birth = true;
      nodes_[static_cast<std::size_t>(a)].parents = b_parents;
    } else if (a_birth && b_birth) {
      unify(a_parents[0], b_parents[0]);
      unify(a_parents[1], b_parents[1]);
      const int root = find(a);
      nodes_[static_cast<std::size_t>(root)].has_birth = true;
      nodes_[static_cast<std::size_t>(root)].parents = a_parents;
    }
    return find(a);
  }

  static void verify_topological_order(TopologyModel& model) {
    const std::size_t size = model.nodes.size();
    if(auto* cache=pap_internal::topology_life_reuse::active;cache&&cache->flat_topology){
      // Identical Kahn traversal and minimum-index order, with one flat edge
      // buffer instead of a separate allocation for every parent's children.
      std::vector<int> degree(size),offset(size+1);
      for(size_t child=0;child<size;++child)if(model.nodes[child].has_birth)for(int p:model.nodes[child].parents){if(p<0||size_t(p)>=size)throw Failure(TopologyStatus::Unsupported,"Unresolved canonical birth slot");++degree[child];++offset[p+1];}
      for(size_t v=1;v<offset.size();++v)offset[v]+=offset[v-1];
      auto cursor=offset;std::vector<int> edges(offset.back());
      for(size_t child=0;child<size;++child)if(model.nodes[child].has_birth)for(int p:model.nodes[child].parents)edges[cursor[p]++]=child;
      std::priority_queue<int,std::vector<int>,std::greater<int>> ready;
      for(size_t v=0;v<size;++v)if(!degree[v])ready.push(v);
      while(!ready.empty()){int p=ready.top();ready.pop();model.topological_order.push_back(p);for(int e=offset[p];e<offset[p+1];++e)if(--degree[edges[e]]==0)ready.push(edges[e]);}
      if(model.topological_order.size()!=size)inapplicable("Directed ancestry cycle after ordinary slot unification");return;
    }
    std::vector<int> degree(size, 0);
    std::vector<std::vector<int> > children(size);
    std::set<int> ready;
    const bool compact_order = pap_internal::topology_life_reuse::active != nullptr;
    std::priority_queue<int, std::vector<int>, std::greater<int>> ready_heap;
    auto insert_ready = [&](int v) { if (compact_order) ready_heap.push(v); else ready.insert(v); };
    for (std::size_t child = 0; child < size; ++child) {
      const TopologyNode& node = model.nodes[child];
      if (!node.has_birth) insert_ready(static_cast<int>(child));
      else for (int side = 0; side < 2; ++side) {
        const int parent = node.parents[static_cast<std::size_t>(side)];
        if (parent < 0 || static_cast<std::size_t>(parent) >= size)
          throw Failure(TopologyStatus::Unsupported, "Unresolved canonical birth slot");
        ++degree[child];
        children[static_cast<std::size_t>(parent)].push_back(static_cast<int>(child));
      }
    }
    while (compact_order ? !ready_heap.empty() : !ready.empty()) {
      const int parent = compact_order ? ready_heap.top() : *ready.begin();
      if (compact_order) ready_heap.pop(); else ready.erase(ready.begin());
      model.topological_order.push_back(parent);
      const std::vector<int>& direct = children[static_cast<std::size_t>(parent)];
      for (std::size_t i = 0; i < direct.size(); ++i) {
        const int child = direct[i];
        if (--degree[static_cast<std::size_t>(child)] == 0) insert_ready(child);
      }
    }
    if (model.topological_order.size() != size)
      inapplicable("Directed ancestry cycle after ordinary slot unification");
  }

  struct DifferenceConstraint {
    int from;
    int to;
    long long upper;
    std::string reason;
    DifferenceConstraint(int f, int t, long long u, const std::string& r)
        : from(f), to(t), upper(u), reason(r) {}
  };

  static long long integer_age(double value, bool lower) {
    const double rounded = lower ? std::ceil(value) : std::floor(value);
    const double safe = static_cast<double>(std::numeric_limits<long long>::max() / 2048);
    if (!std::isfinite(rounded) || rounded < -safe || rounded > safe)
      throw Failure(TopologyStatus::Unsupported, "Nonrepresentable integer age bound");
    return static_cast<long long>(rounded);
  }

  static bool safe_difference(long long a, long long b, long long& result) {
    if ((b > 0 && a < std::numeric_limits<long long>::min() + b) ||
        (b < 0 && a > std::numeric_limits<long long>::max() + b)) return false;
    result = a - b;
    return true;
  }

  static bool check_life_witness(const TopologyModel& model, const Life& life) {
    auto& cache = *pap_internal::topology_life_reuse::active;
    if (!cache.limit) return false;
    ++cache.witness_checks;
    std::vector<long long> birth(model.nodes.size());
    std::vector<int> roles(model.nodes.size());
    for (size_t v = 0; v < model.nodes.size(); ++v) {
      const auto& node = model.nodes[v];
      const auto old = cache.witness.find(node.label);
      long long value = old == cache.witness.end() ? 0 : old->second;
      const int lo = life_value(life.birth_lo, life.borrowed_birth_lo, node.typed_id);
      const int hi = life_value(life.birth_hi, life.borrowed_birth_hi, node.typed_id);
      if (lo != INT_MIN) value = std::max(value, static_cast<long long>(lo));
      if (hi != INT_MIN) value = std::min(value, static_cast<long long>(hi));
      if ((lo != INT_MIN && value < lo) || (hi != INT_MIN && value > hi)) return false;
      birth[v] = value;
    }
    for (size_t v = 0; v < model.nodes.size(); ++v) {
      if (!model.nodes[v].has_birth) continue;
      for (int side = 0; side < 2; ++side) {
        const int p = model.nodes[v].parents[side];
        roles[p] |= 1 << side;
        const double lo = life.min_age[side], hi = life.max_age[side];
        if (std::isnan(lo) || std::isnan(hi) || hi < lo) return false;
        long long age;
        if (!safe_difference(birth[v], birth[p], age) || age < integer_age(lo, true)) return false;
        if (std::isfinite(hi)) {
          if (age > integer_age(hi, false)) return false;
        } else if (hi < 0) return false;
        const int last = life_value(life.year_last, life.borrowed_year_last,
                                    model.nodes[p].typed_id);
        if (last != INT_MIN && birth[v] > last) return false;
      }
    }
    for (size_t v = 0; v < model.nodes.size(); ++v) {
      const int sex = life_sex_value(life, model.nodes[v].typed_id);
      if (roles[v] == 3 || (sex == 1 && (roles[v] & 2)) ||
          (sex == 2 && (roles[v] & 1))) return false;
    }
    // Labels propose birthdays, never identify a scientific result. Every
    // interval, NEW edge, parent orientation and Year.last was checked above.
    ++cache.witness_hits;
    return true;
  }

  static void remember_life_witness(const TopologyModel& model,
      const std::vector<long long>& distances, int zero) {
    auto* cache = pap_internal::topology_life_reuse::active;
    if (!cache || !cache->limit) return;
    for (size_t v = 0; v < model.nodes.size(); ++v) {
      long long birthday;
      if (!safe_difference(distances[v], distances[zero], birthday)) continue;
      const auto& label = model.nodes[v].label;
      const auto at = cache->witness.find(label);
      if (at != cache->witness.end()) {
        at->second = birthday;
        continue;
      }
      const size_t charge = sizeof(std::pair<const std::string, long long>) +
                            label.size() + 1 + 4 * sizeof(void*);
      if (charge <= cache->limit - cache->bytes) {
        try {
          cache->witness.emplace(label, birthday);
          cache->bytes += charge;
        } catch (const std::bad_alloc&) { ++cache->refusals; return; }
          catch (const std::length_error&) { ++cache->refusals; return; }
      }
    }
  }

  static void verify_life_components(const TopologyModel& model, const Life& life) {
    auto& cache = *pap_internal::topology_life_reuse::active;
    using Clock = std::chrono::steady_clock;
    try {
      if (check_life_witness(model, life)) return;
    } catch (const Failure&) {
      verify_life(model, life);
      return;
    }
    const size_t n = model.nodes.size();
    std::map<int, std::vector<int>> components;
    {
    std::vector<int> parent(n);
    for (size_t v = 0; v < n; ++v) parent[v] = v;
    auto root = [&](int v) {
      while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; }
      return v;
    };
    // Recompute components on the NEW motif graph, before optional birth cuts.
    // Newly attached ancestors or co-parents necessarily alter this closure.
    for (size_t v = 0; v < n; ++v) if (model.nodes[v].has_birth) {
      for (int p : model.nodes[v].parents) {
        const int a = root(v), b = root(p);
        if (a != b) parent[a] = b;
      }
    }
    for (size_t v = 0; v < n; ++v) components[root(v)].push_back(v);
    }
    try {
      for (const auto& [r, nodes] : components) {
        if(cache.singleton_life&&nodes.size()==1&&!model.nodes[nodes[0]].has_birth){
          const auto& node=model.nodes[nodes[0]];const int lo=life_value(life.birth_lo,life.borrowed_birth_lo,node.typed_id),hi=life_value(life.birth_hi,life.borrowed_birth_hi,node.typed_id);
          // The NEW component closure proves no incident birth/parent role.
          // Its entire constraint system is just this one birth interval.
          // Contradictions retain the original full verifier and reason.
          if(lo==INT_MIN||hi==INT_MIN||lo<=hi){++cache.singleton_proofs;continue;}
        }
        auto start = Clock::now();
        std::vector<int> index(n, -1);
        for (size_t i = 0; i < nodes.size(); ++i) index[nodes[i]] = i;
        std::string key;
        auto add = [&](const auto& x) {
          key.append(reinterpret_cast<const char*>(&x), sizeof(x));
        };
        add(life.min_age); add(life.max_age); add(nodes.size());
        for (int v : nodes) {
          const auto& node = model.nodes[v];
          add(node.has_birth);
          for (int p : node.parents) { const int local = p < 0 ? -1 : index[p]; add(local); }
          for (int value : {life_sex_value(life, node.typed_id),
              life_value(life.birth_lo, life.borrowed_birth_lo, node.typed_id),
              life_value(life.birth_hi, life.borrowed_birth_hi, node.typed_id),
              life_value(life.year_last, life.borrowed_year_last, node.typed_id)}) add(value);
        }
        ++cache.requests;
        const bool hit = cache.valid.count(key);
        cache.lookup_seconds += std::chrono::duration<double>(Clock::now() - start).count();
        if (hit) { ++cache.hits; continue; }
        TopologyModel component;
        for (int v : nodes) {
          auto node = model.nodes[v];
          if (node.has_birth) for (int& p : node.parents) p = index[p];
          component.nodes.push_back(std::move(node));
        }
        start = Clock::now();
        verify_life(component, life);
        cache.verify_seconds += std::chrono::duration<double>(Clock::now() - start).count();
        ++cache.builds;
        const size_t charge = sizeof(std::string) + key.capacity() + 1 + 4 * sizeof(void*);
        if (charge <= cache.limit - cache.bytes) {
          cache.valid.insert(std::move(key));
          cache.bytes += charge;
        } else ++cache.refusals;
      }
    } catch (const Failure&) {
      // Recover the original whole-model classification and exact reason.
      verify_life(model, life);
    }
  }

  static bool verify_relative_life(const TopologyModel& model,const Life& life) {
    namespace rl=pap_internal::relative_life;
    auto* cache=rl::current;if(!cache)return false;if(!cache->limit){return false;}
    const size_t n=model.nodes.size();
    if(!n||n>1024){return false;}
    std::array<long long,2> minimum,maximum;std::array<bool,2> bounded;
    for(int side=0;side<2;++side){double lo=life.min_age[side],hi=life.max_age[side];if(!std::isfinite(lo)||std::fabs(lo)>INT_MAX||std::isnan(hi)||hi<lo||(std::isfinite(hi)&&std::fabs(hi)>INT_MAX)){return false;}
      minimum[side]=integer_age(lo,true);bounded[side]=std::isfinite(hi);maximum[side]=bounded[side]?integer_age(hi,false):0;}
    std::vector<int> years(n,INT_MIN),last(n),roles(n),parent(n);bool anchored=false;
    for(size_t v=0;v<n;++v){parent[v]=v;const auto& node=model.nodes[v];last[v]=life_value(life.year_last,life.borrowed_year_last,node.typed_id);if(node.typed_id>0){int lo=life_value(life.birth_lo,life.borrowed_birth_lo,node.typed_id),hi=life_value(life.birth_hi,life.borrowed_birth_hi,node.typed_id);if(lo==INT_MIN||lo!=hi){return false;}years[v]=lo;anchored=true;}}
    if(!anchored){return false;}

    auto root=[&](int v){while(parent[v]!=v){parent[v]=parent[parent[v]];v=parent[v];}return v;};bool typed_feasible=true;
    // Substitution of fixed typed years leaves independent latent components.
    // Every typed--typed constraint and every sex requirement is still checked.
    for(size_t v=0;v<n;++v)if(model.nodes[v].has_birth)for(int side=0;side<2;++side){int p=model.nodes[v].parents[side];roles[p]|=1<<side;
      if(years[v]!=INT_MIN&&years[p]!=INT_MIN){
        const int child_id=model.nodes[v].typed_id,parent_id=model.nodes[p].typed_id;
        if(child_id>0&&size_t(child_id)<=cache->source_parents.size()&&cache->source_parents[child_id-1][side]==parent_id){continue;}
        long long age=static_cast<long long>(years[v])-years[p];bool edge_ok=age>=minimum[side]&&(!bounded[side]||age<=maximum[side])&&(last[p]==INT_MIN||years[v]<=last[p]);typed_feasible&=edge_ok;}
      else if(years[v]==INT_MIN&&years[p]==INT_MIN){int a=root(v),b=root(p);if(a!=b)parent[a]=b;}
    }
    for(size_t v=0;v<n;++v){int sex=life_sex_value(life,model.nodes[v].typed_id);if(roles[v]==3){inapplicable("Ordinary parent-role conflict for "+model.nodes[v].label);}if((sex==1&&(roles[v]&2))||(sex==2&&(roles[v]&1))){inapplicable("Known sex incompatible with parent side for "+model.nodes[v].label);}}

    std::map<int,std::vector<int>> groups;for(size_t v=0;v<n;++v)if(years[v]==INT_MIN)groups[root(v)].push_back(v);
    for(size_t v=0;v<n;++v)if(model.nodes[v].has_birth)for(int p:model.nodes[v].parents){if(years[v]==INT_MIN&&years[p]!=INT_MIN)groups[root(v)].push_back(p);if(years[v]!=INT_MIN&&years[p]==INT_MIN)groups[root(p)].push_back(v);}
    if(groups.empty())groups[-1]={};
    auto add=[](std::string& s,const auto& value){s.append(reinterpret_cast<const char*>(&value),sizeof(value));};
    struct Missing {rl::Template* shape;std::string key;std::vector<int> nodes;};std::vector<Missing> missing;bool all_ready=true;
    for(auto& item:groups){auto& nodes=item.second;std::sort(nodes.begin(),nodes.end());nodes.erase(std::unique(nodes.begin(),nodes.end()),nodes.end());
      // A free founder with one fixed-year child has one nonempty interval.
      // It contributes no coupled feasibility question and no cache key.
      if(nodes.size()==2){int latent=years[nodes[0]]==INT_MIN?nodes[0]:nodes[1],child=latent==nodes[0]?nodes[1]:nodes[0];
        if(years[latent]==INT_MIN&&years[child]!=INT_MIN&&!model.nodes[latent].has_birth&&model.nodes[child].has_birth){int count=0,side_used=-1;for(int side=0;side<2;++side)if(model.nodes[child].parents[side]==latent){++count;side_used=side;}
          if(count==1&&(!bounded[side_used]||minimum[side_used]<=maximum[side_used])){continue;}}}
      std::vector<int> index(n,-1),local_roles(n);for(size_t i=0;i<nodes.size();++i)index[nodes[i]]=i;
      for(int v:nodes)if(model.nodes[v].has_birth)for(int side=0;side<2;++side){int p=model.nodes[v].parents[side];if(index[p]>=0&&(years[v]==INT_MIN||years[p]==INT_MIN))local_roles[p]|=1<<side;}
      std::string shape,key;add(shape,nodes.size());add(shape,life.min_age);add(shape,life.max_age);long long anchor=0;for(int v:nodes)if(years[v]!=INT_MIN){anchor=years[v];break;}
      for(int v:nodes){bool typed=years[v]!=INT_MIN;add(shape,typed);for(int p:model.nodes[v].parents){int mapped=p<0||typed&&years[p]!=INT_MIN?-1:index[p];add(shape,mapped);}
        long long birth=typed?static_cast<long long>(years[v])-anchor:0;add(key,birth);bool dated=local_roles[v]&&last[v]!=INT_MIN;add(key,dated);if(dated){long long relative=static_cast<long long>(last[v])-anchor;add(key,relative);}}

      auto found=cache->shapes.find(shape);if(found==cache->shapes.end()){size_t charge=shape.capacity()+nodes.size()*sizeof(int)+sizeof(rl::Template)+128;if(charge>cache->limit-cache->bytes){return false;}rl::Template t;for(int v:nodes)t.roles.push_back(roles[v]);found=cache->shapes.emplace(std::move(shape),std::move(t)).first;cache->bytes+=charge;}
      auto& t=found->second;
      if(t.outcomes.count(key)){}else{all_ready=false;missing.push_back({&t,std::move(key),nodes});}

    }

    if(typed_feasible&&all_ready){return true;}
    // Failed diagnostics retain the complete absolute numeric system and its
    // original constraint order; only display labels are rebound on lookup.

    std::string negative;add(negative,n);add(negative,life.min_age);add(negative,life.max_age);
    for(size_t v=0;v<n;++v){const auto& node=model.nodes[v];add(negative,node.has_birth);add(negative,node.parents);bool typed=node.typed_id>0;add(negative,typed);add(negative,years[v]);add(negative,last[v]);}
    auto format_failure=[&](const rl::Outcome& result){std::ostringstream text;text<<"Joint birth/age/Year.last constraints infeasible:";for(const auto& edge:result.cycle){const auto& node=model.nodes[edge.node];std::string reason;
      if(edge.kind==1)reason=node.label+" birth >= "+integer_string(life_value(life.birth_lo,life.borrowed_birth_lo,node.typed_id));else if(edge.kind==2)reason=node.label+" birth <= "+integer_string(life_value(life.birth_hi,life.borrowed_birth_hi,node.typed_id));else if(edge.kind==3||edge.kind==4)reason=model.nodes[edge.parent].label+" -> "+node.label+(edge.kind==3?" minimum ":" maximum ")+side_name(edge.side)+" age";else reason=model.nodes[edge.parent].label+" Year.last >= birth("+node.label+")";text<<" ["<<reason<<"]";}return text.str();};

    if(auto at=cache->negative.find(negative);at!=cache->negative.end()){inapplicable(format_failure(at->second));}

    // Fixed typed years separate the latent variables into independent blocks.
    // Solve only missing blocks, with typed--typed source constraints removed.
    // The original full solver is retained solely for authoritative diagnostics
    // when a block (or a newly proposed typed edge) is infeasible.

    bool feasible=typed_feasible;
    for(auto& m:missing){
      if(m.shape->outcomes.count(m.key))continue;
      TopologyModel component;std::vector<int> index(n,-1);
      for(int v:m.nodes){index[v]=component.nodes.size();component.nodes.push_back(model.nodes[v]);}
      std::vector<std::array<bool,2>> edges(component.nodes.size(),{{false,false}});
      for(size_t i=0;i<m.nodes.size();++i){int v=m.nodes[i];auto& node=component.nodes[i];
        if(!node.has_birth)continue;
        for(int side=0;side<2;++side){int p=model.nodes[v].parents[side];bool active=index[p]>=0&&(years[v]==INT_MIN||years[p]==INT_MIN);edges[i][side]=active;node.parents[side]=active?index[p]:-1;}}

      try{verify_life(component,life,nullptr,&edges);}catch(const Failure&){feasible=false;}

      if(!feasible)break;

      size_t charge=m.key.capacity()+sizeof(rl::Outcome)+128;
      if(charge<=cache->limit-cache->bytes){m.shape->outcomes.emplace(std::move(m.key),rl::Outcome{true,true,{}});cache->bytes+=charge;}

    }
    if(feasible)return true;

    rl::Outcome result;
    try{verify_life(model,life,&result);}catch(const Failure&){if(result.complete&&!result.feasible){size_t charge=negative.capacity()+result.cycle.capacity()*sizeof(rl::Constraint)+sizeof(rl::Outcome)+128;if(charge<=cache->limit-cache->bytes){cache->negative.emplace(std::move(negative),std::move(result));cache->bytes+=charge;}}throw;}

    return true;
  }

  static void verify_life(const TopologyModel& model, const Life& life, pap_internal::relative_life::Outcome* captured=nullptr, const std::vector<std::array<bool,2>>* enabled_edges=nullptr) {

    // Both successful return and exception immediately decide applicability.

    const std::size_t size = model.nodes.size();
    const int zero = static_cast<int>(size);
    std::vector<int> roles(size, 0);
    std::vector<DifferenceConstraint> constraints;
    std::vector<pap_internal::relative_life::Constraint> descriptors;if(captured)descriptors.reserve(8*size);
    for (std::size_t child = 0; child < size; ++child) {
      const TopologyNode& node = model.nodes[child];
      const int lo = life_value(life.birth_lo, life.borrowed_birth_lo, node.typed_id);
      const int hi = life_value(life.birth_hi, life.borrowed_birth_hi, node.typed_id);
      if (lo != INT_MIN)
        constraints.push_back(DifferenceConstraint(static_cast<int>(child), zero,
          -static_cast<long long>(lo), node.label + " birth >= " + integer_string(lo)));
      if(captured&&lo!=INT_MIN)descriptors.push_back({1,static_cast<int>(child),-1,-1});
      if (hi != INT_MIN)
        constraints.push_back(DifferenceConstraint(zero, static_cast<int>(child),
          static_cast<long long>(hi), node.label + " birth <= " + integer_string(hi)));
      if(captured&&hi!=INT_MIN)descriptors.push_back({2,static_cast<int>(child),-1,-1});
      if (!node.has_birth) continue;
      for (int side = 0; side < 2; ++side) {
        if(enabled_edges&&!(*enabled_edges)[child][side])continue;
        const int parent = node.parents[static_cast<std::size_t>(side)];
        roles[static_cast<std::size_t>(parent)] |= 1 << side;
        const TopologyNode& source = model.nodes[static_cast<std::size_t>(parent)];
        const double lower = life.min_age[static_cast<std::size_t>(side)];
        const double upper = life.max_age[static_cast<std::size_t>(side)];
        if (std::isnan(lower) || std::isnan(upper) || upper < lower)
          throw Failure(TopologyStatus::Unsupported, "Invalid side-specific age bounds");
        const long long minimum = integer_age(lower, true);
        constraints.push_back(DifferenceConstraint(static_cast<int>(child), parent, -minimum,
          source.label + " -> " + node.label + " minimum " + side_name(side) + " age"));
        if(captured)descriptors.push_back({3,static_cast<int>(child),parent,side});
        if (std::isfinite(upper)) {
          const long long maximum = integer_age(upper, false);
          constraints.push_back(DifferenceConstraint(parent, static_cast<int>(child), maximum,
            source.label + " -> " + node.label + " maximum " + side_name(side) + " age"));
          if(captured)descriptors.push_back({4,static_cast<int>(child),parent,side});
        } else if (upper < 0.0) {
          throw Failure(TopologyStatus::Unsupported, "Invalid infinite maximum age");
        }
        const int last = life_value(life.year_last, life.borrowed_year_last, source.typed_id);
        if (last != INT_MIN)
          constraints.push_back(DifferenceConstraint(zero, static_cast<int>(child),
            static_cast<long long>(last), source.label + " Year.last >= birth(" + node.label + ")"));
        if(captured&&last!=INT_MIN)descriptors.push_back({5,static_cast<int>(child),parent,side});
      }
    }
    for (std::size_t node = 0; node < size; ++node) {
      const int sex = life_sex_value(life, model.nodes[node].typed_id);
      if (roles[node] == 3)
        inapplicable("Ordinary parent-role conflict for " + model.nodes[node].label);
      if ((sex == 1 && (roles[node] & 2) != 0) ||
          (sex == 2 && (roles[node] & 1) != 0))
        inapplicable("Known sex incompatible with parent side for " + model.nodes[node].label);
    }
    // All nodes start at zero, equivalent to adding a super-source. This
    // detects ANY infeasible cycle, including constraints on latent bridges.
    const std::size_t count = size + 1U;
    std::vector<long long> distances(count, 0);
    std::vector<int> previous(count, -1);
    int changed = -1;
    for (std::size_t step = 0; step < count; ++step) {
      changed = -1;
      for (std::size_t edge = 0; edge < constraints.size(); ++edge) {
        const DifferenceConstraint& item = constraints[edge];
        const long long origin = distances[static_cast<std::size_t>(item.from)];
        if ((item.upper < 0 && origin < std::numeric_limits<long long>::min() - item.upper) ||
            (item.upper > 0 && origin > std::numeric_limits<long long>::max() - item.upper))
          throw Failure(TopologyStatus::Unsupported, "Life-history integer work exceeds safe range");
        const long long candidate = origin + item.upper;
        if (candidate < distances[static_cast<std::size_t>(item.to)]) {
          distances[static_cast<std::size_t>(item.to)] = candidate;
          previous[static_cast<std::size_t>(item.to)] = static_cast<int>(edge);
          changed = item.to;
        }
      }
      if (changed < 0) {
        if(captured){captured->complete=true;captured->feasible=true;}
        remember_life_witness(model, distances, zero);
        return;
      }
    }
    int node = changed;
    for (std::size_t step = 0; step < count; ++step)
      node = constraints[static_cast<std::size_t>(previous[static_cast<std::size_t>(node)])].from;
    const int start = node;
    std::ostringstream reason;
    reason << "Joint birth/age/Year.last constraints infeasible:";
    if(captured){captured->complete=true;captured->feasible=false;}
    do {
      const DifferenceConstraint& item = constraints[
        static_cast<std::size_t>(previous[static_cast<std::size_t>(node)])];
      if(captured)captured->cycle.push_back(descriptors[static_cast<std::size_t>(previous[static_cast<std::size_t>(node)])]);
      reason << " [" << item.reason << "]";
      node = item.from;
    } while (node != start);
    inapplicable(reason.str());
  }
};

// Restore literal source founders after a generic configuration builder has
// completed every parent slot. An unchanged source founder remains a founder
// in the raw model; its otherwise unused independent latent parents integrate
// to its HWE prior. Then trim construction-only nodes and refresh the order.
inline void restore_source_founders(TopologyModel& model, const Input& source,
                                    bool focal_has_proposed_parent) {
  for (std::size_t original = 0U; original < source.nodes.size(); ++original) {
    const auto& node = source.nodes[original];
    if (node.has_birth || node.parents[0] >= 0 || node.parents[1] >= 0 ||
        (static_cast<int>(original) == source.focal && focal_has_proposed_parent))
      continue;
    const int canonical = model.original_to_canonical[original];
    if (canonical < 0) continue;
    auto& restored = model.nodes[static_cast<std::size_t>(canonical)];
    restored.has_birth = false;
    restored.parents = {{-1, -1}};
  }

  std::vector<unsigned char> keep(model.nodes.size(), 0U);
  std::vector<int> pending;
  for (std::size_t i = 0U; i < model.nodes.size(); ++i)
    if (model.nodes[i].typed_id > 0 || model.nodes[i].observed)
      pending.push_back(static_cast<int>(i));
  while (!pending.empty()) {
    const int at = pending.back(); pending.pop_back();
    if (keep[static_cast<std::size_t>(at)]) continue;
    keep[static_cast<std::size_t>(at)] = 1U;
    const auto& node = model.nodes[static_cast<std::size_t>(at)];
    if (node.has_birth)
      for (int parent : node.parents) pending.push_back(parent);
  }

  std::vector<int> remap(model.nodes.size(), -1);
  std::vector<TopologyNode> nodes;
  for (std::size_t i = 0U; i < model.nodes.size(); ++i)
    if (keep[i]) {
      remap[i] = static_cast<int>(nodes.size());
      nodes.push_back(std::move(model.nodes[i]));
    }
  for (auto& node : nodes)
    if (node.has_birth)
      for (int& parent : node.parents)
        parent = remap[static_cast<std::size_t>(parent)];
  model.focal = model.focal >= 0 ? remap[static_cast<std::size_t>(model.focal)] : -1;
  model.candidate = model.candidate >= 0 ?
    remap[static_cast<std::size_t>(model.candidate)] : -1;
  for (int& canonical : model.original_to_canonical)
    if (canonical >= 0) canonical = remap[static_cast<std::size_t>(canonical)];
  model.nodes = std::move(nodes);

  std::vector<int> degree(model.nodes.size(), 0);
  std::vector<std::vector<int>> children(model.nodes.size());
  std::set<int> ready;
  for (std::size_t i = 0U; i < model.nodes.size(); ++i) {
    const auto& node = model.nodes[i];
    if (node.has_birth)
      for (int parent : node.parents) {
        if (parent < 0 || static_cast<std::size_t>(parent) >= model.nodes.size())
          throw Failure(TopologyStatus::Unsupported,
            "Restored source founder has an incomplete parent identity");
        ++degree[i];
        children[static_cast<std::size_t>(parent)].push_back(static_cast<int>(i));
      }
    if (degree[i] == 0) ready.insert(static_cast<int>(i));
  }
  model.topological_order.clear();
  while (!ready.empty()) {
    const int parent = *ready.begin();
    ready.erase(ready.begin());
    model.topological_order.push_back(parent);
    for (int child : children[static_cast<std::size_t>(parent)])
      if (--degree[static_cast<std::size_t>(child)] == 0) ready.insert(child);
  }
  if (model.topological_order.size() != model.nodes.size())
    throw Failure(TopologyStatus::Inapplicable,
                  "Restored source topology contains a directed cycle");
  model.provenance.push_back("Unchanged source founder retains its literal HWE factor");
}

inline std::vector<TopologyModel> specifications() {
  std::vector<TopologyModel> result;
  const char* simple[] = {"PO", "U", "FS"};
  for (int i = 0; i < 3; ++i) {
    TopologyModel model;
    model.name = simple[i]; model.relationship = simple[i];
    result.push_back(model);
  }
  for (int side = 1; side <= 2; ++side) {
    TopologyModel model;
    model.relationship = "HS";
    model.name = std::string("HS_") + side_name(side - 1);
    model.shared_side = side; model.orientation_index = side - 1;
    result.push_back(model);
  }
  int orientation = 2;
  const char* families[] = {"GP", "FA", "HA"};
  for (int family = 0; family < 3; ++family)
    for (int direction = 0; direction < 2; ++direction)
      for (int bridge = 1; bridge <= 2; ++bridge)
        for (int shared = 1; shared <= (family == 1 ? 1 : 2); ++shared) {
          TopologyModel model;
          model.relationship = families[family];
          model.focal_is_descendant = direction == 0;
          model.bridge_side = bridge;
          model.shared_side = family == 1 ? 0 : shared;
          model.orientation_index = orientation++;
          model.name = std::string(families[family]) + "_" +
            (direction == 0 ? "focal" : "candidate") + "_" + side_name(bridge - 1);
          if (family != 1) model.name += std::string("_") + side_name(shared - 1);
          result.push_back(model);
        }
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i].model_order = static_cast<int>(i);
  return result;
}

}  // namespace topology_detail

// Ordered: PO,U,FS,HS2,GP8,FA4,HA8. Status is independent of genetic likelihood.
inline std::vector<TopologyModel> build_models(const Input& input) {
  std::vector<TopologyModel> models = topology_detail::specifications();
  // One immutable, validated source for this exact current Input. It is never
  // retained across a graph mutation or a different proposal boundary.
  auto* cache=pap_internal::topology_life_reuse::active;
  std::unique_ptr<topology_detail::Builder> source;
  if(cache&&cache->reuse_builder&&cache->limit)try{source=std::make_unique<topology_detail::Builder>(input);++cache->source_builds;}
    catch(const topology_detail::Failure&){}catch(const std::bad_alloc&){}catch(const std::length_error&){}
  for (std::size_t i = 0; i < models.size(); ++i) {
    try {
      if(topology_detail::Builder::try_bound(input,models[i]))continue;
      topology_detail::Builder builder=source?*source:topology_detail::Builder(input);
      if(source)++cache->source_copies;
      builder.apply(models[i]);
      builder.finish(models[i]);
      topology_detail::Builder::remember_bound(input,models[i]);
    } catch (const topology_detail::Failure& error) {
      models[i].status = error.status;
      models[i].reason = error.what();
      models[i].nodes.clear();
      models[i].topological_order.clear();
      models[i].original_to_canonical.clear();
      models[i].focal = models[i].candidate = -1;
    }
  }
  return models;
}

}  // namespace relationship_topology

#endif
