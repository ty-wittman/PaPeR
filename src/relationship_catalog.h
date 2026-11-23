#ifndef PAP_RELATIONSHIP_CATALOG_H
#define PAP_RELATIONSHIP_CATALOG_H

// Literal reviewed Stage-1 full DAGs on the same immutable bounded ordinary
// Input. The catalog owns relationship motifs/biological slot variants only;
// the existing factor compiler owns their unaltered raw log partitions.
// Missing ancestor roles are latent, retained typed identities are fixed and
// all observations remain model-owned exactly once. No ancestry is published.

#include "relationship_topology.h"
#include "input_checks.h"
#include "motif_templates.h"
#include <utility>

namespace relationship_topology {
namespace full_catalog {

struct Birth {
  std::string child;
  std::array<std::string, 2> parents;
};
struct Definition {
  std::string relationship;
  bool parent_numerator = false;
  bool registered_class = true;
  bool reverse_allowed = false;
  std::vector<Birth> births;
  // A marginally symmetric registry row can still have a directed motif.
  // This flag exposes an additional biologically required endpoint swap
  // without pretending that its production orientation mask already owns it.
  bool reverse_review_required = false;
};
struct Member {
  std::string relationship;
  bool parent_numerator = false;
  bool registered_class = true;
  bool reverse = false;
  bool ownership_review_required = false;
  bool registered_endpoint_orientation = true;
  bool orientation_review_required = false;
  unsigned int parent_slot_bits = 0U;
  TopologyModel topology;
};

inline std::vector<Definition> definitions(bool include_unregistered = false,
                                           bool include_s37_parent_forms = false) {
  // A/B are the canonical focal endpoints. Parent role permutations below
  // are essential: canonical genetic motifs are not sex-oriented pedigrees.
  std::vector<Definition> out{
    {"HS_PO", true, true, true,
      {{"B", {{"P", "Q"}}}, {"A", {{"B", "P"}}}}},
    {"PO_HA", true, true, false,
      {{"B", {{"G", "H"}}}, {"C", {{"G", "J"}}}, {"A", {{"B", "C"}}}}},
    {"HS_HA", false, true, true,
      {{"B", {{"P", "R"}}}, {"Q", {{"R", "S"}}}, {"A", {{"P", "Q"}}}}, true},
    {"HS_F1C", false, true, false,
      {{"Q", {{"G", "H"}}}, {"R", {{"G", "H"}}},
       {"A", {{"P", "Q"}}}, {"B", {{"P", "R"}}}}},
    {"HS_FA", false, true, true,
      {{"B", {{"P", "R"}}}, {"Q", {{"P", "R"}}}, {"A", {{"P", "Q"}}}}},
    {"HS_GP", false, true, true,
      {{"B", {{"P", "S"}}}, {"Q", {{"B", "R"}}}, {"A", {{"P", "Q"}}}}},
    {"HS_GGG", false, true, true,
      {{"B", {{"P", "S"}}}, {"C", {{"B", "R"}}},
       {"Q", {{"C", "T"}}}, {"A", {{"P", "Q"}}}}},
    {"FA_HA", false, true, true,
      {{"B", {{"G", "H"}}}, {"P", {{"G", "H"}}},
       {"Q", {{"G", "J"}}}, {"A", {{"P", "Q"}}}}},
    {"FS_HA", false, true, false,
      {{"Q", {{"P", "R"}}}, {"A", {{"P", "Q"}}}, {"B", {{"P", "Q"}}}}},
    {"FS_H1C", false, true, false,
      {{"P", {{"G", "R"}}}, {"Q", {{"G", "S"}}},
       {"A", {{"P", "Q"}}}, {"B", {{"P", "Q"}}}}},
    {"FS_PARENTS_FS", false, true, false,
      {{"P", {{"G", "H"}}}, {"Q", {{"G", "H"}}},
       {"A", {{"P", "Q"}}}, {"B", {{"P", "Q"}}}}}
  };
  // Reviewed S37 parent form, appended to keep all older integer identities.
  // Existing default catalog consumers retain their original question set.
  // The public PO_FA alias of HS_PO is NOT this graph. This form supplies
  // an actual tested parent edge and is never a reversed/rival hypothesis.
  if (include_s37_parent_forms || include_unregistered)
    out.push_back({"PO_FA_co_parent_FS", true, true, false,
      {{"B", {{"G", "H"}}}, {"C", {{"G", "H"}}}, {"A", {{"B", "C"}}}}});
  if (include_unregistered) {
    out.push_back({"PO_GP_through_co_parent", true, false, false,
      {{"C", {{"B", "Q"}}}, {"A", {{"B", "C"}}}}});
  }
  return out;
}

namespace detail {

inline void fail(TopologyStatus status, const std::string& reason) {
  throw topology_detail::Failure(status, reason);
}

inline std::vector<std::string> roles(const Definition& definition) {
  std::set<std::string> names{"A", "B"};
  for (const Birth& birth : definition.births) {
    names.insert(birth.child);
    names.insert(birth.parents[0]);
    names.insert(birth.parents[1]);
  }
  return std::vector<std::string>(names.begin(), names.end());
}

inline const pap_internal::motif_templates::Program* compiled_motif(const Definition& d,unsigned bits){
 namespace mt=pap_internal::motif_templates;auto* cache=mt::active;if(!cache)return nullptr;++cache->requests;
 if(!cache->limit){++cache->refusals;return nullptr;}
 try{
  std::string key;auto integer=[&](std::size_t v){key.append(reinterpret_cast<const char*>(&v),sizeof(v));};auto text=[&](const std::string& v){integer(v.size());key+=v;};
  text(d.relationship);integer(d.parent_numerator);integer(d.registered_class);integer(d.reverse_allowed);integer(d.reverse_review_required);integer(bits);integer(d.births.size());
  for(const auto& b:d.births){text(b.child);text(b.parents[0]);text(b.parents[1]);}
  auto at=cache->entries.find(key);if(at!=cache->entries.end()){++cache->hits;return at->second.get();}
  auto p=std::make_shared<mt::Program>();const auto names=roles(d);p->roles=names.size();
  auto id=[&](const std::string& name){return int(std::lower_bound(names.begin(),names.end(),name)-names.begin());};p->a=id("A");p->b=id("B");
  std::map<std::string,int> parental;for(size_t row=0;row<d.births.size();++row){int flip=(bits>>row)&1U;parental[d.births[row].parents[0]]|=1<<flip;parental[d.births[row].parents[1]]|=1<<(1-flip);}
  for(const auto& r:parental)if(r.second==3){p->pre_failure="Full motif requires both parental roles for canonical "+r.first;break;}
  std::vector<unsigned char> bound(names.size()),done(d.births.size());bound[p->a]=bound[p->b]=1;size_t remaining=done.size();
  while(remaining&&!p->pre_failure.size()){
   bool progress=false;for(size_t row=0;row<d.births.size();++row){const auto& b=d.births[row];if(done[row]||!bound[id(b.child)])continue;
    int flip=(bits>>row)&1U;for(int role=0;role<2;++role){int parent=id(b.parents[role]);p->operations.push_back({id(b.child),parent,1+(role^flip)});bound[parent]=1;}
    done[row]=1;--remaining;progress=true;
   }if(!progress){p->post_failure="Full motif contains an unanchored nonfounder role";break;}
  }
  // An unanchored failure is reported AFTER preceding dynamic unifications,
  // preserving their original typed-conflict/error precedence.
  if(p->pre_failure.empty()&&p->post_failure.empty())for(size_t a=0;a<names.size();++a)for(size_t b=a+1;b<names.size();++b)
   p->distinct.push_back({int(a),int(b),"Full motif distinct canonical roles coincide: "+names[a]+","+names[b]});
  size_t bytes=key.capacity()+sizeof(mt::Program)+sizeof(void*)*8+p->operations.capacity()*sizeof(mt::Operation)+p->distinct.capacity()*sizeof(mt::Distinct)+p->pre_failure.capacity()+p->post_failure.capacity();
  for(const auto& c:p->distinct)bytes+=c.reason.capacity();if(bytes>cache->limit-cache->bytes){++cache->refusals;return nullptr;}
  auto result=cache->entries.emplace(std::move(key),p);cache->bytes+=bytes;++cache->builds;return result.first->second.get();
 }catch(const std::bad_alloc&){++cache->refusals;return nullptr;}catch(const std::length_error&){++cache->refusals;return nullptr;}
}
inline void apply_motif(const pap_internal::motif_templates::Program& p,topology_detail::Builder& builder,const Input& input,bool reverse){
 std::vector<int> bound(p.roles,-1);bound[p.a]=reverse?input.candidate:input.focal;bound[p.b]=reverse?input.focal:input.candidate;
 for(const auto& op:p.operations){int slot=builder.parent_slot(bound[op.child],op.side);if(bound[op.parent]<0)bound[op.parent]=slot;else bound[op.parent]=builder.unify_nodes(slot,bound[op.parent]);}
 if(!p.post_failure.empty())fail(TopologyStatus::Unsupported,p.post_failure);
 for(const auto& c:p.distinct)builder.require_distinct_nodes(bound[c.first],bound[c.second],c.reason);
}

// Check source slots before the motif materializes its own latent parents.
// Joint configuration callers also supply the proposed opposite parent.
inline void require_hspo_parentless(const Input& input, bool opposite_parent = false) {
  if (opposite_parent)
    fail(TopologyStatus::Inapplicable, "HS_PO requires neither endpoint to have assigned parents");
  for (int endpoint : {input.focal, input.candidate}) {
    if (endpoint < 0 || static_cast<std::size_t>(endpoint) >= input.nodes.size())
      fail(TopologyStatus::Unsupported, "HS_PO endpoint outside source boundary");
    const auto& node = input.nodes[static_cast<std::size_t>(endpoint)];
    if (node.parents[0] >= 0 || node.parents[1] >= 0)
      fail(TopologyStatus::Inapplicable, "HS_PO requires neither endpoint to have assigned parents");
  }
}

inline Member build(const Input& input, const Definition& definition,
                    bool reverse, unsigned int parent_slot_bits,
                    int model_order, const topology_detail::Builder* prepared_source=nullptr,
                    bool parentless_hspo=false) {
  Member member;
  member.relationship = definition.relationship;
  // A registered reverse HS_PO has candidate as offspring of focal. It is
  // a valid nonparent alternative to the actual tested candidate-parent PO,
  // not a stronger parent numerator or an impossible relationship class.
  member.parent_numerator = definition.parent_numerator && !reverse;
  member.registered_class = definition.registered_class;
  member.reverse = reverse;
  member.ownership_review_required = definition.parent_numerator && reverse;
  member.orientation_review_required = definition.reverse_review_required && reverse;
  member.registered_endpoint_orientation = definition.registered_class &&
                                          !member.orientation_review_required;
  member.parent_slot_bits = parent_slot_bits;
  TopologyModel& model = member.topology;
  model.relationship = definition.relationship;
  model.name = definition.relationship + (reverse ? "_reverse_slots_" : "_forward_slots_") +
               std::to_string(parent_slot_bits);
  model.model_order = model_order;
  model.focal_is_descendant = !reverse;
  // Admission precedes both structural-cache lookup and motif application.
  if (parentless_hspo && definition.relationship == "HS_PO" && !reverse) {
    try { require_hspo_parentless(input); }
    catch (const topology_detail::Failure& error) {
      model.status = error.status;
      model.reason = error.what();
      return member;
    }
  }
  std::string bound_tag;bool cache_tag=true;try{if(TopologyTemplateCache::active||FinishedModelTransfer::active){using pap_internal::topology_templates::append;using pap_internal::topology_templates::append_string;append(bound_tag,reverse);append(bound_tag,parent_slot_bits);append(bound_tag,definition.parent_numerator);append(bound_tag,definition.registered_class);append(bound_tag,definition.reverse_allowed);append(bound_tag,definition.reverse_review_required);for(const auto& birth:definition.births){append_string(bound_tag,birth.child);for(const auto& parent:birth.parents)append_string(bound_tag,parent);}}}catch(const std::bad_alloc&){cache_tag=false;if(TopologyTemplateCache::active)++TopologyTemplateCache::active->refusals;}catch(const std::length_error&){cache_tag=false;if(TopologyTemplateCache::active)++TopologyTemplateCache::active->refusals;}
  if(cache_tag&&(!pap_internal::input_checks::enabled()||!prepared_source||prepared_source->bound_to(input))&&topology_detail::Builder::try_bound(input,model,bound_tag))return member;
  try {
    if (input.proposed_parent_side != 1 && input.proposed_parent_side != 2)
      fail(TopologyStatus::Unsupported, "Full common boundary requires a tested dam or sire side");
    if (input.focal < 0 || static_cast<std::size_t>(input.focal) >= input.nodes.size())
      fail(TopologyStatus::Unsupported, "Full focal index outside common boundary");
    if (input.nodes[static_cast<std::size_t>(input.focal)].parents[
          static_cast<std::size_t>(input.proposed_parent_side - 1)] >= 0)
      fail(TopologyStatus::Unsupported, "Full common boundary requires a masked tested parent slot");
    if (definition.births.size() >= sizeof(unsigned int) * CHAR_BIT ||
        parent_slot_bits >= (1U << definition.births.size()))
      fail(TopologyStatus::Unsupported, "Full motif parent-slot variant outside finite catalog");
    if (reverse && !definition.reverse_allowed)
      fail(TopologyStatus::Inapplicable, "Full motif reverse endpoint orientation is not registered");
    // Stage 1 has one sex-compatible parental role per individual. A genetic
    // motif can require several canonical parent-slot permutations before
    // this condition is satisfied. Keep invalid variants explicitly tagged.
    const auto* motif=compiled_motif(definition,parent_slot_bits);
    if(motif&&!motif->pre_failure.empty())fail(TopologyStatus::Inapplicable,motif->pre_failure);
    if(!motif){
    std::map<std::string, int> parent_roles;
    for (std::size_t row = 0U; row < definition.births.size(); ++row) {
      const int flip = (parent_slot_bits >> row) & 1U;
      parent_roles[definition.births[row].parents[0]] |= 1 << flip;
      parent_roles[definition.births[row].parents[1]] |= 1 << (1 - flip);
    }
    for (const auto& role : parent_roles) if (role.second == 3)
      fail(TopologyStatus::Inapplicable,
           "Full motif requires both parental roles for canonical " + role.first);
    }
    if(pap_internal::input_checks::enabled()){if(prepared_source&&!prepared_source->bound_to(input))throw std::invalid_argument("full motif source belongs to another input");}
    auto builder=[&]{return prepared_source?*prepared_source:topology_detail::Builder(input);}();
    if(motif)apply_motif(*motif,builder,input,reverse);else {
    std::map<std::string, int> bound;
    bound["A"] = reverse ? input.candidate : input.focal;
    bound["B"] = reverse ? input.focal : input.candidate;
    std::vector<unsigned char> done(definition.births.size(), 0);
    std::size_t remaining = definition.births.size();
    while (remaining > 0U) {
      bool progress = false;
      for (std::size_t row = 0U; row < definition.births.size(); ++row) {
        if (done[row]) continue;
        const Birth& birth = definition.births[row];
        const auto child = bound.find(birth.child);
        if (child == bound.end()) continue;
        const int flip = (parent_slot_bits >> row) & 1U;
        for (int role = 0; role < 2; ++role) {
          const int slot = builder.parent_slot(child->second, 1 + (role ^ flip));
          const auto parent = bound.find(birth.parents[static_cast<std::size_t>(role)]);
          if (parent == bound.end()) bound[birth.parents[static_cast<std::size_t>(role)]] = slot;
          else parent->second = builder.unify_nodes(slot, parent->second);
        }
        done[row] = 1;
        --remaining;
        progress = true;
      }
      if (!progress)
        fail(TopologyStatus::Unsupported, "Full motif contains an unanchored nonfounder role");
    }
    const std::vector<std::string> motif_roles = roles(definition);
    for (const auto& role : motif_roles) if (!bound.count(role))
      fail(TopologyStatus::Unsupported, "Full motif role was not materialized");
    // The literal registry has distinct individual roles, except for its
    // explicitly repeated parent symbols. A retained typed contradiction or
    // implicit additional collapse defines a different hypothesis.
    for (std::size_t a = 0U; a < motif_roles.size(); ++a)
      for (std::size_t b = a + 1U; b < motif_roles.size(); ++b)
        builder.require_distinct_nodes(bound.at(motif_roles[a]), bound.at(motif_roles[b]),
          "Full motif distinct canonical roles coincide: " + motif_roles[a] + "," + motif_roles[b]);
    }
    if (member.parent_numerator &&
        builder.parent_slot(input.focal, input.proposed_parent_side) !=
          builder.unify_nodes(input.candidate, input.candidate))
      fail(TopologyStatus::Inapplicable,
           "Full parent numerator does not retain the actual tested direct PO orientation");
    builder.finish(model);
    model.provenance.push_back("Literal full catalog DAG with original masked common Input");
    model.provenance.push_back("All canonical role identities distinct except explicitly shared symbols");
    model.provenance.push_back("All observations/retained typed births remain boundary-owned; no ancestry publication");
    if (!definition.registered_class)
      model.provenance.push_back("Additional unregistered mathematical DAG; no production registry ownership inferred");
    if (member.ownership_review_required)
      model.provenance.push_back("Reverse HS_PO is a nonparent hypothesis for actual tested PO; historical row ownership requires scientific review before production promotion");
    if (member.orientation_review_required)
      model.provenance.push_back("HS_HA reversed directed DAG is a reviewed mathematical alternative; marginal registry symmetry does not establish typed-context endpoint closure");
    if(cache_tag)topology_detail::Builder::remember_bound(input,model,bound_tag);
  } catch (const topology_detail::Failure& error) {
    model.status = error.status;
    model.reason = error.what();
    model.nodes.clear();
    model.original_to_canonical.clear();
    model.topological_order.clear();
    model.focal = model.candidate = -1;
  }
  return member;
}

}  // namespace detail

// Return every registered direction, explicitly marked scientifically
// reviewed reverse HS_HA direction, and EVERY parental slot permutation,
// including biological Inapplicable and resource Unsupported variants. No
// numerical availability or marginal-coefficient identity is inferred here.
// Caller must score a fixed whole-panel topology before taking any maximum.
inline std::size_t member_count(const std::vector<Definition>& catalog) {
  std::size_t count = 0U;
  for (const auto& definition : catalog)
    count += (1U << definition.births.size()) *
      (definition.reverse_allowed ? 2U : 1U);
  return count;
}

template<class Visitor>
inline void visit_models(const Input& input,
                         const std::vector<Definition>& catalog,
                         Visitor&& visitor) {
  int model_order = 28;
  for (const Definition& definition : catalog) {
    const unsigned int variants = 1U << definition.births.size();
    for (int direction = 0; direction <= (definition.reverse_allowed ? 1 : 0); ++direction)
      for (unsigned int bits = 0U; bits < variants; ++bits)
        visitor(detail::build(input, definition, direction != 0,
                              bits, model_order++));
  }
}

template<class Visitor>
inline void visit_models(const Input& input, Visitor&& visitor,
                         bool include_unregistered = false) {
  const auto catalog = definitions(include_unregistered);
  visit_models(input, catalog, std::forward<Visitor>(visitor));
}

// Retain the complete-vector interface for diagnostic/reference consumers.
// Production preparation can consume and release one literal graph at a time.
inline std::vector<Member> build_models(const Input& input,
                                       bool include_unregistered = false) {
  const auto catalog = definitions(include_unregistered);
  std::vector<Member> result;
  result.reserve(member_count(catalog));
  visit_models(input, catalog, [&](Member&& member) {
    result.push_back(std::move(member));
  });
  return result;
}

}  // namespace full_catalog
}  // namespace relationship_topology

#endif
