#pragma once
#include "forest_messages.h"
namespace pap_internal::pedigree_likelihood {
std::vector<int> choose_order(const std::vector<std::vector<int>>& scopes, const std::vector<unsigned char>& pairing,
                                     const std::vector<int>& dom, int retained, int* width);

// ---------------------------------------------------------------------------
// Per-locus transmission tables shared by every question and visit.
struct Transmission {
  bool implicit = false;
  bool two_support = false;
  std::vector<int> sparse_allele;
  std::vector<double> sparse_weight;
  std::vector<std::array<int,2>> alleles;
  std::vector<double> pairing;  // [c = genotype of alleles (a, b)], index c + G a + G k b
  std::vector<double> one;      // T(c | x) with a population mate, index c + G x
  std::vector<int> partner;     // the allele a with genotype(a, b) = c, index c + G b; -1 none
  double one_value(const LocusView& L, std::size_t index) const { pap_internal::removed_path("one_value"); }
};

void build_transmission(const LocusView& L, Transmission& T, bool implicit = false);

// Parent-pair weights of a family's non-local children, W(d, s) =
// prod_c sum_g T(g | d, s) m_{c->f}(g). One slot per (family, children) key,
// computed lazily per locus and shared by every question of the visit.
struct PairCache {
  struct Slot { int family = -1; bool dam = false, sire = false; std::vector<int> edges; std::vector<Scaled> values; std::vector<unsigned char> ready; };
  std::map<std::vector<int>, int> index;
  std::vector<Slot> slots;
  std::size_t computed = 0, requests = 0;
  void clear() { index.clear(); slots.clear(); computed = requests = 0; }
  int slot(int family, bool dam, bool sire, const std::vector<int>& edges);
  const Scaled& get(int id, int locus, int loci, const Forest& F, const LocusView& L, const LocusMessages& M, Context& ctx) {
    ++requests;
    Slot& s = slots[id];
    if (s.ready.empty()) { s.values.resize(loci); s.ready.assign(loci, 0); }
    Scaled& out = s.values[locus];
    if (s.ready[locus]) return out;
    const bool dam = s.dam, sire = s.sire;
    const std::size_t cells = static_cast<std::size_t>(dam ? L.G : 1) * (sire ? L.G : 1);
    out.log_cells.clear();out.cells.assign(cells, 1.0); out.scale = 0.0;
    if(ctx.latent_brood_exclusion && (F.families[s.family].dam==0 || F.families[s.family].sire==0)) {
      s.ready[locus]=1;++computed;return out;
    }
    Scaled table;
    for (int e : s.edges) {
      const Scaled& c = ctx.to_family(e, locus, L, M);
      child_table_scaled(L,c,dam,sire,table);table.scale+=c.scale;
      multiply_into(out,table);normalize(out);
      if(out.scale==-INFINITY)break;
    }
    s.ready[locus] = 1; ++computed;
    return out;
  }
};

// ---------------------------------------------------------------------------
// Locus-independent plan for one hypothesis model.
struct Program;
namespace resource_plan {struct Execution;}

struct Plan {
  ResourceLimits limits;
  bool factorized_one = false;
  bool ancestral = false;
  bool local_family_pairs = false;
  bool latent_brood_exclusion = false;
  int restored_family_pairs = 0;
  bool high_allele_programs = false;
  enum Source : unsigned char { Emission, Prior, Incoming, Outgoing };
  struct Unary { int var; Source source; int index; };   // typed id (Emission) or edge
  struct Birth { std::array<int, 3> vars; std::array<int, 2> gametes{{-1, -1}}; bool excluded_brood=false; };  // child, dam, sire
  struct Pair { int slot; std::array<int, 2> vars; };     // dam, sire (-1 absent)
  struct PieceRef { int slot; std::vector<int> vars; };   // cached forest piece over these variables
  int barren_removed=0,barren_prior_normalizers=0,barren_population_normalizers=0;
  int off_scope_observations_before=0,descendant_observations_before=0,off_scope_messages=0;
  int retained = 0;                      // piece plans: the first `retained` variables are kept
  std::vector<PieceRef> pieces;
  int variables = 0, genotype_variables = 0, width = 0, family_terms = 0, steiner_individuals = 0, redundant_births = 0;
  std::vector<Unary> unaries;
  std::vector<Birth> births;
  std::vector<Pair> pairs;
  std::vector<unsigned char> allele;   // per variable: gamete variable with k states
  std::vector<int> touched_components;   // forest components holding local individuals
  bool family_side_untouched = false;    // the focal's former birth family side, when split off
  mutable std::vector<std::shared_ptr<const Program>> programs;  // by allele count k
  mutable std::vector<int> order;   // elimination order, chosen once per plan
  mutable int order_width = 0;
  bool ready = false;
};

// Directed inheritance cycles are biologically inapplicable, not a
// numerical refusal or a reason to change the key's interpretation.
bool cyclic_plan(const Plan& p);
void require_acyclic_plan(const Plan& p);

// Forest pieces: connected forest-only parts of a plan, eliminated onto the
// model individuals they touch. Hypotheses of a visit share them by key; each
// is computed once per locus.
struct PieceCache {
  struct Slot { Plan plan; std::vector<Scaled> values; std::vector<unsigned char> ready; std::vector<double> arena; };
  std::map<std::vector<long long>, int> index;
  std::vector<std::unique_ptr<Slot>> slots;
  // Owned by this numerical bundle/worker, never shared across workers.
  // Cached outputs own their cells separately; only temporary scratch is reused.
  std::vector<double> shared_arena;
  bool workspace_busy = false;
  struct Workspace {
    PieceCache& owner;bool shared;std::vector<double>& arena;
    Workspace(PieceCache& o,Slot& s):owner(o),shared(s.plan.limits.shared_piece_scratch),arena(shared?o.shared_arena:s.arena){
      if(!shared)return;
      if(owner.workspace_busy||!s.plan.pieces.empty())throw std::logic_error("nested M_ped2 piece workspace");
      owner.workspace_busy=true;
    }
    ~Workspace(){if(shared)owner.workspace_busy=false;}
    Workspace(const Workspace&)=delete;Workspace& operator=(const Workspace&)=delete;
  };
  std::size_t requests = 0, computed = 0;
  double operations = 0.0;
  void clear() { index.clear(); slots.clear(); requests = computed = 0; operations = 0.0; }
  int slot(std::vector<long long> key, Plan&& plan);
  const Scaled& get(int id, int locus, int loci, const Forest& F, const LocusView& L, const LocusMessages& M,
                    const Transmission& T, Context& ctx, PairCache& pairs);
  // Every lane's value of piece id, computed together where possible.
  template <class LaneInput, class LaneTables>
  void get_lanes(int id, int W, const LaneInput* in, const LaneTables& fixed, int loci, const Forest& F, Context& ctx, PairCache& pairs);
};

// A key for everything evaluation reads from a plan, invariant to how the
// plan's genotype variables are numbered and to the order of a birth's two
// parents (transmission is symmetric in them). Variables are ordered by
// colour refinement over the factor graph, ties by their original number;
// the key serializes the relabelled plan, so equal keys are isomorphic plans
// and have equal values in one visit. Gamete variables are internal to their
// birth and do not enter the key.
std::uint64_t mix_hash(std::uint64_t h, std::uint64_t x);
// Set per run only when every locus has Hardy-Weinberg priors consistent with
// its allele frequencies and Mendelian gametes (see reversible_loci).
bool& reversible_keys();

std::vector<long long> plan_key(const Plan& p, std::size_t* cyclic_plans = nullptr);

void assign_gametes(Plan& plan, int count);

} // namespace pap_internal::pedigree_likelihood
#include "barren_variables.h"
namespace pap_internal::pedigree_likelihood {

// Common existing pedigree likelihood piece partition and genotype/allele program backend.
Plan finish_plan(Plan plan,const std::vector<long long>& name,int count,
                       std::size_t model_unaries,std::size_t model_births,PieceCache& pieces,
                       bool split,int max_touch,int focal_var);

// A question owns its endpoints and only their ancestors. Literal hypothesis
// births and latent roles take precedence over current assigned births. Shared
// ancestors are one variable/factor, never independent endpoint marginals.
struct AncestralClosure {
  std::vector<::relationship_topology::TopologyNode> nodes;
  std::vector<unsigned char> active;
};
// Shared structural expansion and cycle check; no numerical plan is required.
AncestralClosure ancestral_closure(const ::relationship_topology::TopologyModel& source,
                                         const Forest& F,const std::vector<int>& roots);

Plan ancestral_plan(const ::relationship_topology::TopologyModel& source,const Forest& F,
                           PieceCache& pieces,const std::vector<int>& roots,bool split,
                           int max_touch,bool retain_focal,bool factorized_one,
                           bool high_allele_programs,ResourceLimits limits);

Plan plan_model(const ::relationship_topology::TopologyModel& source, const Forest& F, const Context& ctx, PairCache& cache,
                       PieceCache& pieces, bool defer_cut_births = false, bool split = true, bool retain_focal = false,
                       int max_touch = 2, bool factorized_one = false, bool high_allele_programs = false, ResourceLimits limits = {},
                       bool local_family_pairs = false, bool latent_brood_exclusion = false,
                       bool ancestral = false,const std::vector<int>& ancestral_roots = {},bool ancestral_reference = false);

// Elimination compiled for one plan and one allele count: table shapes,
// strides and arena offsets are fixed, so a locus only runs arithmetic.
struct Program {
  enum Kind : unsigned char { UnaryTable, GameteTable, PairingTable, OneTable, OneGameteTable, PairTable, PieceTable, Produced };
  struct Slot { Kind kind; int index; std::size_t offset; };
  struct Step {
    int variable=-1;
    std::vector<int> inputs;
    int dim = 0;                           // states of the eliminated variable
    std::vector<int> dims;                 // remaining axes
    std::vector<std::size_t> axis0;        // per input: stride of the eliminated variable
    std::vector<std::size_t> strides;      // per input x remaining axis
    std::size_t cells = 1, offset = 0;
    int out = -1;
    // A pairing indicator input leaves one term of the sum: the eliminated
    // child genotype (role 0) or gamete (roles 1, 2) is read off the others.
    int select = -1, select_role = 0;
  };
  ResourceLimits limits;
  std::shared_ptr<const resource_plan::Execution> exact;
  std::vector<std::vector<int>> audit_scopes;
  std::vector<int> audit_domains;
  // Last attempted step: retained on refusal, without allocating its table.
  int resource_step=-1, eliminated_dim=0, input_count=0, genotype_axes=0, allele_axes=0, selection=-1;
  double required_cells=0, required_arena=0, step_work=0, required_operations=0;
  std::vector<int> output_dims;
  int G = 0, k = 0, width = 0;
  double operations = 0.0;   // multiply-adds per evaluation
  std::vector<Slot> slots;
  std::vector<Step> steps;
  std::vector<int> scalars;
  int result = -1;                         // piece plans: table over the retained variables
  std::vector<std::size_t> unary_offset;   // per variable, or npos
  std::size_t arena = 0;
  live_arena::Stats arena_stats;
  double log_uncovered = 0.0;
  std::uint64_t profile_identity=0;
  int model_variables=0,model_births=0,model_pieces=0,retained_axes=0;
  bool ready = false,coalesced_inputs=false;
  bool search_capped=false;std::size_t search_attempts=0;
  const char* refusal = nullptr;
};

} // namespace pap_internal::pedigree_likelihood
#include "program_profile.h"
namespace pap_internal::pedigree_likelihood {

std::shared_ptr<const Program> compile_scalar_program(const Plan& plan, int G, int k);

} // namespace pap_internal::pedigree_likelihood
#include "resource_plan.h"
namespace pap_internal::pedigree_likelihood {
std::array<double,12>& program_cache_census();
bool& program_cache_census_enabled();
std::shared_ptr<const Program> compile_program(const Plan& plan,int G,int k);

// A program depends only on the plan's shape: its variables and domains and
// the scopes and order of its tables, not on which messages or individuals
// fill them. Compiled programs are shared by shape across questions and
// visits; the registry evicts least-recently-used entries within a byte budget.
struct ShapeHash {
  std::size_t operator()(const std::vector<int>& v) const;
};
#include "program_registry.h"
bool& program_registry_enabled();
std::vector<int> plan_shape(const Plan& plan, int G, int k);

// Resource diagnostic only. Counts and elapsed caller preparation; never an admission input.
std::array<double,12>& program_cache_census();
bool& program_cache_census_enabled();
std::string& program_cache_trace();
const Program& program_for(const Plan& plan, int G, int k);

} // namespace pap_internal::pedigree_likelihood
namespace pap_internal::pedigree_likelihood {

// Runs a compiled program on one locus. Every table is normalized to unit
// maximum and its scale accumulated, so products stay in range.
// Everything outside the plan's components: the stored total minus the
// touched components, plus the split-off family side when untouched.

}
