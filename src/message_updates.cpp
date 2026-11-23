#include "message_updates.h"
#include "scoring.h"
namespace pap_internal::search {
static std::vector<int> topological_order(const VisitGraph& pass);
static std::vector<double> transmit_alleles(
    const LocusTables& locus, const std::vector<double>& message);
static std::vector<double> child_from_alleles(
    const LocusTables& locus, const std::vector<double>& dam,
    const std::vector<double>& sire);
static std::vector<double> parent_joint_native(
    const LocusTables& locus,
    const std::vector<double>& first, const std::vector<double>& second,
    const IbdClass& relation);
static std::vector<double> parent_joint(
    VisitGraph& pass, const LocusTables& locus,
    const std::vector<double>& first, const std::vector<double>& second,
    const IbdClass& relation);
static std::vector<double> child_from_joint(
    const LocusTables& locus, const std::vector<double>& joint);
static MessagePublication publish_messages_full_serial(
    EngineState& compiled, const IntegerMatrix& parents,
    int graph_revision);
static void add_signature_contribution(
    const LocusTables& locus, std::size_t cell_offset, int child,
    double& total, double& squared, double& indexed);
static std::string signature_string(
    std::size_t count, double total, double squared, double indexed);
static void recompute_signature(
    const EngineState& compiled, std::size_t& count,
    double& total, double& squared, double& indexed);
static List active_candidate_side(
    const CandidateSide& source, const IntegerMatrix& parents,
    int parent_column);

static std::vector<int> topological_order(const VisitGraph& pass) {
  const int n = pass.compiled.n;
  std::vector<int> indegree(static_cast<std::size_t>(n + 1), 0);
  std::queue<int> queue;
  for (int child = 1; child <= n; ++child) {
    indegree[child] = (pass.dam[child] > 0 ? 1 : 0) +
      (pass.sire[child] > 0 ? 1 : 0);
    if (indegree[child] == 0) queue.push(child);
  }
  std::vector<int> order;
  order.reserve(static_cast<std::size_t>(n));
  while (!queue.empty()) {
    const int node = queue.front();
    queue.pop();
    order.push_back(node);
    for (int child : pass.reverse[node])
      if (--indegree[child] == 0) queue.push(child);
  }
  if (static_cast<int>(order.size()) != n)
    engine_stop("native Stage-1 publication requires an acyclic parent graph");
  return order;
}

static std::vector<double> transmit_alleles(
    const LocusTables& locus, const std::vector<double>& message) {
  std::vector<double> out(static_cast<std::size_t>(locus.k), R_NegInf);
  std::vector<double> terms(static_cast<std::size_t>(locus.G));
  for (int allele = 0; allele < locus.k; ++allele) {
    for (int genotype = 0; genotype < locus.G; ++genotype)
      terms[genotype] = message[genotype] +
        locus.transmission_ln[genotype + allele * locus.G];
    out[allele] = log_sum(terms, false, "native Stage-1 transmitted message");
  }
  return out;
}

static std::vector<double> child_from_alleles(
    const LocusTables& locus, const std::vector<double>& dam,
    const std::vector<double>& sire) {
  std::vector<double> out(static_cast<std::size_t>(locus.G), R_NegInf);
  for (int a = 0; a < locus.k; ++a) for (int b = a; b < locus.k; ++b) {
    const int genotype = genotype_index_zero(a, b, locus.k);
    out[genotype] = dam[a] + sire[b];
    if (a != b) out[genotype] = log_add(out[genotype], dam[b] + sire[a]);
  }
  return out;
}

static std::vector<double> parent_joint_native(
    const LocusTables& locus,
    const std::vector<double>& first, const std::vector<double>& second,
    const IbdClass& relation) {
  const std::size_t cells = static_cast<std::size_t>(locus.G) * locus.G;
  std::vector<double> probability(cells, 0.0);
  std::vector<double> prior(static_cast<std::size_t>(locus.G));
  std::vector<double> allele(static_cast<std::size_t>(locus.k));
  for (int g = 0; g < locus.G; ++g) prior[g] = std::exp(locus.prior_ln[g]);
  for (int a = 0; a < locus.k; ++a) allele[a] = std::exp(locus.log_p[a]);
  for (int sire = 0; sire < locus.G; ++sire)
    for (int dam = 0; dam < locus.G; ++dam)
      probability[dam + sire * locus.G] =
        relation.k0 * prior[dam] * prior[sire];
  if (relation.k1 > 0.0)
    for (int shared = 0; shared < locus.k; ++shared)
      for (int left = 0; left < locus.k; ++left)
        for (int right = 0; right < locus.k; ++right) {
          const int g1 = genotype_index_zero(shared, left, locus.k);
          const int g2 = genotype_index_zero(shared, right, locus.k);
          probability[g1 + g2 * locus.G] += relation.k1 * allele[shared] *
            allele[left] * allele[right];
        }
  if (relation.k2 > 0.0)
    for (int g = 0; g < locus.G; ++g)
      probability[g + g * locus.G] += relation.k2 * prior[g];
  std::vector<double> out(cells, R_NegInf);
  for (int sire = 0; sire < locus.G; ++sire)
    for (int dam = 0; dam < locus.G; ++dam) {
      const double mass = probability[dam + sire * locus.G];
      if (mass > 0.0) out[dam + sire * locus.G] = std::log(mass) +
        first[dam] - locus.prior_ln[dam] +
        second[sire] - locus.prior_ln[sire];
  }
  normalize_in_place(out);
  return out;
}

static std::vector<double> parent_joint(
    VisitGraph& pass, const LocusTables& locus,
    const std::vector<double>& first, const std::vector<double>& second,
    const IbdClass& relation) {
  std::vector<double> out = parent_joint_native(
    locus, first, second, relation);

  return out;
}

static std::vector<double> child_from_joint(
    const LocusTables& locus, const std::vector<double>& joint) {
  std::vector<double> allele_joint(
    static_cast<std::size_t>(locus.k) * locus.k, R_NegInf);
  std::vector<double> terms;
  terms.reserve(static_cast<std::size_t>(locus.G) * locus.G);
  for (int sire_allele = 0; sire_allele < locus.k; ++sire_allele)
    for (int dam_allele = 0; dam_allele < locus.k; ++dam_allele) {
      terms.clear();
      for (int sire = 0; sire < locus.G; ++sire)
        for (int dam = 0; dam < locus.G; ++dam)
          terms.push_back(joint[dam + sire * locus.G] +
            locus.transmission_ln[dam + dam_allele * locus.G] +
            locus.transmission_ln[sire + sire_allele * locus.G]);
      allele_joint[dam_allele + sire_allele * locus.k] = log_sum(
        terms, false, "native Stage-1 parent-joint child contraction");
    }
  std::vector<double> out(static_cast<std::size_t>(locus.G), R_NegInf);
  for (int a = 0; a < locus.k; ++a) for (int b = a; b < locus.k; ++b) {
    const int genotype = genotype_index_zero(a, b, locus.k);
    out[genotype] = allele_joint[a + b * locus.k];
    if (a != b) out[genotype] = log_add(
      out[genotype], allele_joint[b + a * locus.k]);
  }
  return out;
}

static MessagePublication publish_messages_full_serial(
    EngineState& compiled, const IntegerMatrix& parents,
    int graph_revision) {
  const auto started = std::chrono::steady_clock::now();
  VisitGraph pass(compiled, 1, graph_revision, parents);
  pass.bounded_messages = true;
  const std::vector<int> order = topological_order(pass);
  MessagePublication publication;
  publication.total_nodes = static_cast<std::size_t>(compiled.n);
  publication.dirty_nodes = static_cast<std::size_t>(compiled.n);
  publication.nodes_visited = static_cast<std::size_t>(compiled.n) *
    static_cast<std::size_t>(compiled.loci);
  publication.nodes_recomputed = publication.nodes_visited;
  publication.message_vectors_rebuilt = publication.nodes_visited;
  publication.likelihood_cells = publication.nodes_visited;
  publication.node_likelihood_ln.assign(
    static_cast<std::size_t>(compiled.n), 0.0);
  publication.locus_likelihood_ln.assign(
    static_cast<std::size_t>(compiled.loci) * compiled.n, 0.0);
  std::size_t count = 0;
  double total = 0.0, squared = 0.0, indexed = 0.0;
  for (int locus_index = 0; locus_index < compiled.loci; ++locus_index) {
    LocusTables& locus = compiled.locus[static_cast<std::size_t>(locus_index)];
    std::vector<std::vector<double>> messages(
      static_cast<std::size_t>(compiled.n + 1));
    for (int child : order) {
      const int dam = pass.dam[child];
      const int sire = pass.sire[child];
      std::vector<double> inherited;
      if (dam > 0 && sire > 0) {
        const IbdClass relation = pass.ibd(dam, sire);
        if (relation.k1 > 0.0 || relation.k2 > 0.0) {
          inherited = child_from_joint(
            locus, parent_joint(pass, locus, messages[dam],
                                       messages[sire], relation));
        } else inherited = child_from_alleles(
          locus, transmit_alleles(locus, messages[dam]),
          transmit_alleles(locus, messages[sire]));
      } else {
        const std::vector<double> dam_allele = dam > 0 ?
          transmit_alleles(locus, messages[dam]) : locus.log_p;
        const std::vector<double> sire_allele = sire > 0 ?
          transmit_alleles(locus, messages[sire]) : locus.log_p;
        inherited = child_from_alleles(locus, dam_allele, sire_allele);
      }
      std::vector<double>& message = messages[child];
      message.resize(static_cast<std::size_t>(locus.G));
      for (int genotype = 0; genotype < locus.G; ++genotype)
        message[genotype] = inherited[genotype] +
          locus.child_emission(child - 1, genotype);
      const double likelihood = log_sum(
        message, false, "native Stage-1 individual likelihood");
      if (!R_finite(likelihood))
        engine_stop("native Stage-1 publication produced an invalid individual likelihood");
      publication.node_likelihood_ln[child - 1] += likelihood;
      publication.locus_likelihood_ln[
        static_cast<std::size_t>(locus_index) * compiled.n + child - 1] =
          likelihood;
      normalize_in_place(message);
      for (int genotype = 0; genotype < locus.G; ++genotype)
        locus.message_ln[genotype + (child - 1) * locus.G] = message[genotype];
      publication.genotype_state_elements +=
        static_cast<std::size_t>(locus.G);
      const std::vector<double> transmitted = transmit_alleles(locus, message);
      for (int allele = 0; allele < locus.k; ++allele)
        locus.transmitted_message_ln[(child - 1) * locus.k + allele] =
          transmitted[allele];
      publication.transmitted_state_elements +=
        static_cast<std::size_t>(locus.k);
    }
    for (int child = 1; child <= compiled.n; ++child) {
      long double message_total = 0.0L;
      long double message_squared = 0.0L;
      long double message_indexed = 0.0L;
      for (int genotype = 0; genotype < locus.G; ++genotype) {
        const double value = locus.message_ln[
          genotype + (child - 1) * locus.G];
        message_total += static_cast<long double>(value);
        message_squared += static_cast<long double>(value) * value;
        message_indexed += static_cast<long double>(value) *
          static_cast<long double>(count + genotype + 1);
      }
      total += static_cast<double>(message_total);
      squared += static_cast<double>(message_squared);
      indexed += static_cast<double>(message_indexed);
      count += static_cast<std::size_t>(locus.G);
    }
  }
  publication.total_likelihood_ln = std::accumulate(
    publication.node_likelihood_ln.begin(),
    publication.node_likelihood_ln.end(), 0.0);
  char signature[256];
  std::snprintf(signature, sizeof(signature),
                "%zu %.17g:%zu %.17g:%zu %.17g",
                count, total, count, squared, count, indexed);
  publication.evidence_signature = signature;
  publication.signature_count = count;
  publication.signature_total = total;
  publication.signature_squared = squared;
  publication.signature_indexed = indexed;
  publication.aggregate_updates = static_cast<std::size_t>(compiled.n) *
    static_cast<std::size_t>(compiled.loci);
  publication.signature_global_scans = 1;
  publication.allocations = static_cast<std::size_t>(compiled.loci) *
    (static_cast<std::size_t>(compiled.n) + 1);
  publication.seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - started).count();
  return publication;
}

static void add_signature_contribution(
    const LocusTables& locus, std::size_t cell_offset, int child,
    double& total, double& squared, double& indexed) {
  long double message_total = 0.0L;
  long double message_squared = 0.0L;
  long double message_indexed = 0.0L;
  const std::size_t offset = cell_offset +
    static_cast<std::size_t>(child - 1) * static_cast<std::size_t>(locus.G);
  for (int genotype = 0; genotype < locus.G; ++genotype) {
    const double value = locus.message_ln[
      genotype + (child - 1) * locus.G];
    message_total += static_cast<long double>(value);
    message_squared += static_cast<long double>(value) * value;
    message_indexed += static_cast<long double>(value) *
      static_cast<long double>(offset + static_cast<std::size_t>(genotype) + 1);
  }
  total += static_cast<double>(message_total);
  squared += static_cast<double>(message_squared);
  indexed += static_cast<double>(message_indexed);
}

static std::string signature_string(
    std::size_t count, double total, double squared, double indexed) {
  char signature[256];
  std::snprintf(signature, sizeof(signature),
                "%zu %.17g:%zu %.17g:%zu %.17g",
                count, total, count, squared, count, indexed);
  return signature;
}

static void recompute_signature(
    const EngineState& compiled, std::size_t& count,
    double& total, double& squared, double& indexed) {
  count = 0;
  total = 0.0;
  squared = 0.0;
  indexed = 0.0;
  std::size_t cell_offset = 0;
  for (int locus_index = 0; locus_index < compiled.loci; ++locus_index) {
    const LocusTables& locus = compiled.locus[
      static_cast<std::size_t>(locus_index)];
    for (int child = 1; child <= compiled.n; ++child) {
      double child_total = 0.0, child_squared = 0.0, child_indexed = 0.0;
      add_signature_contribution(
        locus, cell_offset, child,
        child_total, child_squared, child_indexed);
      total += child_total;
      squared += child_squared;
      indexed += child_indexed;
      count += static_cast<std::size_t>(locus.G);
    }
    cell_offset += static_cast<std::size_t>(compiled.n) *
      static_cast<std::size_t>(locus.G);
  }
}

MessagePublication publish_messages_full(
    EngineState& compiled, const IntegerMatrix& parents,
    int graph_revision) {
  {
    MessagePublication publication = publish_messages_full_serial(
      compiled, parents, graph_revision);
    publication.requested_workers = 1;
    publication.effective_publication_workers = 1;
    publication.worker_thread_count = 0;
    publication.publication_job_count = 0;
    publication.worker_r_api_calls = 0;
    publication.worker_seconds.assign(1, publication.seconds);
    publication.locus_owner.assign(
      static_cast<std::size_t>(compiled.loci), 1);
    return publication;
  }

}

void publish_messages_selective(
    EngineState& compiled, const IntegerMatrix& parents,
    int graph_revision, MessagePublication& publication,
    const DirtyRegion& dirty_region) {
  const auto started = std::chrono::steady_clock::now();
  if(pap_internal::input_checks::enabled()){if (dirty_region.graph_revision_after != graph_revision)engine_stop("selective Stage-1 dirty region has the wrong graph revision");}
  const std::vector<int>& dirty_nodes = dirty_region.publication_nodes;
  if (publication.node_likelihood_ln.size() !=
        static_cast<std::size_t>(compiled.n) ||
      publication.locus_likelihood_ln.size() !=
        static_cast<std::size_t>(compiled.n) * compiled.loci)
    engine_stop("selective Stage-1 publication lacks a complete prior publication");
  if (!std::is_sorted(dirty_nodes.begin(), dirty_nodes.end()))
    engine_stop("selective Stage-1 dirty node IDs must be sorted");
  for (int node : dirty_nodes)
    if (node < 1 || node > compiled.n)
      engine_stop("selective Stage-1 dirty node ID is out of range");

  VisitGraph pass(compiled, 1, graph_revision, parents);
  pass.bounded_messages = true;
  const std::vector<int> order = topological_order(pass);
  const std::size_t dirty_count = dirty_nodes.size();
  std::vector<double> old_node_likelihood(dirty_count);
  for (std::size_t i = 0; i < dirty_count; ++i) {
    const int child = dirty_nodes[i];
    old_node_likelihood[i] = publication.node_likelihood_ln[child - 1];
  }

  publication.selective = true;
  publication.dirty_nodes = dirty_count;
  publication.total_nodes = static_cast<std::size_t>(compiled.n);
  publication.nodes_visited = dirty_count *
    static_cast<std::size_t>(compiled.loci);
  publication.nodes_recomputed = publication.nodes_visited;
  publication.message_vectors_rebuilt = publication.nodes_visited;
  publication.likelihood_cells = publication.nodes_visited;
  publication.genotype_state_elements = 0;
  publication.transmitted_state_elements = 0;
  publication.aggregate_updates = dirty_count;
  publication.signature_global_scans = 0;
  publication.copies = 0;
  publication.allocations = 0;

  for (int locus_index = 0; locus_index < compiled.loci; ++locus_index) {
    LocusTables& locus = compiled.locus[static_cast<std::size_t>(locus_index)];
    for (int child : order) {
      if (!std::binary_search(dirty_nodes.begin(), dirty_nodes.end(), child))
        continue;
      const int rank = static_cast<int>(std::lower_bound(
        dirty_nodes.begin(), dirty_nodes.end(), child) - dirty_nodes.begin());
      const int dam = pass.dam[child];
      const int sire = pass.sire[child];
      std::vector<double> inherited;
      if (dam > 0 && sire > 0) {
        const IbdClass relation = pass.ibd(dam, sire);
        if (relation.k1 > 0.0 || relation.k2 > 0.0) {
          std::vector<double> dam_message(
            locus.message_ln.begin() + (dam - 1) * locus.G,
            locus.message_ln.begin() + dam * locus.G);
          std::vector<double> sire_message(
            locus.message_ln.begin() + (sire - 1) * locus.G,
            locus.message_ln.begin() + sire * locus.G);
          publication.copies += static_cast<std::size_t>(2 * locus.G);
          inherited = child_from_joint(
            locus, parent_joint(pass, locus, dam_message,
                                       sire_message, relation));
        } else {
          std::vector<double> dam_transmitted(
            locus.transmitted_message_ln.begin() + (dam - 1) * locus.k,
            locus.transmitted_message_ln.begin() + dam * locus.k);
          std::vector<double> sire_transmitted(
            locus.transmitted_message_ln.begin() + (sire - 1) * locus.k,
            locus.transmitted_message_ln.begin() + sire * locus.k);
          publication.copies += static_cast<std::size_t>(2 * locus.k);
          inherited = child_from_alleles(
            locus, dam_transmitted, sire_transmitted);
        }
      } else {
        const std::vector<double> dam_allele = dam > 0 ?
          std::vector<double>(
            locus.transmitted_message_ln.begin() + (dam - 1) * locus.k,
            locus.transmitted_message_ln.begin() + dam * locus.k) :
          locus.log_p;
        const std::vector<double> sire_allele = sire > 0 ?
          std::vector<double>(
            locus.transmitted_message_ln.begin() + (sire - 1) * locus.k,
            locus.transmitted_message_ln.begin() + sire * locus.k) :
          locus.log_p;
        publication.copies += static_cast<std::size_t>(
          (dam > 0 ? locus.k : 0) + (sire > 0 ? locus.k : 0));
        inherited = child_from_alleles(
          locus, dam_allele, sire_allele);
      }
      std::vector<double> message(static_cast<std::size_t>(locus.G));
      ++publication.allocations;
      for (int genotype = 0; genotype < locus.G; ++genotype)
        message[genotype] = inherited[genotype] +
          locus.child_emission(child - 1, genotype);
      const double likelihood = log_sum(
        message, false, "native Stage-1 selective individual likelihood");
      if (!R_finite(likelihood))
        engine_stop("native Stage-1 selective publication produced an invalid individual likelihood");
      const std::size_t likelihood_index =
        static_cast<std::size_t>(locus_index) * compiled.n + child - 1;
      const double old_locus_likelihood =
        publication.locus_likelihood_ln[likelihood_index];
      publication.locus_likelihood_ln[likelihood_index] = likelihood;
      publication.node_likelihood_ln[child - 1] +=
        likelihood - old_locus_likelihood;
      normalize_in_place(message);
      for (int genotype = 0; genotype < locus.G; ++genotype)
        locus.message_ln[genotype + (child - 1) * locus.G] = message[genotype];
      publication.genotype_state_elements += static_cast<std::size_t>(locus.G);
      const std::vector<double> transmitted = transmit_alleles(locus, message);
      ++publication.allocations;
      for (int allele = 0; allele < locus.k; ++allele)
        locus.transmitted_message_ln[(child - 1) * locus.k + allele] =
          transmitted[allele];
      publication.transmitted_state_elements += static_cast<std::size_t>(locus.k);
      (void)rank;
    }
  }

  double signature_total = publication.signature_total;
  double signature_squared = publication.signature_squared;
  double signature_indexed = publication.signature_indexed;
  recompute_signature(
    compiled, publication.signature_count, signature_total,
    signature_squared, signature_indexed);
  publication.signature_global_scans = 1;
  publication.aggregate_updates += static_cast<std::size_t>(compiled.n) *
    static_cast<std::size_t>(compiled.loci);
  for (std::size_t i = 0; i < dirty_count; ++i)
    publication.total_likelihood_ln +=
      publication.node_likelihood_ln[dirty_nodes[i] - 1] -
      old_node_likelihood[i];
  publication.signature_total = signature_total;
  publication.signature_squared = signature_squared;
  publication.signature_indexed = signature_indexed;
  publication.evidence_signature = signature_string(
    publication.signature_count, signature_total, signature_squared,
    signature_indexed);

  // Cache invalidation is a graph-revision side effect of every full
  // publication and is intentionally retained exactly on the selective path.
  publication.seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - started).count();
}

IntegerVector focal_order(
    const IntegerMatrix& parents, const NumericVector& expected_year,
    const IntegerVector& id_order_rank) {
  const int n = parents.nrow();
  std::vector<int> generation(static_cast<std::size_t>(n + 1), 0);
  EngineState dummy;
  dummy.n = n;
  VisitGraph pass(dummy, 1, 0, parents);
  const std::vector<int> topo = topological_order(pass);
  for (int child : topo) {
    const int dam = parents(child - 1, 0);
    const int sire = parents(child - 1, 1);
    if (dam > 0) generation[child] = std::max(
      generation[child], generation[dam] + 1);
    if (sire > 0) generation[child] = std::max(
      generation[child], generation[sire] + 1);
  }
  std::vector<int> order(static_cast<std::size_t>(n));
  std::iota(order.begin(), order.end(), 1);
  std::stable_sort(order.begin(), order.end(), [&](int left, int right) {
    const double ly = expected_year[left - 1];
    const double ry = expected_year[right - 1];
    if (ly != ry) return ly < ry;
    if (generation[left] != generation[right])
      return generation[left] < generation[right];
    return id_order_rank[left - 1] < id_order_rank[right - 1];
  });
  return wrap(order);
}

std::string parent_signature(
    const CharacterVector& ids, const IntegerMatrix& parents) {
  std::ostringstream out;
  out << "ids=";
  for (int i = 0; i < ids.size(); ++i) {
    if (i) out << ',';
    out << as<std::string>(ids[i]);
  }
  out << ";edges=";
  bool first = true;
  for (int side = 0; side < 2; ++side) for (int child = 0;
       child < parents.nrow(); ++child) if (parents(child, side) > 0) {
    if (!first) out << ';';
    first = false;
    out << as<std::string>(ids[child]) << '|' <<
      (side == 0 ? "dam" : "sire") << '|' <<
      as<std::string>(ids[parents(child, side) - 1]);
  }
  return out.str();
}

static List active_candidate_side(
    const CandidateSide& source, const IntegerMatrix& parents,
    int parent_column) {
  const int n = parents.nrow();
  std::vector<int> offsets(static_cast<std::size_t>(n + 1), 0);
  std::vector<int> ids;
  std::vector<double> quick;
  std::vector<int> forced(static_cast<std::size_t>(n), 0);
  ids.reserve(source.ids.size());
  quick.reserve(source.quick.size());
  for (int child = 0; child < n; ++child) {
    offsets[child] = static_cast<int>(ids.size());
    const int begin = source.offsets[child];
    const int end = source.offsets[child + 1];
    for (int at = begin; at < end; ++at) if (source.permanent[at]) {
      ids.push_back(source.ids[at]);
      quick.push_back(source.quick[at]);
    }
    const int incumbent = parents(child, parent_column);
    if (incumbent <= 0) continue;
    const bool retained = std::find(ids.begin() + offsets[child], ids.end(),
                                    incumbent) != ids.end();
    if (retained) continue;
    int incumbent_at = -1;
    for (int at = begin; at < end; ++at) if (source.ids[at] == incumbent) {
      incumbent_at = at;
      break;
    }
    if (incumbent_at < 0)
      engine_stop("native Stage-1 current incumbent lacks static candidate metadata");
    ids.push_back(incumbent);
    quick.push_back(source.quick[incumbent_at]);
    forced[child] = 1;
  }
  offsets[n] = static_cast<int>(ids.size());
  return List::create(
    _["offsets"] = offsets, _["ids"] = ids,
    _["quick"] = quick,
    _["forced"] = forced);
}

List refresh_active_candidates(
    List& inputs, const IntegerMatrix& parents,
    const CandidateSide& dam_source,
    const CandidateSide& sire_source) {
  const List dam = active_candidate_side(dam_source, parents, 0);
  const List sire = active_candidate_side(sire_source, parents, 1);
  inputs["dam_offsets"] = dam["offsets"];
  inputs["dam_ids"] = dam["ids"];
  inputs["dam_quick"] = dam["quick"];
  inputs["sire_offsets"] = sire["offsets"];
  inputs["sire_ids"] = sire["ids"];
  inputs["sire_quick"] = sire["quick"];
  inputs["forced"] = List::create(
    _["dam"] = dam["forced"], _["sire"] = sire["forced"]);
  return List::create(
    _["dam_offsets"] = dam["offsets"], _["dam_ids"] = dam["ids"],
    _["dam_quick"] = dam["quick"],
    _["sire_offsets"] = sire["offsets"], _["sire_ids"] = sire["ids"],
    _["sire_quick"] = sire["quick"],
    _["forced"] = inputs["forced"]);
}

List message_matrices(const EngineState& compiled) {
  List messages(compiled.loci);
  for (int locus_index = 0; locus_index < compiled.loci; ++locus_index) {
    const LocusTables& locus = compiled.locus[
      static_cast<std::size_t>(locus_index)];
    NumericMatrix message(locus.G, compiled.n);
    std::copy(locus.message_ln.begin(), locus.message_ln.end(), message.begin());
    messages[locus_index] = message;
  }
  return messages;
}

NumericMatrix locus_likelihood_matrix(
    const MessagePublication& publication, int loci, int n) {
  NumericMatrix out(loci, n);
  for (int locus = 0; locus < loci; ++locus)
    for (int child = 0; child < n; ++child)
      out(locus, child) = publication.locus_likelihood_ln[
        static_cast<std::size_t>(locus) * n + child];
  return out;
}

List publication_work_list(const MessagePublication& publication) {
  return List::create(
    _["mode"] = publication.selective ? "selective" : "full",
    _["dirty_nodes"] = static_cast<double>(publication.dirty_nodes),
    _["total_nodes"] = static_cast<double>(publication.total_nodes),
    _["nodes_visited"] = static_cast<double>(publication.nodes_visited),
    _["nodes_recomputed"] = static_cast<double>(publication.nodes_recomputed),
    _["message_vectors_rebuilt"] = static_cast<double>(
      publication.message_vectors_rebuilt),
    _["genotype_state_elements"] = static_cast<double>(
      publication.genotype_state_elements),
    _["transmitted_state_elements"] = static_cast<double>(
      publication.transmitted_state_elements),
    _["likelihood_cells"] = static_cast<double>(publication.likelihood_cells),
    _["aggregate_updates"] = static_cast<double>(
      publication.aggregate_updates),
    _["signature_global_scans"] = static_cast<double>(
      publication.signature_global_scans),
    _["copies"] = static_cast<double>(publication.copies),
    _["allocations"] = static_cast<double>(publication.allocations),
    _["requested_workers"] = publication.requested_workers,
    _["effective_publication_workers"] =
      publication.effective_publication_workers,
    _["worker_thread_count"] = publication.worker_thread_count,
    _["publication_job_count"] = publication.publication_job_count,
    _["worker_r_api_calls"] = publication.worker_r_api_calls,
    _["worker_rcpp_api_calls"] = static_cast<int>(0),
    _["worker_seconds"] = publication.worker_seconds,
    _["locus_owner"] = publication.locus_owner,
    _["assignment_policy"] =
      "deterministic_lpt_G2_descending_original_locus_index_lightest_worker",
    _["reduction_policy"] = "serial_original_locus_order",
    _["finalization_seconds"] = publication.finalization_seconds,
    _["dispatch_seconds"] = publication.dispatch_seconds,
    _["publication_seconds"] = publication.seconds);
}

CandidateSide::CandidateSide(const List& inputs, const char* prefix, int n) {
    const std::string root(prefix);
    const IntegerVector offsets_r = inputs[root + "_offsets"];
    const IntegerVector ids_r = inputs[root + "_ids"];
    const NumericVector quick_r = inputs[root + "_quick"];
    const LogicalVector permanent_r = inputs[root + "_permanent"];
    offsets.assign(offsets_r.begin(), offsets_r.end());
    ids.assign(ids_r.begin(), ids_r.end());
    quick.assign(quick_r.begin(), quick_r.end());
    permanent.reserve(static_cast<std::size_t>(permanent_r.size()));
    for (int value : permanent_r) permanent.push_back(value == TRUE ? 1 : 0);
    if (static_cast<int>(offsets.size()) != n + 1 ||
        ids.size() != quick.size() ||
        ids.size() != permanent.size() || offsets.front() != 0 ||
        offsets.back() != static_cast<int>(ids.size()))
      engine_stop("native Stage-1 static candidate side is malformed");
  }

}
