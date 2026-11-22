#include "engine_state.h"
namespace pap_internal::search {

}

#include "genotype_checks.h"
using namespace pap_internal::search;
// [[Rcpp::export]]
SEXP build_snapshot_cpp(List inputs) {
  const auto started = std::chrono::steady_clock::now();
  if (!inputs.containsElementNamed("n_individuals") ||
      !inputs.containsElementNamed("message_revision"))
    engine_stop("simp snapshot inputs are missing dimensions or message revision");
  const int n = as<int>(inputs["n_individuals"]);
  const int message_revision = as<int>(inputs["message_revision"]);
  if (n < 1) engine_stop("simp snapshot requires at least one individual");
  const List prior_list = inputs["prior_ln"];
  const List log_p_list = inputs["log_p"];
  const List transmission_list = inputs["transmission_log_weight"];
  const List child_index_list = inputs["child_index"];
  const List observation_list = inputs["observations_ln"];
  const List message_list = inputs["messages_ln"];
  const bool captured_transmitted = inputs.containsElementNamed("transmitted_messages_ln");
  const List transmitted_list = captured_transmitted ?
    as<List>(inputs["transmitted_messages_ln"]) : List();
  const int loci = prior_list.size();
  if (loci < 1 || log_p_list.size() != loci || transmission_list.size() != loci ||
      child_index_list.size() != loci || observation_list.size() != loci ||
      message_list.size() != loci ||
      (captured_transmitted && transmitted_list.size() != loci))
    engine_stop("simp snapshot locus inputs must have equal nonzero lengths");
  XPtr<EngineState> pointer(new EngineState(), true);
  EngineState& compiled = *pointer;
  compiled.n = n;
  compiled.loci = loci;
  compiled.message_revision = message_revision;
  compiled.locus.reserve(loci);
  if (inputs.containsElementNamed("genotype_alleles")) {
    const IntegerMatrix genotype_alleles = inputs["genotype_alleles"];
    if (genotype_alleles.nrow() != n || genotype_alleles.ncol() != 2 * loci)
      engine_stop("simp snapshot genotype alleles have incompatible dimensions");
    compiled.genotype_alleles.assign(genotype_alleles.begin(), genotype_alleles.end());
  }
  for (int l = 0; l < loci; ++l) {
    NumericVector prior = prior_list[l];
    NumericVector log_p = log_p_list[l];
    NumericMatrix transmission = transmission_list[l];
    IntegerMatrix child_index = child_index_list[l];
    NumericMatrix observations = observation_list[l];
    NumericMatrix messages = message_list[l];
    const int k = log_p.size();
    const int G = prior.size();
    if (k < 1 || G != k * (k + 1) / 2)
      engine_stop("simp snapshot locus has invalid k/G dimensions");
    if (transmission.nrow() != G || transmission.ncol() != k)
      engine_stop("simp snapshot transmission dimensions are invalid");
    if (child_index.nrow() != k || child_index.ncol() != k)
      engine_stop("simp snapshot child-index dimensions are invalid");
    validate_locus_matrix(observations, G, n, "observations");
    validate_locus_matrix(messages, G, n, "messages");
    LocusTables locus;
    locus.k = k; locus.G = G;
    locus.prior_ln.assign(prior.begin(), prior.end());
    locus.log_p.assign(log_p.begin(), log_p.end());
    locus.transmission_ln.assign(transmission.begin(), transmission.end());
    locus.child_index.assign(child_index.begin(), child_index.end());
    locus.message_ln.assign(messages.begin(), messages.end());
    locus.genotype_first.resize(static_cast<std::size_t>(G));
    locus.genotype_second.resize(static_cast<std::size_t>(G));
    for (int a = 0; a < k; ++a) for (int b = a; b < k; ++b) {
      const int expected = genotype_index_zero(a, b, k) + 1;
      if (child_index(a, b) != expected)
        engine_stop("simp snapshot child-index values must use the canonical one-based order");
      const int genotype = expected - 1;
      locus.genotype_first[genotype] = a;
      locus.genotype_second[genotype] = b;
    }
    locus.child_kernel_ln.resize(static_cast<std::size_t>(n) * k * k, R_NegInf);
    for (int child = 0; child < n; ++child)
      for (int a = 0; a < k; ++a) for (int b = 0; b < k; ++b) {
        const int genotype = child_index(std::min(a, b), std::max(a, b)) - 1;
        locus.child_kernel_ln[(child * k + a) * k + b] =
          observations[genotype + child * G];
      }
    locus.transmitted_message_ln.resize(static_cast<std::size_t>(n) * k, R_NegInf);
    for (int individual = 0; individual < n; ++individual)
      for (int allele = 0; allele < k; ++allele) {
        LogSumExpAccumulator accumulator;
        for (int genotype = 0; genotype < G; ++genotype)
          accumulator.add_log(
            locus.message_ln[genotype + individual * G] +
              locus.transmission_ln[genotype + allele * G],
            "simp snapshot transmitted message");
        locus.transmitted_message_ln[individual * k + allele] = accumulator.value();
      }
    if (captured_transmitted) {
      const NumericMatrix transmitted = transmitted_list[l];
      validate_locus_matrix(transmitted, k, n, "captured transmitted messages");
      locus.transmitted_message_ln.assign(transmitted.begin(), transmitted.end());
    }
    compiled.locus.push_back(std::move(locus));
  }
  std::size_t bytes = sizeof(EngineState);
  for (const LocusTables& locus : compiled.locus) {
    bytes += sizeof(LocusTables);
    bytes += sizeof(double) * (locus.prior_ln.size() + locus.log_p.size() +
      locus.transmission_ln.size() + locus.message_ln.size() +
      locus.child_kernel_ln.size() +
      locus.transmitted_message_ln.size());
    bytes += sizeof(int) * (locus.child_index.size() +
      locus.genotype_first.size() + locus.genotype_second.size());
  }
  bytes += sizeof(int) * compiled.genotype_alleles.size();
  compiled.bytes = static_cast<double>(bytes);

  pointer.attr("search_snapshot_bytes") = compiled.bytes;
  pointer.attr("search_snapshot_message_revision") = compiled.message_revision;
  return pointer;
}
