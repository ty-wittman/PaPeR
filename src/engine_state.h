#pragma once
#include "engine_counters.h"
namespace pap_internal::search {
struct EvaluationState { EngineCounters telemetry; };
struct EngineState : EvaluationState {
  std::vector<pap_internal::pair_rescue::Visit> pair_rescue_visits;
  std::vector<std::array<double,14>> raw_joint_checks;
  std::vector<pap_internal::pair_rescue::Diagnostic> pair_rescue_scores;
  std::vector<pap_internal::pair_rescue::PairDiagnostic> pair_rescue_pairs;
  std::vector<pap_internal::pair_rescue::Requirement> pair_rescue_requirements;
  int n = 0;
  int loci = 0;
  int message_revision = 0;
  std::vector<LocusTables> locus;
  std::vector<int> genotype_alleles;
  double bytes = 0.0;
};

}
