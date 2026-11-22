#pragma once
#include "score_records.h"
namespace pap_internal::search {
struct VisitScorer {
  double partition_memo_pass_peak_bytes = 0.0;
  double partition_memo_pass_peak_entries = 0.0;
  int requested_workers = 1;
  int effective_workers = 1;
  int worker_thread_count = 0;
  std::vector<std::unique_ptr<LocusContext>> contexts;
  std::vector<std::vector<int>> assignments;
  std::vector<int> owner_by_locus;
  std::vector<int> local_position_by_locus;
  std::vector<double> weight_by_locus;

 private:
  std::atomic<bool> unavailable_message_trial{false};


 public:
  VisitScorer(EngineState& compiled,
                       const SearchState& authoritative,
                       const SearchInputs& input);

  LocusContext& primary() { return *contexts.front(); }
  const LocusContext& primary() const { return *contexts.front(); }


  std::size_t result_metadata_bytes() const;


  void capture_result_storage(const SearchInputs& input, bool completed_wave = false);


  // Raw and relationship phases share the existing per-locus Q cache. preparation
  // keeps one result row; Wn reuses its already-owned slots. Only this caller
  // reduces scientific totals, always in original marker order.
  void score_raw_panel(const SearchInputs& input, int child,
      const std::vector<RelationshipContexts>& side_contexts,
      std::vector<ConfigurationScore>& configurations,
      const std::function<void(int)>& progress = {});

  // Stage two preserves the original locus ownership and strict ascending
  // reduction; the same bounded result wave is reused after stage one ends.


  // The raw-Q kernel may change operation order. Refine complete retained
  // panels before ANY raw support certificate or identity/chooser use, for
  // every subset of competitors that later relationship gates might retain.
  void refine_raw_panels(const SearchInputs& input,
      const std::vector<RelationshipContexts>& side_contexts,
      std::vector<ConfigurationScore>& configurations);

  void configure_raw_arithmetic(bool fast);

  void configure_message_scoring(bool selected, bool published);

  // Check one genuinely required raw question before speculative relationship
  // panels. The largest locus is a work-order hint only: backend selection
  // requires the original scorer's actual unavailable partition. Available
  // answers stay in the existing locus cache and are reduced in marker order.
  void probe_required_raw_backend(const SearchInputs& input,
      const std::vector<RelationshipContexts>& side_contexts,
      const std::vector<ConfigurationScore>& configurations);

  bool message_trial_aborted() const;

  // A changed floating-point kernel must never change a threshold decision.
  // Resolve overlapping numerical intervals using the retained log kernel,
  // over the complete panel and exact proposal boundary. This is a numeric
  // refinement, not an assignment tolerance or a second search lifecycle.
  void refine_relationship_panels(
      const SearchInputs& input, int child,
      const std::vector<RelationshipContexts>& side_contexts,
      std::vector<ConfigurationScore>& configurations);


  List provenance() const;
};


void reduce_configuration_locus(
    const ConfigurationResult& result,
    ConfigurationScore& configuration);


DataFrame raw_joint_output(const EngineState& compiled);


List pair_rescue_output(const EngineState& compiled);





}
