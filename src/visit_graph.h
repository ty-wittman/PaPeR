#pragma once
#include "engine_state.h"
#include "execution_config.h"
namespace pap_internal::search {

struct VisitGraph {
  pap_internal::EvaluationOptions options;
  EngineState& compiled;
  // Each worker owns its availability state and shares immutable locus data.
  EngineCounters& telemetry;
  int focal = 0;
  int graph_revision = 0;

  // Use compiled local-pedigree programs for configuration scoring.
  bool bounded_messages = false;

  bool raw_scaled_fast = false;
  bool message_scoring = false;

  std::vector<int> dam;
  std::vector<int> sire;
  std::vector<std::vector<int>> reverse;
  // Graph identity is independent of observation storage. Stage 1 initializes
  // the current one-to-one typed mapping once per pass; reusable graph views
  // may instead carry kNoEvidenceRow for latent nodes.
  std::vector<int> node_to_evidence_row;
  std::set<int> evidence_exclusions;
  std::map<std::pair<int, int>, IbdClass> ibd_cache;
  VisitGraph(EngineState& compiled_, int focal_,
           int graph_revision_, const IntegerMatrix& parents,
           EngineCounters* telemetry_override = nullptr);

  VisitGraph(EngineState& compiled_, int focal_,
           int graph_revision_, const std::vector<int>& parents_column_major,
           EngineCounters& telemetry_override);

  int parent(int child, int side) const;

  void set_focal(int child, int revision);

  void remove_reverse_edge(int parent_id, int child);

  void add_reverse_edge(int parent_id, int child);

  void apply_parent_move(int child, int new_dam, int new_sire);

  std::vector<int> ancestors(int node) const;

  IbdClass ibd(int first, int second);

 private:
  void initialize_parents(const int* parents, bool scheduled_worker);
};
}
