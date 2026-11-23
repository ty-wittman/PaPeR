#pragma once
#include <Rcpp.h>
#include <vector>
#include <string>
namespace pap_internal { struct DirtyRegion; }
namespace pap_internal::search {
struct EngineState;
using namespace Rcpp;
using pap_internal::DirtyRegion;

struct MessagePublication {
  std::vector<double> node_likelihood_ln;
  std::vector<double> locus_likelihood_ln;
  double total_likelihood_ln = 0.0;
  std::string evidence_signature;
  std::size_t signature_count = 0;
  double signature_total = 0.0;
  double signature_squared = 0.0;
  double signature_indexed = 0.0;
  bool selective = false;
  std::size_t dirty_nodes = 0;
  std::size_t total_nodes = 0;
  std::size_t nodes_visited = 0;
  std::size_t nodes_recomputed = 0;
  std::size_t message_vectors_rebuilt = 0;
  std::size_t genotype_state_elements = 0;
  std::size_t transmitted_state_elements = 0;
  std::size_t likelihood_cells = 0;
  std::size_t aggregate_updates = 0;
  std::size_t signature_global_scans = 0;
  std::size_t copies = 0;
  std::size_t allocations = 0;
  int requested_workers = 1;
  int effective_publication_workers = 1;
  int worker_thread_count = 0;
  int publication_job_count = 0;
  int worker_r_api_calls = 0;
  std::vector<double> worker_seconds;
  std::vector<int> locus_owner;
  double finalization_seconds = 0.0;
  double dispatch_seconds = 0.0;
  double seconds = 0.0;
};


MessagePublication publish_messages_full(
    EngineState& compiled, const IntegerMatrix& parents,
    int graph_revision);


void publish_messages_selective(
    EngineState& compiled, const IntegerMatrix& parents,
    int graph_revision, MessagePublication& publication,
    const DirtyRegion& dirty_region);


IntegerVector focal_order(
    const IntegerMatrix& parents, const NumericVector& expected_year,
    const IntegerVector& id_order_rank);


std::string parent_signature(
    const CharacterVector& ids, const IntegerMatrix& parents);


struct CandidateSide {
  std::vector<int> offsets;
  std::vector<int> ids;
  std::vector<double> quick;
  std::vector<int> permanent;

  CandidateSide(const List& inputs, const char* prefix, int n);
};


List refresh_active_candidates(
    List& inputs, const IntegerMatrix& parents,
    const CandidateSide& dam_source,
    const CandidateSide& sire_source);


List message_matrices(const EngineState& compiled);


NumericMatrix locus_likelihood_matrix(
    const MessagePublication& publication, int loci, int n);


List publication_work_list(const MessagePublication& publication);


}
