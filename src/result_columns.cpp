#include "result_columns.h"
#include <functional>
namespace pap_internal::search {
namespace {
template<class Vector,class Rows,class Accessor>
Vector row_column(const Rows& rows,Accessor access) {
 Vector values(rows.size());
 for(std::size_t i=0;i<rows.size();++i)values[i]=access(rows[i]);
 return values;
}

// A column's R type, name and row accessor are declared together, in output order.
template<class Row>
class RowColumns {
 const std::vector<Row>& rows;
 Rcpp::List columns;
 public:
 explicit RowColumns(const std::vector<Row>& source):rows(source) {}
 template<class Vector,class Accessor>
 void add(const std::string& name,Accessor access) {
   Vector values=row_column<Vector>(rows,[&](const Row& row){return std::invoke(access,row);});
   columns.push_back(values,name);
 }
 Rcpp::DataFrame frame()const{return Rcpp::DataFrame(columns);}
};
}

DataFrame accepted_move_rows(const std::vector<AcceptedMove>& rows) {
  RowColumns<AcceptedMove> columns(rows);
  columns.add<IntegerVector>("child", &AcceptedMove::child);
  columns.add<IntegerVector>("old_dam", &AcceptedMove::old_dam);
  columns.add<IntegerVector>("old_sire", &AcceptedMove::old_sire);
  columns.add<IntegerVector>("new_dam", &AcceptedMove::new_dam);
  columns.add<IntegerVector>("new_sire", &AcceptedMove::new_sire);
  columns.add<NumericVector>("current_objective_ln", &AcceptedMove::current_objective_ln);
  columns.add<NumericVector>("selected_objective_ln", &AcceptedMove::selected_objective_ln);
  columns.add<NumericVector>("objective_delta_ln", &AcceptedMove::objective_delta_ln);
  columns.add<NumericVector>("dam_support_ln", &AcceptedMove::dam_support_ln);
  columns.add<NumericVector>("sire_support_ln", &AcceptedMove::sire_support_ln);
  columns.add<NumericVector>("pair_support_ln", &AcceptedMove::pair_support_ln);
  columns.add<NumericVector>("identity_margin_ln", &AcceptedMove::identity_margin_ln);
  columns.add<IntegerVector>("trio_me_count", &AcceptedMove::trio_me);
  columns.add<IntegerVector>("graph_revision_before", &AcceptedMove::graph_revision_before);
  columns.add<IntegerVector>("graph_revision_after", &AcceptedMove::graph_revision_after);
  return columns.frame();
}

DataFrame focal_summary_rows(const std::vector<VisitSummary>& rows) {
  RowColumns<VisitSummary> columns(rows);
  columns.add<IntegerVector>("focal", &VisitSummary::focal);
  columns.add<IntegerVector>("starting_dam", &VisitSummary::starting_dam);
  columns.add<IntegerVector>("starting_sire", &VisitSummary::starting_sire);
  columns.add<IntegerVector>("selected_dam", &VisitSummary::selected_dam);
  columns.add<IntegerVector>("selected_sire", &VisitSummary::selected_sire);
  columns.add<NumericVector>("current_objective_ln", &VisitSummary::current_objective_ln);
  columns.add<NumericVector>("selected_objective_ln", &VisitSummary::selected_objective_ln);
  columns.add<NumericVector>("dam_support_ln", &VisitSummary::dam_support_ln);
  columns.add<NumericVector>("sire_support_ln", &VisitSummary::sire_support_ln);
  columns.add<NumericVector>("pair_support_ln", &VisitSummary::pair_support_ln);
  columns.add<NumericVector>("identity_margin_ln", &VisitSummary::identity_margin_ln);
  columns.add<NumericVector>("full_relationship_ln", &VisitSummary::full_relationship_ln);
  columns.add<IntegerVector>("trio_me_count", &VisitSummary::trio_me);
  columns.add<IntegerVector>("trio_me_pass", &VisitSummary::trio_me_pass);
  columns.add<IntegerVector>("near_tie", &VisitSummary::near_tie);
  columns.add<IntegerVector>("accepted", &VisitSummary::accepted);
  columns.add<IntegerVector>("message_scoring_model", &VisitSummary::message_scoring_model);
  columns.add<IntegerVector>("raw_message_approximate", &VisitSummary::raw_message_approximate);
  columns.add<IntegerVector>("rejection_reason", &VisitSummary::rejection_reason);
  columns.add<IntegerVector>("graph_revision_before", &VisitSummary::graph_revision_before);
  columns.add<IntegerVector>("graph_revision_after", &VisitSummary::graph_revision_after);
  columns.add<IntegerVector>("dam_candidate_count", &VisitSummary::dam_candidate_count);
  columns.add<IntegerVector>("sire_candidate_count", &VisitSummary::sire_candidate_count);
  columns.add<IntegerVector>("configurations_total", &VisitSummary::configurations_total);
  columns.add<IntegerVector>("single_configurations_total", &VisitSummary::single_configurations_total);
  columns.add<IntegerVector>("pair_configurations_total", &VisitSummary::pair_configurations_total);
  columns.add<IntegerVector>("single_hard_feasible", &VisitSummary::single_hard_feasible);
  columns.add<IntegerVector>("pair_hard_feasible", &VisitSummary::pair_hard_feasible);
  columns.add<IntegerVector>("single_trio_me_pass", &VisitSummary::single_trio_me_pass);
  columns.add<IntegerVector>("pair_trio_me_pass", &VisitSummary::pair_trio_me_pass);
  columns.add<IntegerVector>("single_relationship_pass", &VisitSummary::single_relationship_pass);
  columns.add<IntegerVector>("pair_relationship_pass", &VisitSummary::pair_relationship_pass);
  columns.add<IntegerVector>("single_eligible", &VisitSummary::single_eligible);
  columns.add<IntegerVector>("pair_eligible", &VisitSummary::pair_eligible);
  columns.add<IntegerVector>("configurations_eligible", &VisitSummary::configurations_eligible);
  columns.add<IntegerVector>("current_hard_valid", &VisitSummary::current_hard_valid);
  columns.add<IntegerVector>("current_assignment_confident", &VisitSummary::current_assignment_confident);
  columns.add<IntegerVector>("current_retention_eligible", &VisitSummary::current_retention_eligible);
  columns.add<IntegerVector>("recheck_evaluation", &VisitSummary::recheck_evaluation);
  columns.add<IntegerVector>("complete_pair_skipped", &VisitSummary::complete_pair_skipped);
  columns.add<NumericVector>("focal_wall_seconds", &VisitSummary::focal_wall_seconds);
  columns.add<NumericVector>("chooser_seconds", &VisitSummary::chooser_seconds);
  columns.add<NumericVector>("recursive_family_to_cavity_calls", [](const auto& row){return 0.0;});
  return columns.frame();
}

DataFrame configuration_rows(const std::vector<ConfigurationDiagnostic>& rows, const SearchInputs& input) {
  RowColumns<ConfigurationDiagnostic> columns(rows);
  columns.add<IntegerVector>("focal_position", &ConfigurationDiagnostic::focal_position);
  columns.add<IntegerVector>("focal", &ConfigurationDiagnostic::focal);
  columns.add<IntegerVector>("graph_revision", &ConfigurationDiagnostic::graph_revision);
  columns.add<IntegerVector>("configuration_position", &ConfigurationDiagnostic::configuration_position);
  columns.add<IntegerVector>("dam", &ConfigurationDiagnostic::dam);
  columns.add<IntegerVector>("sire", &ConfigurationDiagnostic::sire);
  columns.add<IntegerVector>("nparents", &ConfigurationDiagnostic::nparents);
  columns.add<IntegerVector>("status", &ConfigurationDiagnostic::status);
  columns.add<IntegerVector>("reason", &ConfigurationDiagnostic::reason);
  columns.add<IntegerVector>("route", &ConfigurationDiagnostic::route);
  columns.add<IntegerVector>("trio_me_count", &ConfigurationDiagnostic::trio_me);
  columns.add<NumericVector>("objective_ln", &ConfigurationDiagnostic::objective_ln);
  columns.add<IntegerVector>("raw_comparison_available", &ConfigurationDiagnostic::raw_comparison_available);
  columns.add<CharacterVector>("raw_status", [](const auto& row){return row.raw_status == ::relationship_topology::configuration::PredictiveStatus::DefinedFinite ? "defined_finite" : row.raw_status == ::relationship_topology::configuration::PredictiveStatus::DefinedGeneticZero ? "defined_genetic_zero" : row.raw_status == ::relationship_topology::configuration::PredictiveStatus::ConditioningZero ? "conditioning_zero" : "unavailable";});
  columns.add<NumericVector>("dam_support_ln", &ConfigurationDiagnostic::dam_support_ln);
  columns.add<NumericVector>("sire_support_ln", &ConfigurationDiagnostic::sire_support_ln);
  columns.add<NumericVector>("pair_support_ln", &ConfigurationDiagnostic::pair_support_ln);
  columns.add<NumericVector>("identity_margin_ln", &ConfigurationDiagnostic::identity_margin_ln);
  columns.add<IntegerVector>("hard_feasible", &ConfigurationDiagnostic::hard_feasible);
  columns.add<IntegerVector>("hard_valid", &ConfigurationDiagnostic::hard_valid);
  columns.add<IntegerVector>("trio_me_pass", &ConfigurationDiagnostic::trio_me_pass);
  columns.add<IntegerVector>("relationship_pass", &ConfigurationDiagnostic::relationship_pass);
  columns.add<IntegerVector>("relationship_comparison_available", &ConfigurationDiagnostic::relationship_comparison_available);
  columns.add<IntegerVector>("full_relationship_gate_pass", &ConfigurationDiagnostic::full_relationship_gate_pass);
  columns.add<NumericVector>("full_relationship_ln", &ConfigurationDiagnostic::full_relationship_ln);
  columns.add<NumericVector>("full_operational_ln", &ConfigurationDiagnostic::full_operational_ln);
  columns.add<NumericVector>("identity_objective_ln", &ConfigurationDiagnostic::identity_objective_ln);
  columns.add<NumericVector>("full_alternative_ln", &ConfigurationDiagnostic::full_alternative_ln);
  columns.add<IntegerVector>("assignment_confident", &ConfigurationDiagnostic::assignment_confident);
  columns.add<IntegerVector>("eligible_for_retention", &ConfigurationDiagnostic::eligible_for_retention);
  columns.add<IntegerVector>("eligible_for_addition", &ConfigurationDiagnostic::eligible_for_addition);
  columns.add<IntegerVector>("eligible_for_replacement", &ConfigurationDiagnostic::eligible_for_replacement);
  columns.add<IntegerVector>("eligible_for_removal", &ConfigurationDiagnostic::eligible_for_removal);
  columns.add<IntegerVector>("eligible", &ConfigurationDiagnostic::eligible);
  columns.add<IntegerVector>("near_tie", &ConfigurationDiagnostic::near_tie);
  columns.add<IntegerVector>("is_current", &ConfigurationDiagnostic::is_current);
  columns.add<IntegerVector>("is_selected", &ConfigurationDiagnostic::is_selected);
  columns.add<IntegerVector>("current_dam_before_test", &ConfigurationDiagnostic::current_dam_before_test);
  columns.add<IntegerVector>("current_sire_before_test", &ConfigurationDiagnostic::current_sire_before_test);
  columns.add<IntegerVector>("focal_rejection_reason", &ConfigurationDiagnostic::focal_rejection_reason);
  DataFrame configuration_frame = columns.frame();

    configuration_frame["raw_reference_error"] = row_column<NumericVector>(rows, [](const auto& row){return row.raw_reference_error;});
    configuration_frame["side_reference_uncertainty_dam"] = row_column<NumericVector>(rows, [](const auto& row){return row.side_reference_uncertainty[0];});
    configuration_frame["side_reference_uncertainty_sire"] = row_column<NumericVector>(rows, [](const auto& row){return row.side_reference_uncertainty[1];});
    static constexpr std::array<const char*, kRelationshipCount>
      relationship_id{{"U", "FS", "HS", "GP", "FA", "HA"}};
    for (int side = 0; side < 2; ++side)
      for (int relationship = 0;
           relationship < kRelationshipCount; ++relationship)
        configuration_frame[std::string("relationship_") +
          (side == 0 ? "dam_" : "sire_") +
          relationship_id[static_cast<std::size_t>(relationship)] + "_ln"] =
            row_column<NumericVector>(rows, [=](const auto& row){return row.full_relationship_member_ln[static_cast<std::size_t>(
              side * kRelationshipCount + relationship)];});
    static constexpr std::array<const char*, kFullParentNumeratorCount>
      parent_id{{"PO", "HS_PO", "PO_HA"}};
    for (int member = 0; member < kFullParentNumeratorCount; ++member)
      configuration_frame[std::string("relationship_parent_") +
        parent_id[static_cast<std::size_t>(member)] + "_ln"] =
          row_column<NumericVector>(rows, [=](const auto& row){return row.full_parent_numerator_ln[static_cast<std::size_t>(member)];});
    if (input.relationship_mode == 1) {
      configuration_frame["full_common_boundary"] = row_column<IntegerVector>(rows, [](const auto& row){return row.full_common_boundary;});
      configuration_frame["full_limiting_side"] = row_column<IntegerVector>(rows, [](const auto& row){return row.full_limiting_side;});
      configuration_frame["full_dam_relationship_ln"] = row_column<NumericVector>(rows, [](const auto& row){return row.full_side_relationship_ln[0];});
      configuration_frame["full_sire_relationship_ln"] = row_column<NumericVector>(rows, [](const auto& row){return row.full_side_relationship_ln[1];});
      configuration_frame["full_parent_model"] = row_column<CharacterVector>(rows, [](const auto& row){return row.full_parent_model;});
      configuration_frame["full_alternative_model"] = row_column<CharacterVector>(rows, [](const auto& row){return row.full_alternative_model;});
      for (int side = 0; side < 2; ++side)
        for (int relationship = 0;
             relationship < kRelationshipCount; ++relationship)
          configuration_frame[std::string("full_") +
            (side == 0 ? "dam_" : "sire_") +
            relationship_id[static_cast<std::size_t>(relationship)] + "_ln"] =
              row_column<NumericVector>(rows, [=](const auto& row){return row.full_relationship_member_ln[static_cast<std::size_t>(
                side * kRelationshipCount + relationship)];});
      for (int member = 0; member < kFullParentNumeratorCount; ++member)
        configuration_frame[std::string("full_parent_") +
          parent_id[static_cast<std::size_t>(member)] + "_ln"] =
            row_column<NumericVector>(rows, [=](const auto& row){return row.full_parent_numerator_ln[static_cast<std::size_t>(member)];});
      static constexpr std::array<const char*, kFullContextualRelationshipCount>
        contextual_id{{"HS_HA", "HS_F1C", "HS_FA", "HS_GP", "HS_GGG", "FA_HA",
                       "FS_HA", "FS_H1C", "FS_PARENTS_FS"}};
      for (int side = 0; side < 2; ++side)
        for (int relationship = 0;
             relationship < kFullContextualRelationshipCount; ++relationship) {
          const int registry_index = input.full_contextual_relationship_index[
            static_cast<std::size_t>(relationship)];
          const unsigned char allowed = input.full_relationship_orientation_mask[
            static_cast<std::size_t>(registry_index)];
          for (int orientation = 0; orientation < kFullContextualOrientationCount;
               ++orientation) {
            const unsigned char bit = orientation == 0 ?
              FULL_ORIENTATION_FORWARD : FULL_ORIENTATION_REVERSE;
            // HS_HA's directed reverse is distinct on contextual evidence
            // and is included even when the registry mask is symmetric.
            const bool contextual_reverse = relationship == 0 && orientation == 1;
            if (!(allowed & bit) && !contextual_reverse) continue;
            const std::string direction = orientation == 1 ? "reverse" :
              ((allowed & FULL_ORIENTATION_REVERSE) ? "forward" : "symmetric");
            configuration_frame[std::string("full_") +
              (side == 0 ? "dam_" : "sire_") +
              contextual_id[static_cast<std::size_t>(relationship)] + "_" +
              direction + "_ln"] = row_column<NumericVector>(rows, [=](const auto& row){return row.full_contextual_member_ln[
                contextual_member_offset(side, relationship, orientation)];});
          }
        }
    }

  return configuration_frame;
}

DataFrame locus_diagnostic_rows(const std::vector<ConfigurationLocusDiagnostic>& rows) {
  RowColumns<ConfigurationLocusDiagnostic> columns(rows);
  columns.add<IntegerVector>("focal_position", &ConfigurationLocusDiagnostic::focal_position);
  columns.add<IntegerVector>("focal", &ConfigurationLocusDiagnostic::focal);
  columns.add<IntegerVector>("graph_revision", &ConfigurationLocusDiagnostic::graph_revision);
  columns.add<IntegerVector>("configuration_position", &ConfigurationLocusDiagnostic::configuration_position);
  columns.add<IntegerVector>("locus", &ConfigurationLocusDiagnostic::locus);
  columns.add<NumericVector>("objective_ln", &ConfigurationLocusDiagnostic::objective_ln);
  return columns.frame();
}


DataFrame visit_diagnostic_rows(const std::vector<DiagnosticVisitRow>& rows) {
  RowColumns<DiagnosticVisitRow> columns(rows);
  columns.add<IntegerVector>("visit_index", &DiagnosticVisitRow::index);
  columns.add<IntegerVector>("focal", &DiagnosticVisitRow::focal);
  columns.add<IntegerVector>("graph_revision_before", &DiagnosticVisitRow::graph_before);
  columns.add<IntegerVector>("graph_revision_after", &DiagnosticVisitRow::graph_after);
  columns.add<LogicalVector>("accepted", &DiagnosticVisitRow::accepted);
  columns.add<NumericVector>("wall_seconds", &DiagnosticVisitRow::wall);
  columns.add<NumericVector>("preparation_seconds", &DiagnosticVisitRow::preparation);
  columns.add<NumericVector>("generic_work_cells", &DiagnosticVisitRow::generic_work);
  columns.add<NumericVector>("factor_preparation_seconds", &DiagnosticVisitRow::factor_preparation);
  columns.add<NumericVector>("absorption_eligible", &DiagnosticVisitRow::absorption_eligible);
  columns.add<NumericVector>("absorption_distinct", &DiagnosticVisitRow::absorption_distinct);
  columns.add<NumericVector>("absorption_reused", &DiagnosticVisitRow::absorption_reused);
  return columns.frame();
}

DataFrame panel_diagnostic_rows(const std::vector<DiagnosticPanelRow>& rows) {
  RowColumns<DiagnosticPanelRow> columns(rows);
  columns.add<IntegerVector>("visit_index", &DiagnosticPanelRow::visit);
  columns.add<IntegerVector>("focal", &DiagnosticPanelRow::focal);
  columns.add<IntegerVector>("graph_revision", &DiagnosticPanelRow::graph);
  columns.add<IntegerVector>("dam", &DiagnosticPanelRow::dam);
  columns.add<IntegerVector>("sire", &DiagnosticPanelRow::sire);
  columns.add<NumericVector>("objective_ln", &DiagnosticPanelRow::objective);
  columns.add<NumericVector>("dam_support_ln", &DiagnosticPanelRow::dam_support);
  columns.add<NumericVector>("sire_support_ln", &DiagnosticPanelRow::sire_support);
  columns.add<NumericVector>("pair_support_ln", &DiagnosticPanelRow::pair_support);
  columns.add<NumericVector>("full_relationship_ln", &DiagnosticPanelRow::full_relationship);
  columns.add<NumericVector>("identity_margin_ln", &DiagnosticPanelRow::identity);
  columns.add<NumericVector>("identity_objective_ln", &DiagnosticPanelRow::joint);
  columns.add<LogicalVector>("identity_competitor", &DiagnosticPanelRow::competitor);
  columns.add<LogicalVector>("is_current", &DiagnosticPanelRow::current);
  columns.add<LogicalVector>("is_selected", &DiagnosticPanelRow::selected);
  columns.add<LogicalVector>("eligible", &DiagnosticPanelRow::eligible);
  columns.add<LogicalVector>("raw_available", &DiagnosticPanelRow::raw_available);
  columns.add<LogicalVector>("full_available", &DiagnosticPanelRow::full_available);
  columns.add<CharacterVector>("strongest_rival", &DiagnosticPanelRow::strongest_rival);
  return columns.frame();
}

DataFrame relationship_provenance_rows(const std::vector<ConfigurationDiagnostic>& configurations) {
 struct ProvenanceRow {
  const ConfigurationDiagnostic* configuration;
  const RelationshipProvenanceDiagnostic* source;
 };
 std::vector<ProvenanceRow> rows;
 for(const auto& row:configurations)for(int side=0;side<2;++side)
  if(row.relationship_provenance_seen[side])rows.push_back({&row,&row.relationship_provenance[side]});
 RowColumns<ProvenanceRow> columns(rows);
  columns.add<IntegerVector>("focal_position", [](const auto& row){return row.configuration->focal_position;});
  columns.add<IntegerVector>("configuration_position", [](const auto& row){return row.configuration->configuration_position;});
  columns.add<IntegerVector>("graph_revision", [](const auto& row){return row.configuration->graph_revision;});
  columns.add<IntegerVector>("tested_parent_side", [](const auto& row){return row.source->tested_parent_side;});
  columns.add<IntegerVector>("focal", [](const auto& row){return row.source->focal;});
  columns.add<IntegerVector>("candidate", [](const auto& row){return row.source->candidate;});
  columns.add<IntegerVector>("proposed_dam", [](const auto& row){return row.source->proposed_dam;});
  columns.add<IntegerVector>("proposed_sire", [](const auto& row){return row.source->proposed_sire;});
  columns.add<IntegerVector>("masked_focal_dam", [](const auto& row){return row.source->focal_dam;});
  columns.add<IntegerVector>("masked_focal_sire", [](const auto& row){return row.source->focal_sire;});
  columns.add<IntegerVector>("masked_candidate_dam", [](const auto& row){return row.source->candidate_dam;});
  columns.add<IntegerVector>("masked_candidate_sire", [](const auto& row){return row.source->candidate_sire;});
  columns.add<IntegerVector>("fixed_other_parent", [](const auto& row){return row.source->fixed_other_parent;});
  columns.add<IntegerVector>("has_relevant_context", [](const auto& row){return row.source->has_relevant_context;});
  columns.add<IntegerVector>("focal_endpoint_excluded_child", [](const auto& row){return row.source->focal_endpoint_excluded_child;});
  columns.add<IntegerVector>("candidate_endpoint_excluded_child", [](const auto& row){return row.source->candidate_endpoint_excluded_child;});
  columns.add<IntegerVector>("fs_orientation_mask", [](const auto& row){return row.source->fs_orientation_mask;});
  columns.add<IntegerVector>("hs_orientation_mask", [](const auto& row){return row.source->hs_orientation_mask;});
  columns.add<IntegerVector>("gp_orientation_mask", [](const auto& row){return row.source->gp_orientation_mask;});
  columns.add<IntegerVector>("fa_orientation_mask", [](const auto& row){return row.source->fa_orientation_mask;});
  columns.add<IntegerVector>("ha_orientation_mask", [](const auto& row){return row.source->ha_orientation_mask;});
  columns.add<IntegerVector>("loci_accumulated", [](const auto& row){return row.source->loci_accumulated;});
  columns.add<NumericVector>("context_signature", [](const auto& row){return static_cast<double>(row.source->context_signature);});
  columns.add<NumericVector>("proposed_PO_partition_ln", [](const auto& row){return row.source->proposed_po_partition_ln;});
 DataFrame provenance_frame = columns.frame();

    static constexpr std::array<const char*, kRelationshipCount>
      provenance_relationship_id{{"U", "FS", "HS", "GP", "FA", "HA"}};
    for (int relationship = 0;
         relationship < kRelationshipCount; ++relationship)
      provenance_frame[std::string("partition_") +
        provenance_relationship_id[static_cast<std::size_t>(relationship)] +
        "_ln"] = row_column<NumericVector>(rows, [=](const auto& row){return row.source->partition_ln[static_cast<std::size_t>(relationship)];});
    static constexpr std::array<const char*, 22> orientation_id{{
      "HS_dam", "HS_sire",
      "GP_focal_dam_dam", "GP_focal_dam_sire",
      "GP_focal_sire_dam", "GP_focal_sire_sire",
      "GP_candidate_dam_dam", "GP_candidate_dam_sire",
      "GP_candidate_sire_dam", "GP_candidate_sire_sire",
      "FA_focal_dam", "FA_focal_sire", "FA_candidate_dam", "FA_candidate_sire",
      "HA_focal_dam_dam", "HA_focal_dam_sire", "HA_focal_sire_dam",
      "HA_focal_sire_sire", "HA_candidate_dam_dam", "HA_candidate_dam_sire",
      "HA_candidate_sire_dam", "HA_candidate_sire_sire"}};
    for (int orientation = 0; orientation < 22; ++orientation)
      provenance_frame[std::string("orientation_") +
        orientation_id[static_cast<std::size_t>(orientation)] + "_ln"] =
          row_column<NumericVector>(rows, [=](const auto& row){return row.source->orientation_partition_ln[static_cast<std::size_t>(orientation)];});

 return provenance_frame;
}

List graph_diagnostic_rows(const std::vector<GraphStateDiagnostic>& rows, int n) {
    const std::size_t nstates = rows.size();
    IntegerVector state_revision(nstates);
    IntegerMatrix state_dam(nstates, n), state_sire(nstates, n);
    for (std::size_t state_index = 0; state_index < nstates; ++state_index) {
      const GraphStateDiagnostic& snapshot =
        rows[state_index];
      state_revision[state_index] = snapshot.graph_revision;
      for (int individual = 0; individual < n; ++individual) {
        state_dam(static_cast<int>(state_index), individual) =
          snapshot.dam[static_cast<std::size_t>(individual + 1)];
        state_sire(static_cast<int>(state_index), individual) =
          snapshot.sire[static_cast<std::size_t>(individual + 1)];
      }
    }
    return List::create(
      _["graph_revision"] = state_revision,
      _["dam"] = state_dam, _["sire"] = state_sire);
}

}
