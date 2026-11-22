#pragma once
#include "score_records.h"
namespace pap_internal::search {
  struct DiagnosticVisitRow {
    int index=0,focal=0,graph_before=0,graph_after=0;
    bool accepted=false;
    double wall=0,preparation=0,generic_work=0,factor_preparation=0,absorption_eligible=0,absorption_distinct=0,absorption_reused=0;
  };
  struct DiagnosticPanelRow {
    int visit=0,focal=0,graph=0,dam=0,sire=0;
    double objective=R_NegInf,dam_support=R_NegInf,sire_support=R_NegInf,pair_support=R_NegInf,
      full_relationship=R_NegInf,identity=R_NegInf,joint=NA_REAL;
    bool current=false,selected=false,eligible=false,competitor=false;
    bool raw_available=false,full_available=false;std::string strongest_rival;
  };

Rcpp::DataFrame accepted_move_rows(const std::vector<AcceptedMove>& rows);
Rcpp::DataFrame focal_summary_rows(const std::vector<VisitSummary>& rows);
Rcpp::DataFrame configuration_rows(const std::vector<ConfigurationDiagnostic>& rows, const SearchInputs& input);
Rcpp::DataFrame locus_diagnostic_rows(const std::vector<ConfigurationLocusDiagnostic>& rows);
DataFrame visit_diagnostic_rows(const std::vector<DiagnosticVisitRow>& rows);
DataFrame panel_diagnostic_rows(const std::vector<DiagnosticPanelRow>& rows);
DataFrame relationship_provenance_rows(const std::vector<ConfigurationDiagnostic>& rows);
List graph_diagnostic_rows(const std::vector<GraphStateDiagnostic>& rows, int n);
}
