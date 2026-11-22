#pragma once
#include <Rcpp.h>
namespace pap_internal::search {
Rcpp::List run_visits(SEXP compiled_sexp, Rcpp::List inputs);
}
