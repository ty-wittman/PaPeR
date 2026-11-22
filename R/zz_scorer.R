# R-level pedigree likelihood configuration boundary, collated after generated Rcpp wrappers.
# Keep generated RcppExports.R intact and preserve the native interface/ABI.
run_parentage_cpp <- function(compiled, compact, inputs) {
  .Call(`_PaPeR_run_parentage_cpp`, compiled, compact,
        .engine_defaults(inputs))
}
