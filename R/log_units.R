# Explicit public/internal log-unit boundaries.

.llr_log10 <- function(llr_ln) {
  if (!is.numeric(llr_ln))
    stop("llr_ln must be numeric", call. = FALSE)
  llr_ln / log(10)
}

.threshold_log10_to_ln <- function(threshold_log10) {
  if (!is.numeric(threshold_log10))
    stop("threshold_log10 must be numeric", call. = FALSE)
  threshold_log10 * log(10)
}

.public_log_base <- 10L
