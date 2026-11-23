# Adaptive Stage-1 candidate shortlists and incremental configuration scoring.

.section4_calibration_version <- "section4-candidate-calibration-v1"

.shortlist_settings <- function(mode, settings = NULL, required_margin_ln = 0) {
  mode <- match.arg(mode, c("adaptive", "sequoia50"))
  supplied <- !is.null(settings)
  if (is.null(settings)) settings <- list()
  if (!is.list(settings))
    stop("shortlist_settings must be NULL or a named list", call. = FALSE)
  defaults <- list(initial_size = 50L, batch_size = 50L,
                   boundary_width = 0L, expansion_margin_ln = -Inf,
                   safety_cap = Inf, required_margin_ln = required_margin_ln,
                   exhaustive = !supplied,
                   calibration_version = .section4_calibration_version)
  unknown <- setdiff(names(settings), names(defaults))
  if (length(unknown))
    stop("unknown shortlist setting(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  x <- utils::modifyList(defaults, settings)
  whole <- function(value, name, positive = FALSE) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value != floor(value) ||
        value < if (positive) 1 else 0)
      stop(name, " must be one ", if (positive) "positive" else "nonnegative",
           " whole number", call. = FALSE)
    as.integer(value)
  }
  x$initial_size <- whole(x$initial_size, "initial_size")
  x$batch_size <- whole(x$batch_size, "batch_size", TRUE)
  x$boundary_width <- whole(x$boundary_width, "boundary_width")
  if (!is.numeric(x$expansion_margin_ln) ||
      length(x$expansion_margin_ln) != 1L || is.na(x$expansion_margin_ln) ||
      is.nan(x$expansion_margin_ln))
    stop("expansion_margin_ln must be one natural-log margin", call. = FALSE)
  if (!is.numeric(x$required_margin_ln) ||
      length(x$required_margin_ln) != 1L || is.na(x$required_margin_ln) ||
      is.nan(x$required_margin_ln))
    stop("required_margin_ln must be one natural-log margin", call. = FALSE)
  if (!is.numeric(x$safety_cap) || length(x$safety_cap) != 1L ||
      is.na(x$safety_cap) || is.nan(x$safety_cap) || x$safety_cap < 1 ||
      (is.finite(x$safety_cap) && x$safety_cap != floor(x$safety_cap)))
    stop("safety_cap must be a positive whole number or Inf", call. = FALSE)
  if (!is.character(x$calibration_version) ||
      length(x$calibration_version) != 1L || is.na(x$calibration_version) ||
      !nzchar(x$calibration_version))
    stop("calibration_version must be one nonempty character value",
         call. = FALSE)
  x$exhaustive <- isTRUE(x$exhaustive)
  if (mode == "sequoia50") {
    x$initial_size <- 50L
    x$safety_cap <- 50L
    x$exhaustive <- FALSE
    x$calibration_version <- "section4-strict-sequoia50-v1"
  }
  x$mode <- mode
  x$calibration_status <- if (mode == "adaptive" && !supplied)
    "calibrated-exhaustive-fallback-v1" else
      if (mode == "adaptive") "explicit-internal-settings" else
        "strict-sequoia50-comparison"
  class(x) <- c("PaPeR_shortlist_settings", "list")
  x
}

