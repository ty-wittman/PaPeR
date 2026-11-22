# Progress display for assign_parents(control = list(verbose = TRUE)).
#
# On a terminal or in an interactive session one status line (search, round,
# a short bar and the percent of the round done) is redrawn in place every
# two seconds; elsewhere (log files, knitr) one line is written per round.
# options(PaPeR.progress_dynamic = TRUE/FALSE) forces either mode. Everything
# goes through message(), so suppressMessages() silences it. A NULL reporter
# (verbose = FALSE) makes every call a no-op.

.pap_progress <- function(verbose) {
  if (!isTRUE(verbose)) return(NULL)
  p <- new.env(parent = emptyenv())
  p$dynamic <- isTRUE(getOption("PaPeR.progress_dynamic", interactive() || isatty(stderr())))
  p$started <- proc.time()[["elapsed"]]
  p$open <- 0L          # width of the status line currently on screen
  p$last_draw <- -Inf
  p$search <- 0L
  p$pass <- ""
  p
}

.pap_clock <- function(p) {
  s <- round(proc.time()[["elapsed"]] - p$started)
  sprintf("%d:%02d", s %/% 60, s %% 60)
}

# A permanent line (closes any status line first).
.pap_say <- function(p, ...) {
  if (is.null(p)) return(invisible(NULL))
  if (p$open > 0L) { message("\r", strrep(" ", p$open), "\r", appendLF = FALSE); p$open <- 0L }
  message(...)
  invisible(NULL)
}

# A transient status line, redrawn in place; ignored when not on a terminal.
.pap_status <- function(p, text) {
  if (is.null(p) || !p$dynamic) return(invisible(NULL))
  text <- paste0(text, "  ", .pap_clock(p))
  # A line wider than the console wraps, and "\r" cannot redraw it in place.
  text <- substr(text, 1L, max(20L, min(getOption("width", 80L), 80L) - 1L))
  message("\r", formatC(text, width = -max(p$open, nchar(text))), appendLF = FALSE)
  p$open <- nchar(text)
  p$last_draw <- proc.time()[["elapsed"]]
  invisible(NULL)
}

.pap_redraw_seconds <- 2

# The native search calls this at every focal visit with
# c(round, recheck, visits done, visits in the pass, moves so far in the pass).
.pap_visit_callback <- function(p) {
  function(v) {
    r <- v[[1L]]; recheck <- v[[2L]] == 1L
    done <- v[[3L]]; total <- v[[4L]]; moves <- v[[5L]]
    label <- sprintf("search %d  round %d%s", p$search, r, if (recheck) " recheck" else "")
    if (!identical(label, p$pass)) { p$pass <- label; p$pass_started <- proc.time()[["elapsed"]] }
    finished <- done >= total
    if (p$dynamic) {
      if (proc.time()[["elapsed"]] - p$last_draw >= .pap_redraw_seconds) {
        fraction <- if (total > 0L) done / total else 1
        filled <- round(10 * fraction)
        .pap_status(p, sprintf("%-22s [%s%s] %3d%%", label, strrep("#", filled),
                               strrep("-", 10L - filled), as.integer(floor(100 * fraction))))
      }
    } else if (finished && !recheck) {
      .pap_say(p, sprintf("%s: %d moves, %.1f s", label, moves,
                          proc.time()[["elapsed"]] - p$pass_started))
    }
    invisible(NULL)
  }
}

.pap_rate_text <- function(rates) {
  paste(names(rates), signif(unname(rates), 3), sep = " ", collapse = ", ")
}

# One line per finished search: rounds, assignments, the fit and the decision.
.pap_report_search <- function(p, record, seconds) {
  if (is.null(p)) return(invisible(NULL))
  ped <- record$parents
  native <- record$native_seconds %||% NA_real_
  prepare <- (record$search_seconds %||% seconds) - native
  fitting <- record$fit_seconds %||% 0
  text <- sprintf("search %d: %d rounds, %d dams, %d sires (%.1f s: %.1f preparing, %.1f searching, %.1f fitting)",
                  record$pass, as.integer(record$rounds), sum(!is.na(ped$dam)),
                  sum(!is.na(ped$sire)), seconds, prepare, native, fitting)
  fit <- record$fit
  if (is.list(fit) && !is.null(fit$estimated_rates)) {
    name <- names(record$rates_used)[1L]
    iv <- fit$interval
    bounds <- if (is.matrix(iv) && name %in% rownames(iv)) iv[name, c("lower", "upper")] else c(NA, NA)
    text <- paste0(text, sprintf("; %s fit %s [%s, %s]", name,
      signif(fit$estimated_rates[[name]], 3), signif(bounds[[1L]], 3), signif(bounds[[2L]], 3)))
  }
  decision <- record$stop_reason %||% NA_character_
  if (identical(decision, "continue")) {
    decision <- paste("next search with", .pap_rate_text(record$next_rates[1L]))
  } else if (!is.na(decision)) {
    decision <- if (decision %in% names(.pap_decision_text)) .pap_decision_text[[decision]] else
      gsub("_", " ", decision)
  }
  if (!is.na(decision)) text <- paste0(text, "; ", decision)
  .pap_say(p, text)
}

# Plain-language calibration decisions for the per-search line.
.pap_decision_text <- c(
  supplied_rate_retained = "rate kept (within its bounds)",
  guarded_lowering_completed = "rate kept after one lowering",
  insufficient_confident_fraction = "rate kept (too few confident trios to lower it)",
  maximum_searches = "stopped after four searches",
  direction_reversal = "stopped (direction reversal)",
  too_few_trios = "rate kept (too few confident trios to fit)",
  unidentifiable_zero_ceiling = "rate kept (fit not identifiable)",
  fit_not_converged = "rate kept (fit did not converge)",
  pass_not_converged = "rate kept (search did not converge)",
  missing_ceiling_or_observations = "rate kept (no Mendelian-error observations)",
  not_calibrated = "not calibrated")
