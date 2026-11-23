#' Assign genotyped parents to genotyped individuals
#'
#' Finds the dam and sire of each genotyped individual among the other
#' genotyped individuals. Allele frequencies are estimated from \code{GenoM},
#' locus tables are built from the error rates, and the error rate is
#' calibrated from the inferred pedigree (see \code{error_calibration} in
#' \code{\link{pap_control}}). Remove duplicate samples first with
#' \code{\link{find_duplicates}}.
#'
#' @param GenoM genotype matrix or data frame. Rows are individuals (row names
#'   are IDs) and columns are loci. Calls are \code{"a/b"} allele pairs;
#'   \code{"0"} is a missing call. SNPs use alleles 1 and 2.
#' @param LifeHistData data frame with \code{ID}, \code{Sex} (1 female, 2 male,
#'   3 or NA unknown) and \code{BirthYear}, optionally \code{BY.min},
#'   \code{BY.max} and \code{Year.last}. Every genotyped individual needs a
#'   birth year or both \code{BY.min} and \code{BY.max}.
#' @param Dropout probability that a true heterozygote is called as one of its
#'   two homozygotes. One positive number for the whole panel.
#' @param Foreign per-allele-copy probability that a copy is replaced by a
#'   random allele drawn from the population. One positive number. Calibration
#'   never changes it.
#' @param control named list of further settings; see \code{\link{pap_control}}.
#'   Unknown names are an error.
#' @return A data frame with columns \code{id}, \code{dam} and \code{sire}
#'   (\code{NA} where no parent is assigned), plus \code{LLRdam},
#'   \code{LLRsire}, \code{LLRpair}, \code{OHdam}, \code{OHsire}, \code{MEpair}
#'   and availability statuses (see \code{\link{calc_llr}}) when
#'   \code{control$CalcLLR} is \code{TRUE}. Attributes:
#'   \describe{
#'     \item{\code{calibration}}{one row per search: rates used, fitted rate and
#'       95\% bounds, confident trios and the decision.}
#'     \item{\code{calibration_passes}}{full per-search records.}
#'     \item{\code{converged}, \code{rounds}}{search status of the final search.}
#'     \item{\code{settings}}{the resolved \code{control} list.}
#'   }
#' @examples
#' ex <- pap_example
#' # A quick run on the first two birth cohorts (224 individuals)
#' first <- ex$LifeHistData$ID[ex$LifeHistData$BirthYear <= 2016]
#' ped <- assign_parents(ex$GenoM[first, ], ex$LifeHistData,
#'                       Dropout = 0.005, Foreign = 5e-4)
#' head(ped[!is.na(ped$dam), ])
#' attr(ped, "calibration")
#' score_pedigree(ped, ex$Pedigree, typed_ids = first)$dam[c("correct", "wrong", "missed")]
#'
#' \donttest{
#' # The whole example population with final LLRs (about a minute)
#' ped <- assign_parents(ex$GenoM, ex$LifeHistData, Dropout = 0.005,
#'                       Foreign = 5e-4, control = list(CalcLLR = TRUE))
#' summary(ped$LLRdam)
#'
#' # SNPs with Sequoia's error model
#' ped_snp <- assign_parents(ex$GenoM_snp, ex$LifeHistData,
#'                           control = list(error_model = "sequoia", Err = 0.005))
#' }
#' @export
assign_parents <- function(GenoM, LifeHistData, Dropout = NULL, Foreign = NULL,
                           control = list()) {
  settings <- pap_control(control)
  started <- proc.time()[["elapsed"]]
  pap_removed_path_hit_cpp()
  progress <- .pap_progress(settings$verbose)
  tables <- .pap_tables(GenoM, Dropout, Foreign, settings)
  .require_positive_error_rates(tables)
  .pap_say(progress, sprintf("PaPeR: %d individuals, %d loci; %s; calibration %s",
    nrow(GenoM), ncol(GenoM), .pap_rate_text(.pap_current_rates(tables)),
    settings$error_calibration))
  if (settings$error_calibration == "none") {
    run <- .pap_search(GenoM, LifeHistData, tables, settings, progress = progress)
    passes <- list(.pap_pass_record(1L, .pap_current_rates(tables), run, NULL,
                                    list(status = "not_calibrated"), NULL))
    .pap_report_search(progress, passes[[1L]], proc.time()[["elapsed"]] - started)
    final <- run$pedigree
    final_tables <- tables
  } else {
    calibrated <- .pap_calibrate(GenoM, LifeHistData, tables, settings, progress)
    passes <- calibrated$passes
    final <- calibrated$final$pedigree
    final_tables <- calibrated$final_tables
    run <- calibrated$final
  }
  if (isTRUE(settings$CalcLLR)) {
    .pap_status(progress, "final LLRs")
    llr_started <- proc.time()[["elapsed"]]
    final <- .pap_final_llr(final, GenoM, LifeHistData, final_tables, settings)
    .pap_say(progress, sprintf("final LLRs (%.1f s)", proc.time()[["elapsed"]] - llr_started))
  }
  .pap_check_removed_path()
  out <- final[, intersect(c("id", "dam", "sire", "LLRdam", "LLRsire", "LLRpair",
                             "OHdam", "OHsire", "MEpair", "DamLLRStatus",
                             "SireLLRStatus", "PairLLRStatus"), names(final)), drop = FALSE]
  for (field in c("final_llr_evidence", "final_llr_hypotheses", "final_llr_unavailable",
                  "log_base", "internal_log_base"))
    attr(out, field) <- attr(final, field, exact = TRUE)
  rownames(out) <- NULL
  attr(out, "calibration") <- .pap_calibration_table(passes)
  attr(out, "calibration_passes") <- passes
  attr(out, "converged") <- isTRUE(run$converged)
  attr(out, "rounds") <- run$rounds
  attr(out, "MaxMismatch") <- run$MaxMismatch
  attr(out, "MaxMendelE") <- run$MaxMendelE
  attr(out, "unavailable") <- run$unavailable
  attr(out, "seconds") <- proc.time()[["elapsed"]] - started
  attr(out, "settings") <- settings
  .pap_say(progress, sprintf("done: %d dams, %d sires in %s", sum(!is.na(out$dam)),
                             sum(!is.na(out$sire)), if (is.null(progress)) "" else .pap_clock(progress)))
  if (!isTRUE(run$converged))
    warning("the final search reached max_rounds before converging", call. = FALSE)
  out
}

#' Settings for assign_parents
#'
#' Returns the complete settings list used by \code{\link{assign_parents}},
#' with any supplied values replacing the defaults. Pass the result, or a
#' partial named list, as \code{control}.
#'
#' @param control named list of settings to change. Names:
#'   \describe{
#'     \item{\code{error_model}}{\code{"dropout_foreign"} (default, uses
#'       \code{Dropout} and \code{Foreign}) or \code{"sequoia"} (Sequoia's SNP
#'       model, uses \code{Err}).}
#'     \item{\code{Err}}{panel-wide genotype error for \code{"sequoia"}.}
#'     \item{\code{error_calibration}}{\code{"raise_guarded"} (default),
#'       \code{"raise_only"} or \code{"none"}. Raise guarded fits the rate on
#'       confident trios, raises it to the lower 95\% bound when higher, may
#'       lower a rate above twice its upper bound once on the first pass, and
#'       restarts the search, at most four searches.}
#'     \item{\code{frequencies}}{optional list of per-locus allele-frequency
#'       vectors; default estimated from \code{GenoM} by \code{MHStats}.}
#'     \item{\code{n_alleles}}{optional per-locus allele counts; default the
#'       largest allele code seen at each locus.}
#'     \item{\code{ForeignQ}}{optional per-locus foreign-allele distributions;
#'       default the allele frequencies.}
#'     \item{\code{Tassign}}{minimum base-10 LLR to assign a parent (0.5).}
#'     \item{\code{margin}}{base-10 identity gap required over the best
#'       configuration with different parents (0.5).}
#'     \item{\code{Tfilter}}{pedigree-free parent-offspring versus unrelated
#'       screen, base 10 (-2); \code{-Inf} disables it.}
#'     \item{\code{prescreen_K}, \code{prescreen_margin_ln}}{support-only
#'       pre-screen: keep up to K candidates per side within the natural-log
#'       margin below Tassign (10, \code{log(100)}). Larger K trades runtime for
#'       recall; \code{Inf} disables the pre-screen.}
#'     \item{\code{MinAgeParent}, \code{MaxAgeParent}}{one or two (dam, sire)
#'       parental age bounds in years (1; default span of the birth years).}
#'     \item{\code{T}, \code{MaxMendelE}}{opposite-homozygote and trio
#'       Mendelian-error ceilings; default calibrated from the tables.}
#'     \item{\code{retain}}{true-parent retention target for those ceilings;
#'       default \code{0.9999^(1/N)}.}
#'     \item{\code{pedigree}}{optional starting pedigree (\code{id}, \code{dam},
#'       \code{sire}) of genotyped parents.}
#'     \item{\code{max_rounds}}{maximum search rounds (41).}
#'     \item{\code{CalcLLR}}{add read-only final parent LLRs (\code{FALSE}).}
#'     \item{\code{evidence}}{\code{"ancestral"} (default) or the experimental
#'       \code{"full"}, which also uses descendant evidence.}
#'     \item{\code{llr_evidence}}{evidence for final LLRs (\code{CalcLLR} and
#'       \code{\link{calc_llr}}): \code{"ancestral"} (default, the evidence
#'       the search decides on) or \code{"full"}, which adds descendants
#'       through forest messages on the fixed pedigree. Independent of
#'       \code{evidence}.}
#'     \item{\code{hypotheses}}{relationship suite: \code{"S7"} (default) or the
#'       experimental \code{"S37"} and \code{"S45"}.}
#'     \item{\code{memory_mib}, \code{work_millions}}{per-calculation resource
#'       limits: memory 64 (default), 128, 256 or 512 MiB; work 400 (default),
#'       800, 1600, 2000 or 5000 million operations. A calculation above them is
#'       unavailable, never zero evidence.}
#'     \item{\code{verbose}}{show progress (\code{FALSE}): a status line with
#'       the search, round and percent of the round done, redrawn every two
#'       seconds, then one line per search with its fit and decision.
#'       Messages go through \code{message()}.}
#'   }
#' @return A named list of all settings.
#' @examples
#' names(pap_control())
#' s <- pap_control(list(Tassign = 1, MaxAgeParent = 4, verbose = TRUE))
#' s[c("Tassign", "MaxAgeParent", "error_calibration")]
#' @export
pap_control <- function(control = list()) {
  defaults <- list(
    error_model = "dropout_foreign", Err = NULL,
    error_calibration = "raise_guarded",
    frequencies = NULL, n_alleles = NULL, ForeignQ = NULL,
    Tassign = 0.5, margin = 0.5, Tfilter = -2,
    prescreen_K = 10L, prescreen_margin_ln = log(100),
    MinAgeParent = 1L, MaxAgeParent = NULL,
    T = NULL, MaxMendelE = NULL, retain = NULL,
    pedigree = NULL, max_rounds = 41L, CalcLLR = FALSE,
    evidence = "ancestral", llr_evidence = "ancestral", hypotheses = "S7",
    memory_mib = 64, work_millions = 400, verbose = FALSE)
  if (is.null(control)) control <- list()
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("control must be a named list", call. = FALSE)
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown))
    stop("unknown control setting(s): ", paste(unknown, collapse = ", "),
         ". Valid names: ", paste(names(defaults), collapse = ", "), call. = FALSE)
  s <- defaults
  s[names(control)] <- control
  # A resolved settings list carries the names the user actually set.
  supplied <- attr(control, "supplied", exact = TRUE) %||% names(control)
  one_of <- function(name, values) {
    v <- s[[name]]
    if (!is.character(v) || length(v) != 1L || !v %in% values)
      stop("control$", name, " must be one of: ", paste(values, collapse = ", "),
           call. = FALSE)
  }
  one_of("error_model", c("dropout_foreign", "sequoia"))
  one_of("error_calibration", c("raise_guarded", "raise_only", "none"))
  one_of("evidence", c("ancestral", "full"))
  one_of("llr_evidence", c("ancestral", "full"))
  one_of("hypotheses", c("S7", "S37", "S45"))
  number <- function(name, lower = -Inf, allow_inf = FALSE) {
    v <- s[[name]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < lower ||
        (!allow_inf && !is.finite(v)))
      stop("control$", name, " must be one number", if (is.finite(lower))
        paste0(" of at least ", lower), call. = FALSE)
  }
  number("Tassign"); number("margin", 0); number("prescreen_margin_ln", 0)
  if (!isTRUE(s$memory_mib %in% c(64, 128, 256, 512)) || length(s$memory_mib) != 1L)
    stop("control$memory_mib must be 64, 128, 256 or 512", call. = FALSE)
  if (!isTRUE(s$work_millions %in% c(400, 800, 1600, 2000, 5000)) || length(s$work_millions) != 1L)
    stop("control$work_millions must be 400, 800, 1600, 2000 or 5000", call. = FALSE)
  if (!is.numeric(s$Tfilter) || length(s$Tfilter) != 1L || is.na(s$Tfilter) ||
      (is.infinite(s$Tfilter) && s$Tfilter > 0))
    stop("control$Tfilter must be one number or -Inf", call. = FALSE)
  K <- s$prescreen_K
  if (!is.numeric(K) || length(K) != 1L || is.na(K) || K <= 0 ||
      (is.finite(K) && (K != floor(K) || K > 100000)))
    stop("control$prescreen_K must be a positive whole number (at most 100000) or Inf",
         call. = FALSE)
  for (name in c("max_rounds")) {
    v <- s[[name]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 1 || v != floor(v))
      stop("control$", name, " must be one positive whole number", call. = FALSE)
  }
  for (name in c("CalcLLR", "verbose"))
    if (!is.logical(s[[name]]) || length(s[[name]]) != 1L || is.na(s[[name]]))
      stop("control$", name, " must be TRUE or FALSE", call. = FALSE)
  if (identical(s$error_model, "sequoia") && is.null(s$Err))
    stop("control$error_model = \"sequoia\" requires control$Err", call. = FALSE)
  if (!identical(s$error_model, "sequoia") && !is.null(s$Err))
    stop("control$Err applies only to error_model = \"sequoia\"", call. = FALSE)
  attr(s, "supplied") <- supplied
  s
}

# Locus tables from the genotypes and the error model. Frequencies come from
# MHStats on all typed individuals unless supplied.
.pap_tables <- function(GenoM, Dropout, Foreign, settings) {
  freqs <- settings$frequencies
  if (is.null(freqs))
    freqs <- MHStats(GenoM, n_alleles = settings$n_alleles,
                     warn_call_rate = 0, warn_Ae = 0)[["ALL"]]
  if (identical(settings$error_model, "sequoia")) {
    if (!is.null(Dropout) || !is.null(Foreign))
      stop("Dropout and Foreign apply only to error_model = \"dropout_foreign\"",
           call. = FALSE)
    return(build_tables(freqs, Err = settings$Err, ErrorModel = "legacy_sequoia"))
  }
  if (is.null(Dropout) || is.null(Foreign))
    stop("Dropout and Foreign are required (or use control$error_model = \"sequoia\")",
         call. = FALSE)
  build_tables(freqs, ErrorModel = "dropout_foreign", Dropout = Dropout,
               Foreign = Foreign, ForeignQ = settings$ForeignQ)
}

# Native functions that no supported input reaches were replaced by stubs that
# record their name (src/pap_removed_path.h). Stop rather than return results
# from a search that reached one, even if the native caller caught the error.
.pap_check_removed_path <- function() {
  hit <- pap_removed_path_hit_cpp()
  if (!is.na(hit))
    stop("internal error: PaPeR reached removed native code (", hit,
         "). Please report it with the input that caused it.", call. = FALSE)
  invisible(NULL)
}

# Locus tables for the pair utilities: only the error-model and frequency
# settings of the control list are used.
.pap_utility_tables <- function(GenoM, Dropout, Foreign, control) {
  pap_removed_path_hit_cpp()
  tables <- .pap_tables(GenoM, Dropout, Foreign, pap_control(control))
  .require_positive_error_rates(tables)
  tables
}

.pap_current_rates <- function(tables) {
  t <- tables[[1L]]
  if (identical(t$ErrorModel, "dropout_foreign"))
    c(Dropout = t$Dropout, Foreign = t$Foreign) else c(Err = t$Err)
}

# pedigree likelihood execution inputs from the user settings. Everything else in the native
# input comes from .engine_profile() and .engine_defaults().
.engine_controls <- function(settings, registry_generation) {
  # Settings left at their defaults are not passed, so each evidence mode keeps
  # its own validated pre-screen and limit defaults.
  supplied <- attr(settings, "supplied") %||% character()
  out <- list(pedigree_hypothesis_suite = settings$hypotheses,
              likelihood_evidence = settings$evidence,
              likelihood_registry_generation = registry_generation)
  if ("prescreen_K" %in% supplied) {
    K <- settings$prescreen_K
    out$pair_rescue <- is.finite(K)
    out$pair_rescue_top_k <- if (is.finite(K)) as.integer(K) else 100000L
  }
  if ("prescreen_margin_ln" %in% supplied)
    out$pair_rescue_margin_ln <- settings$prescreen_margin_ln
  if ("memory_mib" %in% supplied) out$likelihood_memory_mib <- settings$memory_mib
  if ("work_millions" %in% supplied) out$likelihood_work_millions <- settings$work_millions
  # The native planner requires this opt-in to admit expanded resource limits.
  if (settings$memory_mib != 64 || settings$work_millions != 400)
    out$likelihood_high_allele_programs <- TRUE
  out
}

# One parentage search at fixed tables. Returns the final pedigree and the
# records needed for error calibration.
.pap_search <- function(GenoM, LifeHistData, tables, settings,
                        registry_generation = .next_registry_generation(),
                        progress = NULL) {
  if (!is.null(progress)) progress$search <- progress$search + 1L
  .pap_status(progress, sprintf("search %d  preparing candidates",
                                if (is.null(progress)) 0L else progress$search))
  G <- .prepare_for_tables(GenoM, tables)
  ids <- rownames(G)
  n <- length(ids)
  retain <- .validate_retention_target(settings$retain %||% 0.9999^(1 / n))
  lh <- .normalize_life_history(LifeHistData, ids)
  .require_birth_bounds(lh)
  bounds <- .parent_bounds(lh, settings$MinAgeParent, settings$MaxAgeParent)
  minp <- bounds$min; maxp <- bounds$max
  ctx <- .new_inference_context(G, tables, "full")
  pm <- po_mismatch_prob(tables)
  T <- settings$T
  if (is.null(T)) T <- calibrate_maxmismatch(pm$p_fm, retain = retain)$T
  T <- as.integer(T)
  oh_calibration <- .calibrate_screen_threshold(pm$p_fm, retain, "OH", threshold = T)
  if (is.null(settings$MaxMendelE)) {
    me_calibration <- calibrate_max_mendelian_error(pm, retain = retain)
  } else {
    me_calibration <- .calibrate_screen_threshold(
      pm$p_me_true_trio, retain, "MaxMendelE",
      threshold = as.integer(settings$MaxMendelE))
    me_calibration$p_me_true_trio <- pm$p_me_true_trio
  }
  state <- .initialize_parentage_state(
    ids, genotype_rows = seq_along(ids), life_history = lh,
    min_age_parent = minp, max_age_parent = maxp,
    pedigree = settings$pedigree, ctx = ctx)
  universe <- .relationship_universe("full", "PO", operation = "parent")
  candidate_universe <- .relationship_universe("simp", "PO", operation = "parent")
  pool <- .build_candidate_pool(ctx, lh, ids, minp, maxp, T, settings$Tfilter,
                                candidate_universe, oh_calibration, me_calibration)
  .attach_parentage_candidate_pool(state, pool, ctx)
  snapshot <- .build_message_snapshot(ctx, state)
  state$message_snapshot <- snapshot
  state$current_messages <- if (is.null(snapshot)) NULL else snapshot$messages
  Tassign_ln <- .threshold_log10_to_ln(settings$Tassign)
  margin_ln <- .threshold_log10_to_ln(settings$margin)
  Tfilter_ln <- .threshold_log10_to_ln(settings$Tfilter)
  checkpoint <- NULL
  diagnostics <- list(
    complete_recheck = "first_dirty_wave_convergence",
    search_loop_engine = "bounded_messages",
    native_state_callback = function(event) {
      if (is.null(checkpoint) && identical(event$event, "start")) checkpoint <<- event
      invisible(NULL)
    })
  if (!is.null(progress)) diagnostics$visit_progress_callback <- .pap_visit_callback(progress)
  driver <- .run_search(
    ctx, state, pool, universe, snapshot, Tassign_ln, margin_ln,
    Tfilter_ln = Tfilter_ln, tol = 1e-8, max_rounds = as.integer(settings$max_rounds),
    shortlist_settings = .shortlist_settings("sequoia50", NULL, required_margin_ln = margin_ln),
    run_check_pedigree = !is.null(settings$pedigree),
    relationship_mode = "full", pass_diagnostics = diagnostics,
    likelihood_controls = .engine_controls(settings, registry_generation))
  .pap_check_removed_path()
  ped <- driver$pedigree
  ped <- data.frame(id = ped$id, dam = ped$dam, sire = ped$sire, stringsAsFactors = FALSE)
  native <- driver$native_result
  list(pedigree = ped, converged = isTRUE(driver$converged),
       rounds = driver$rounds, MaxMismatch = T, MaxMendelE = pool$MaxMendelE,
       graph_revision = driver$state$graph_revision,
       focal_summaries = c(lapply(driver$pass_results, `[[`, "focal_summary"),
                           lapply(driver$recheck_pass_results %||% list(), `[[`, "focal_summary")),
       complete_trio_record_state = native$complete_trio_record_state,
       unavailable = native$likelihood_availability,
       native_seconds = driver$native_engine_seconds,
       checkpoint = checkpoint)
}
