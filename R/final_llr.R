#' Score the assigned parents of a fixed pedigree
#'
#' Computes read-only pedigree likelihood relationship LLRs with the supplied error rates.
#' No search or error calibration is run. Each assigned parent is compared
#' with the nonparent alternatives in \code{control$hypotheses}, keeping the
#' other assigned parent. \code{control$llr_evidence} selects ancestral evidence
#' (the default) or full forest evidence, including descendants. Forest cuts
#' approximate correlations outside the local question in loopy pedigrees.
#'
#' @param Pedigree data frame with \code{id}, \code{dam}, \code{sire}. All IDs
#'   and assigned parents must be genotyped. Omitted genotyped rows are founders.
#' @inheritParams assign_parents
#' @return The supplied pedigree with base-10 \code{LLRdam}, \code{LLRsire},
#'   \code{LLRpair} (the smaller side, only when both parents are assigned),
#'   \code{OHdam}, \code{OHsire}, \code{MEpair}, and \code{DamLLRStatus},
#'   \code{SireLLRStatus}, \code{PairLLRStatus}. Status is \code{"available"},
#'   \code{"not_assigned"}, or \code{"unavailable"}. Unavailable calculations
#'   and unassigned sides have \code{NA} LLRs. Row order and parent links are
#'   preserved. Attributes \code{final_llr_evidence}, \code{final_llr_hypotheses}
#'   and \code{final_llr_unavailable} record the scoring settings and refusals.
#' @examples
#' ex <- pap_example
#' # Score the true pedigree of the first two cohorts, keeping only
#' # genotyped parents
#' first <- ex$LifeHistData$ID[ex$LifeHistData$BirthYear <= 2016]
#' truth <- ex$Pedigree[ex$Pedigree$id %in% first, c("id", "dam", "sire")]
#' truth$dam[!truth$dam %in% first] <- NA
#' truth$sire[!truth$sire %in% first] <- NA
#' scored <- calc_llr(truth, ex$GenoM[first, ], ex$LifeHistData,
#'                    Dropout = 0.005, Foreign = 5e-4)
#' summary(scored$LLRdam)
#' @export
calc_llr <- function(Pedigree, GenoM, LifeHistData, Dropout = NULL, Foreign = NULL,
                     control = list()) {
  settings <- pap_control(control)
  pap_removed_path_hit_cpp()
  tables <- .pap_tables(GenoM, Dropout, Foreign, settings)
  .require_positive_error_rates(tables)
  out <- .pap_final_llr(Pedigree, GenoM, LifeHistData, tables, settings)
  .pap_check_removed_path()
  out
}

# A fresh context owns the fixed pedigree; no search checkpoint or warmed
# search messages are used. The empty discovery pool makes the shared input
# builder include only incumbents, without performing candidate discovery.
.pap_final_llr <- function(Pedigree, GenoM, LifeHistData, tables, settings,
                           batch_size = 16L) {
  if (!is.data.frame(Pedigree) ||
      !all(c("id", "dam", "sire") %in% names(Pedigree)))
    stop("Pedigree must be a data.frame with id, dam, and sire columns", call. = FALSE)
  G <- .prepare_for_tables(GenoM, tables)
  ids <- rownames(G)
  lh <- .normalize_life_history(LifeHistData, ids)
  .require_birth_bounds(lh)
  ctx <- .new_inference_context(G, tables, "full")
  state <- .initialize_parentage_state(
    ids, life_history = lh, min_age_parent = settings$MinAgeParent,
    max_age_parent = settings$MaxAgeParent, pedigree = Pedigree, ctx = ctx)
  p <- state$parents
  out <- Pedigree
  for (field in c("LLRdam", "LLRsire", "LLRpair")) out[[field]] <- rep(NA_real_, nrow(out))
  for (field in c("OHdam", "OHsire", "MEpair")) out[[field]] <- rep(NA_integer_, nrow(out))
  for (field in c("DamLLRStatus", "SireLLRStatus", "PairLLRStatus"))
    out[[field]] <- rep("not_assigned", nrow(out))
  # Only supplied output rows need scores; other genotyped individuals still
  # contribute their observations through the fixed pedigree's messages.
  at <- match(as.character(out$id), ids)
  focals <- at[rowSums(p[at, , drop = FALSE] > 0L) > 0L]
  availability <- NULL
  if (length(focals)) {
    snapshot <- .build_message_snapshot(ctx, state, purpose = "final-output-scoring")
    # Final LLRs have their own evidence setting; the search setting is not used.
    scoring <- settings; scoring$evidence <- settings$llr_evidence
    universe <- .relationship_universe("full", "PO", operation = "parent")
    pool <- list(rows = .candidate_pool_empty_rows(), shortlist_rows = list(),
                 MaxMendelE = as.integer(ncol(G) / 2L))
    setup <- .prepare_search(
      state, pool, universe, snapshot, .threshold_log10_to_ln(settings$Tassign),
      .threshold_log10_to_ln(settings$margin), -Inf, 1e-8, 1L,
      .shortlist_settings("sequoia50", NULL),
      .engine_controls(scoring, .next_registry_generation()),
      run_check_pedigree = FALSE, relationship_mode = "full")
    input <- setup$input
    # Disable search screens and request complete side contrasts even for a
    # weak supplied assignment. Resource limits and biological admission stay.
    input$pair_rescue <- FALSE
    for (side in c("dam", "sire")) {
      input[[paste0(side, "_quick")]] <- rep(0, length(input[[paste0(side, "_ids")]]))
      input[[paste0(side, "_permanent")]] <- rep(TRUE, length(input[[paste0(side, "_ids")]]))
    }
    input$private_frozen_decision_audit <- TRUE
    input$private_decision_only_diagnostics <- TRUE
    input$configuration_diagnostics <- TRUE
    input$configuration_locus_diagnostics <- FALSE
    input$diagnostic_start_position <- 1L; input$diagnostic_end_position <- 0L
    input$recheck_only <- TRUE
    input$focal_start_position <- input$focal_end_position <- 1L

    input$start_graph_revision <- input$message_revision <- 1L
    input$diagnostic_visit_recheck_order <- NULL
    invisible(publish_messages_cpp(snapshot$compiled, p, 1L, 1L, 1L))
    batches <- split(sort(focals), ceiling(seq_along(focals) / batch_size))
    for (batch in batches) {
      input$diagnostic_visit_order <- input$recheck_focals <- as.integer(batch)
      z <- run_parentage_cpp(snapshot$compiled, setup$compact, input)
      if (!isTRUE(z$frozen_state_audit) ||
          !identical(unname(z$final_parents), unname(p)) || nrow(z$accepted_moves))
        stop("internal error: final LLR audit mutated the pedigree", call. = FALSE)
      .pap_check_removed_path()
      if (is.numeric(z$likelihood_availability)) availability <-
        if (is.null(availability)) z$likelihood_availability else availability + z$likelihood_availability
      d <- as.data.frame(z$configuration_diagnostics)
      for (f in batch) {
        i <- match(f, at)
        r <- d[d$focal == f & d$dam == p[f, 1L] & d$sire == p[f, 2L], , drop = FALSE]
        for (side in seq_len(2L)) {
          if (p[f, side] == 0L) next
          label <- c("dam", "sire")[side]
          status <- c("DamLLRStatus", "SireLLRStatus")[side]
          value <- if (nrow(r) == 1L) r[[paste0("full_", label, "_relationship_ln")]] else NA_real_
          out[[status]][i] <- if (is.na(value)) "unavailable" else "available"
          out[[paste0("LLR", label)]][i] <- .llr_log10(value)
        }
      }
    }
  }
  for (i in seq_len(nrow(out))) {
    f <- at[i]; d <- p[f, 1L]; s <- p[f, 2L]
    counts <- count_oh_me_trio_cpp(G, f - 1L, d - 1L, s - 1L)
    if (d > 0L) out$OHdam[i] <- counts[["OHdam"]]
    if (s > 0L) out$OHsire[i] <- counts[["OHsire"]]
    if (d > 0L && s > 0L) {
      out$MEpair[i] <- counts[["MEpair"]]
      out$LLRpair[i] <- min(out$LLRdam[i], out$LLRsire[i])
      out$PairLLRStatus[i] <- if (is.na(out$LLRpair[i])) "unavailable" else "available"
    }
  }
  attr(out, "final_llr_evidence") <- settings$llr_evidence
  attr(out, "final_llr_hypotheses") <- settings$hypotheses
  attr(out, "final_llr_unavailable") <- availability
  attr(out, "log_base") <- .public_log_base
  attr(out, "internal_log_base") <- "e"
  out
}
