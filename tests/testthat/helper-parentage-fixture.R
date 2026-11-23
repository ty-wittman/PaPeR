# Shared literal-call native scientific fixture. This file executes no tests
# and constructs no assignment until the helper is explicitly called.

.scientific_fixture <- function(calls, dam = NULL, sire = NULL,
                                  dam_candidates = character(),
                                  sire_candidates = character(),
                                  focal = "A", p = c(.5, .5), err = .01,
                                  Complex = "simp", birth_year = NULL) {
  calls <- as.matrix(calls)
  ids <- rownames(calls)
  if (is.null(colnames(calls))) colnames(calls) <- "L1"
  if (is.null(birth_year)) {
    birth_year <- setNames(rep(0L, length(ids)), ids)
    birth_year[[focal]] <- 5L
  }
  parent_vector <- function(x) {
    out <- setNames(rep(NA_character_, length(ids)), ids)
    if (length(x)) out[names(x)] <- as.character(x)
    out
  }
  dam <- parent_vector(dam)
  sire <- parent_vector(sire)
  frequencies <- if (is.list(p)) p else rep(list(p), ncol(calls))
  names(frequencies) <- colnames(calls)
  tables <- build_tables(frequencies, Err = err)
  geno <- PaPeR:::.prepare_for_tables(calls, tables)
  sex <- rep(NA_integer_, length(ids))
  sex[ids %in% unique(na.omit(dam))] <- 1L
  sex[ids %in% unique(na.omit(sire))] <- 2L
  sex[ids %in% dam_candidates] <- 1L
  sex[ids %in% sire_candidates] <- 2L
  life <- data.frame(ID = ids, Sex = sex,
                     BirthYear = unname(as.integer(birth_year[ids])))
  pedigree <- data.frame(id = ids, dam = unname(dam), sire = unname(sire))
  lh <- PaPeR:::.normalize_life_history(life, ids)
  bounds <- PaPeR:::.parent_bounds(lh, 1L, NULL)
  ctx <- PaPeR:::.new_inference_context(geno, tables, Complex)
  state <- PaPeR:::.initialize_parentage_state(
    ids, seq_along(ids), lh, bounds$min, bounds$max, pedigree, ctx)
  universe <- PaPeR:::.relationship_universe(Complex, "PO", operation = "parent")
  mismatch <- po_mismatch_prob(tables)
  oh <- PaPeR:::.calibrate_screen_threshold(
    mismatch$p_fm, retain = .999, threshold = 99L)
  me <- PaPeR:::.calibrate_screen_threshold(
    mismatch$p_me_true_trio, retain = .999, threshold = 99L)
  me$p_me_true_trio <- mismatch$p_me_true_trio
  pool <- PaPeR:::.build_candidate_pool(
    ctx, lh, ids, bounds$min, bounds$max, 99L, -Inf, universe, oh, me)
  PaPeR:::.attach_parentage_candidate_pool(state, pool, ctx)
  # Compact construction preserves structural -Inf entries.
  snapshot <- PaPeR:::.build_message_snapshot(ctx, state)
  threshold <- .5 * log(10)
  input <- PaPeR:::.engine_inputs(
    snapshot, state, pool, ids, universe, threshold, threshold, -Inf,
    1e-8, 1L, PaPeR:::.shortlist_settings(
      "sequoia50", required_margin_ln = threshold), FALSE,
    relationship_mode = Complex, diagnostics = list(
      configuration_diagnostics = TRUE,
      configuration_locus_diagnostics = TRUE))
  focal_index <- as.integer(match(focal, ids))
  focal_position <- as.integer(match(focal_index, input$focal_order))
  input$focal_start_position <- input$focal_end_position <- focal_position
  input$configuration_locus_focal_position <- focal_position
  input$replay_moves <- matrix(integer(), nrow = 0L, ncol = 7L)
  # Restrict this mathematical query's CSR to its declared candidate set.
  # The retained background pedigree remains intact, including co-parents.
  for (side in c("dam", "sire")) {
    candidates <- if (side == "dam") dam_candidates else sire_candidates
    candidates <- as.integer(match(candidates, ids))
    stopifnot(!anyNA(candidates), !anyDuplicated(candidates))
    offsets <- integer(length(ids) + 1L)
    offsets[seq.int(focal_index + 1L, length(ids) + 1L)] <- length(candidates)
    input[[paste0(side, "_offsets")]] <- offsets
    input[[paste0(side, "_ids")]] <- candidates
    input[[paste0(side, "_quick")]] <- rep(0, length(candidates))
    input[[paste0(side, "_guard")]] <- rep(-Inf, length(candidates))
    input[[paste0(side, "_permanent")]] <- rep(TRUE, length(candidates))
  }
  # Also expose only the literal inputs needed to exercise the supported
  # public production routes. The native adapter inputs above remain private
  # to legacy tests and are not passed into assign_parents()/calc_llr().
  public_control <- list(error_model = "sequoia", Err = err,
                         frequencies = frequencies,
                         error_calibration = "none")
  list(calls = calls, ids = ids, tables = tables, ctx = ctx,
       snapshot = snapshot, compact_inputs =
         PaPeR:::.message_compact_snapshot_inputs(snapshot),
       input = input, focal = focal, dam = dam, sire = sire,
       life_history = life, pedigree = pedigree,
       frequencies = frequencies, public_control = public_control)
}

# Public-path companions used by adapter-free scientific regressions. They
# accept only supported pap_control fields and never set native diagnostic or
# probe flags.
.public_assignment <- function(f, control = list()) {
  ctl <- utils::modifyList(f$public_control, control)
  assign_parents(f$calls, f$life_history, control = ctl)
}

.public_fixed_score <- function(f, pedigree = f$pedigree,
                                       control = list()) {
  ctl <- utils::modifyList(f$public_control, control)
  calc_llr(pedigree, f$calls, f$life_history, control = ctl)
}

# Exact/global reference oracle for the supported fixed-pedigree public LLR.
# It scores every admitted rival in the same retained route family; callers
# can separately compare a simple PO-versus-U formula to those route scores.
