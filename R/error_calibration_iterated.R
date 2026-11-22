# Per-parameter direction memory persists through a retained parameter. A
# reversal in either parameter stops the whole proposal before another search.
.parentage_error_iteration_step <- function(fit, current, previous_direction, pass, direction = 'two_sided') {
  proposal <- .parentage_error_bound_decision(fit, current, direction)
  proposed <- setNames(as.numeric(proposal$used_rate), names(current))
  change <- sign(proposed - current)
  reversal <- change != 0 & previous_direction != 0 & change != previous_direction
  reason <- if (!proposal$adjust) proposal$status else if (any(reversal))
    'direction_reversal' else if (pass >= 4L) 'maximum_searches' else 'continue'
  continue <- identical(reason, 'continue')
  updated <- previous_direction
  if (continue) updated[change != 0] <- change[change != 0]
  list(proposal = proposal, proposed_rates = proposed, change = change,
    reversal_parameters = names(current)[reversal], stop_reason = reason,
    continue = continue, next_rates = if (continue) proposed else current,
    previous_direction = updated)
}

# The coverage guard uses all assigned full trios as its denominator, not only
# those admitted to the estimator. No downward proposal is iterated.
.parentage_guarded_lower_step <- function(fit, current, confident, assigned, pass) {
  admitted <- assigned > 0 && confident / assigned >= .7
  bounded <- fit
  if (identical(fit$status, 'fitted')) {
    bounded$interval <- fit$interval
    upper <- if (is.matrix(bounded$interval)) bounded$interval[names(current), 'upper'] else unname(bounded$interval['upper'])
    # The trigger stays strict factor-two; the target is the original upper
    # bound. Mask ineligible parameters so lower_only cannot lower them.
    target <- ifelse(is.finite(upper) & current > 2 * upper, upper, current)
    if (is.matrix(bounded$interval)) bounded$interval[names(current), 'upper'] <- target
    else bounded$interval['upper'] <- unname(target)
  }
  if (!admitted) bounded$status <- 'insufficient_confident_fraction'
  if (pass > 1L) bounded$status <- 'guarded_lowering_completed'
  step <- .parentage_error_iteration_step(bounded, current,
    setNames(rep(0, length(current)), names(current)), pass, 'lower_only')
  # Keep reported confidence bounds as fitted, separate from the guarded target.
  if (identical(fit$status, 'fitted')) {
    upper <- if (is.matrix(fit$interval)) fit$interval[, 'upper'] else unname(fit$interval['upper'])
    step$proposal$upper_95 <- upper
    if (!is.null(step$proposal$parameters)) for (name in names(current))
      step$proposal$parameters[[name]]$upper_95 <- unname(upper[name])
  }

  step
}

# Raising remains available after each fit. A downward proposal is considered
# only on the original pass, with the coverage and factor-two guards. Later
# passes never lower; every restart still uses the original data and pedigree.
.parentage_error_raise_guarded_step <- function(fit, current, confident, assigned, pass) {
  step <- .parentage_error_iteration_step(fit, current,
    setNames(rep(0, length(current)), names(current)), pass, 'raise_only')
  if (pass != 1L) return(step)
  lower <- .parentage_guarded_lower_step(fit, current, confident, assigned, pass)
  proposed <- step$proposed_rates
  lowering <- lower$proposed_rates < current & proposed == current
  proposed[lowering] <- lower$proposed_rates[lowering]
  step$proposal$direction <- 'raise_guarded'
  step$proposal$used_rate <- if (length(current) == 1L) unname(proposed) else proposed
  step$proposal$adjust <- any(proposed != current)
  step$proposal$raise <- any(proposed > current)
  step$proposal$lower <- any(proposed < current)
  if (step$proposal$adjust) step$proposal$status <- 'parameters_adjusted_to_guarded_bounds'
  if (!is.null(step$proposal$parameters)) for (name in names(current)[lowering])
    step$proposal$parameters[[name]] <- lower$proposal$parameters[[name]]
  step$proposed_rates <- proposed; step$change <- sign(proposed - current)
  step$continue <- step$proposal$adjust
  step$stop_reason <- if (step$continue) 'continue' else step$proposal$status
  step$next_rates <- if (step$continue) proposed else current
  step$previous_direction <- pmax(0, step$change)
  step
}

# Error calibration around .pap_search(): fit the rate on the confident trios
# of each final pedigree, adjust it under the selected policy, and restart the
# search from the original inputs. At most four searches.
.pap_calibrate <- function(GenoM, LifeHistData, tables, settings, progress = NULL) {
  policy <- settings$error_calibration
  mh <- identical(tables[[1L]]$ErrorModel, "dropout_foreign")
  if (mh) {
    current <- c(Dropout = tables[[1L]]$Dropout, Foreign = tables[[1L]]$Foreign)
    rebuild <- function(rate) .parentage_mh_rebuild(tables, rate)
  } else {
    if (!all(vapply(tables, function(t) t$k == 2L, logical(1))))
      stop("error calibration with error_model = \"sequoia\" requires biallelic loci",
           call. = FALSE)
    spec <- .parentage_error_scale_spec(tables)
    current <- c(Err = spec$initial)
    rebuild <- function(rate) spec$rebuild(unname(rate))
  }
  basis_cache <- if (mh) {
    e <- new.env(parent = emptyenv()); e$limit <- 8 * 1024^2
    e$hits <- e$builds <- e$declines <- 0L; e$seconds <- e$bytes <- 0; e
  } else NULL
  generation <- .next_registry_generation()
  previous <- setNames(rep(0, length(current)), names(current))
  passes <- list(); run <- NULL; pass_tables <- tables
  for (pass in 1:4) {
    pass_tables <- if (pass == 1L) tables else rebuild(current)
    started <- proc.time()[["elapsed"]]
    run <- .pap_search(GenoM, LifeHistData, pass_tables, settings, generation, progress)
    search_seconds <- proc.time()[["elapsed"]] - started
    fit_started <- proc.time()[["elapsed"]]
    .pap_status(progress, sprintf("search %d  fitting the error rate", pass))
    confidence <- NULL
    if (!run$converged || is.null(run$checkpoint)) {
      fit <- list(status = if (!run$converged) "pass_not_converged" else "missing_native_checkpoint")
    } else {
      confidence <- .parentage_final_confidence(
        run$checkpoint$compact, run$checkpoint$input, run$pedigree,
        run$focal_summaries, TRUE, final_revision = run$graph_revision,
        dependency_state = run$complete_trio_record_state)
      fit_tables <- pass_tables
      if (!is.null(basis_cache)) attr(fit_tables, "pap_category_cache") <- basis_cache
      fit <- estimate_parentage_error(GenoM, fit_tables, run$pedigree, run$MaxMendelE,
                                      confident = confidence$confident)
      fit <- .parentage_error_ceiling_audit(fit, run$MaxMismatch)
    }
    fit$tables <- NULL
    fit_seconds <- proc.time()[["elapsed"]] - fit_started
    parents <- run$pedigree
    assigned <- sum(!is.na(parents$dam) & nzchar(parents$dam) &
                    !is.na(parents$sire) & nzchar(parents$sire))
    confident <- if (is.null(confidence)) 0L else sum(confidence$confident)
    step <- if (policy == "raise_guarded")
      .parentage_error_raise_guarded_step(fit, current, confident, assigned, pass) else
      .parentage_error_iteration_step(fit, current, previous, pass, "raise_only")
    record <- .pap_pass_record(pass, current, run, fit, step, confidence)
    record$assigned_trios <- assigned
    record$confident_fraction <- if (assigned) confident / assigned else NA_real_
    record$search_seconds <- search_seconds
    record$fit_seconds <- fit_seconds
    passes[[pass]] <- record
    .pap_report_search(progress, record, search_seconds + fit_seconds)
    if (!step$continue) break
    previous <- step$previous_direction
    current <- step$next_rates
  }
  list(passes = passes, final = run, final_tables = pass_tables)
}

.pap_pass_record <- function(pass, rates, run, fit, step, confidence) {
  list(pass = pass, rates_used = rates, parents = run$pedigree, fit = fit,
       confident_trios = if (is.null(confidence)) NA_integer_ else sum(confidence$confident),
       stop_reason = step$stop_reason %||% step$status %||% NA_character_,
       next_rates = step$next_rates %||% rates,
       converged = run$converged, rounds = run$rounds,
       max_mendelian_errors = run$MaxMendelE, unavailable = run$unavailable,
       native_seconds = run$native_seconds, step = step, confidence = confidence)
}

# One row per search: the rate used, the fit and its 95% profile bounds.
.pap_calibration_table <- function(passes) {
  rows <- lapply(passes, function(p) {
    rates <- p$rates_used
    fit <- p$fit %||% list()
    est <- fit$estimated_rates %||% fit$estimate
    iv <- fit$interval
    out <- data.frame(search = p$pass, stringsAsFactors = FALSE)
    for (name in names(rates)) {
      out[[paste0(name, "_used")]] <- unname(rates[[name]])
      e <- if (is.null(est)) NA_real_ else if (length(est) == 1L) unname(est) else unname(est[[name]])
      lo <- hi <- NA_real_
      if (is.matrix(iv) && name %in% rownames(iv)) { lo <- iv[name, "lower"]; hi <- iv[name, "upper"] }
      else if (!is.null(iv) && !is.matrix(iv) && length(rates) == 1L) { lo <- iv[["lower"]]; hi <- iv[["upper"]] }
      out[[paste0(name, "_estimate")]] <- e
      out[[paste0(name, "_lower")]] <- unname(lo)
      out[[paste0(name, "_upper")]] <- unname(hi)
    }
    out$confident_trios <- p$confident_trios
    out$decision <- p$stop_reason
    out$converged <- p$converged
    out$rounds <- p$rounds
    out
  })
  do.call(rbind, rows)
}
