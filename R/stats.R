# Multiallelic marker statistics (Sequoia SnpStats analogue for microhaplotypes).

.mhstats_threshold <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0)
    stop(name, " must be one finite nonnegative number", call. = FALSE)
  as.numeric(x)
}

.mhstats_pseudocount <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0)
    stop("pseudocount must be one finite nonnegative number", call. = FALSE)
  as.numeric(x)
}

.mhstats_group_value <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  if (is.integer(x)) return(as.character(x))
  if (is.numeric(x)) {
    if (any(!is.na(x) & (!is.finite(x) | x != floor(x))))
      stop(name, " must contain characters, factors, or whole numbers", call. = FALSE)
    return(ifelse(is.na(x), NA_character_, format(x, scientific = FALSE, trim = TRUE)))
  }
  if (!is.character(x))
    stop(name, " must contain characters, factors, or whole numbers", call. = FALSE)
  x[!is.na(x) & !nzchar(x)] <- NA_character_
  x
}

.mhstats_sex <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  raw <- as.character(x)
  num <- suppressWarnings(as.numeric(raw))
  missing <- is.na(x) | !nzchar(raw)
  bad <- !missing & (is.na(num) | !is.finite(num) | num != floor(num) | !num %in% 1:3)
  if (any(bad))
    stop("Sex must be 1=female, 2=male, or 3/NA=unknown", call. = FALSE)
  sx <- as.integer(num)
  sx[sx == 3L] <- NA_integer_
  sx
}

.mhstats_metadata <- function(by, ids, strata) {
  if (is.null(by)) {
    if (!is.null(strata))
      stop("strata requires an ID-matched by data frame or matrix", call. = FALSE)
    return(list(data = NULL, strata = character(0), excluded = 0L))
  }
  if (!is.data.frame(by) && !is.matrix(by))
    stop("by must be a data frame or matrix", call. = FALSE)
  by <- as.data.frame(by, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(names(by)) || !"ID" %in% names(by))
    stop("by must contain an ID column", call. = FALSE)
  mid <- as.character(by$ID)
  if (anyNA(mid) || any(!nzchar(mid)) || anyDuplicated(mid))
    stop("by IDs must be unique, non-missing strings", call. = FALSE)
  missing_ids <- setdiff(ids, mid)
  extra_ids <- setdiff(mid, ids)
  if (length(missing_ids) || length(extra_ids)) {
    detail <- c(
      if (length(missing_ids))
        paste0("missing: ", paste(utils::head(missing_ids, 5L), collapse = ", ")),
      if (length(extra_ids))
        paste0("extra: ", paste(utils::head(extra_ids, 5L), collapse = ", "))
    )
    stop("by IDs must exactly match GenoM IDs (", paste(detail, collapse = "; "), ")",
         call. = FALSE)
  }
  by <- by[match(ids, mid), , drop = FALSE]

  has_cohort <- "Cohort" %in% names(by)
  has_birth <- "BirthYear" %in% names(by)
  cohort <- birth <- NULL
  if (has_cohort) cohort <- .mhstats_group_value(by$Cohort, "Cohort")
  if (has_birth) birth <- .mhstats_group_value(by$BirthYear, "BirthYear")
  if (has_cohort && has_birth) {
    same <- (is.na(cohort) & is.na(birth)) |
      (!is.na(cohort) & !is.na(birth) & cohort == birth)
    if (!all(same))
      stop("Cohort and BirthYear contain conflicting values", call. = FALSE)
  }
  canonical_cohort <- if (has_cohort) cohort else birth
  sex <- if ("Sex" %in% names(by)) .mhstats_sex(by$Sex) else NULL

  if (is.null(strata)) {
    requested <- c(if (!is.null(canonical_cohort)) "Cohort",
                   if (!is.null(sex)) "Sex")
    if (!length(requested))
      stop("by must contain Cohort, BirthYear, or Sex", call. = FALSE)
  } else {
    if (!is.character(strata) || !length(strata) || anyNA(strata))
      stop("strata must name Cohort, BirthYear, and/or Sex", call. = FALSE)
    bad <- setdiff(strata, c("Cohort", "BirthYear", "Sex"))
    if (length(bad))
      stop("unsupported strata: ", paste(bad, collapse = ", "), call. = FALSE)
    requested <- unique(ifelse(strata == "BirthYear", "Cohort", strata))
    requested <- c(intersect("Cohort", requested), intersect("Sex", requested))
  }
  if ("Cohort" %in% requested && is.null(canonical_cohort))
    stop("requested Cohort/BirthYear stratum is absent from by", call. = FALSE)
  if ("Sex" %in% requested && is.null(sex))
    stop("requested Sex stratum is absent from by", call. = FALSE)

  dat <- data.frame(ID = ids, stringsAsFactors = FALSE)
  if ("Cohort" %in% requested) dat$Cohort <- canonical_cohort
  if ("Sex" %in% requested) dat$Sex <- sex
  complete <- if (length(requested)) stats::complete.cases(dat[, requested, drop = FALSE])
              else rep(TRUE, length(ids))
  list(data = dat, strata = requested, excluded = sum(!complete), complete = complete)
}

.mhstats_pic <- function(p) {
  sq <- p^2
  pair_term <- sum(vapply(seq_len(max(0L, length(p) - 1L)), function(i) {
    if (i >= length(p)) return(0)
    sum(2 * sq[i] * sq[(i + 1L):length(p)])
  }, numeric(1)))
  value <- 1 - sum(sq) - pair_term
  if (abs(value) < 1e-15) 0 else value
}

.mhstats_slice <- function(G, rows, group, cohort, sex, n_alleles_all,
                           pseudocount, warn_call_rate, warn_Ae, strata) {
  loci <- attr(G, "locus_names")
  kvals <- attr(G, "n_alleles")
  nind <- length(rows)
  n_loci <- length(loci)
  empirical <- likelihood <- vector("list", n_loci)
  n_typed <- n_observed <- n_zero <- integer(n_loci)
  call_rate <- Ae <- PIC <- rep(NA_real_, n_loci)
  zero_adjusted <- logical(n_loci)

  for (l in seq_len(n_loci)) {
    a <- G[rows, 2L * l - 1L]
    b <- G[rows, 2L * l]
    typed <- a > 0L & b > 0L
    n_typed[l] <- sum(typed)
    counts <- tabulate(c(a[typed], b[typed]), nbins = kvals[l])
    n_observed[l] <- sum(counts > 0L)
    n_zero[l] <- sum(counts == 0L)
    call_rate[l] <- n_typed[l] / nind
    allele_names <- as.character(seq_len(kvals[l]))
    if (n_typed[l] == 0L) {
      empirical[[l]] <- stats::setNames(rep(NA_real_, kvals[l]), allele_names)
      likelihood[[l]] <- stats::setNames(rep(NA_real_, kvals[l]), allele_names)
      next
    }
    p <- counts / (2 * n_typed[l])
    empirical[[l]] <- stats::setNames(p, allele_names)
    Ae[l] <- 1 / sum(p^2)
    PIC[l] <- .mhstats_pic(p)
    if (n_zero[l] > 0L && pseudocount > 0) {
      q <- (counts + pseudocount * (counts == 0L)) /
        (2 * n_typed[l] + pseudocount * n_zero[l])
      likelihood[[l]] <- stats::setNames(q, allele_names)
      zero_adjusted[l] <- TRUE
    } else {
      likelihood[[l]] <- empirical[[l]]
    }
  }

  out <- data.frame(
    group = rep(group, n_loci),
    locus = loci,
    n_individuals = rep(nind, n_loci),
    n_typed = n_typed,
    call_rate = call_rate,
    n_alleles = n_observed,
    n_alleles_all = n_alleles_all,
    n_alleles_universe = kvals,
    Ae = Ae,
    PIC = PIC,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if ("Sex" %in% strata)
    out <- cbind(out[1L], Sex = rep(sex, n_loci), out[-1L],
                 stringsAsFactors = FALSE)
  if ("Cohort" %in% strata)
    out <- cbind(out[1L], Cohort = rep(cohort, n_loci), out[-1L],
                 stringsAsFactors = FALSE)
  out$allele_frequencies <- I(empirical)
  out$likelihood_frequencies <- I(likelihood)
  out$n_zero_alleles <- n_zero
  out$zero_adjusted <- zero_adjusted
  out$low_call_rate <- !is.na(call_rate) & call_rate < warn_call_rate
  out$low_Ae <- !is.na(Ae) & Ae < warn_Ae
  class(out) <- c("PaPeR_MHStats_slice", "data.frame")
  out
}

.mhstats_warning <- function(rows, label, threshold) {
  if (!nrow(rows)) return(invisible(NULL))
  keys <- paste0(rows$group, "/", rows$locus)
  shown <- paste(utils::head(keys, 5L), collapse = ", ")
  more <- if (length(keys) > 5L) paste0(" (and ", length(keys) - 5L, " more)") else ""
  warning(length(keys), " group-locus combination(s) have ", label, " below ",
          threshold, ": ", shown, more, call. = FALSE)
}

#' Microhaplotype allele-frequency and marker statistics
#'
#' The multiallelic analogue of Sequoia's \code{SnpStats}. It reports exact
#' empirical haplotype-allele frequencies and marker QC for the full sample and,
#' optionally, ID-matched Cohort (or BirthYear), Sex, or Cohort-by-Sex groups.
#' Cohort and BirthYear are aliases for the same grouping information.
#'
#' Empirical frequencies, effective allele number, and PIC are never
#' pseudocount-adjusted. The separate \code{likelihood_frequencies} list-column
#' adds \code{pseudocount} only to zero-count panel alleles, followed by
#' normalization, so likelihood construction does not assign zero probability
#' to a retained panel allele.
#'
#' @param GenoM public character genotype matrix or data frame, one column per
#'   locus, with calls \code{"a/b"} and missing genotypes \code{"0"}.
#' @param n_alleles optional positive integer vector with one value per locus.
#' @param pseudocount nonnegative amount added only to zero-count panel alleles
#'   in \code{likelihood_frequencies}. It never changes empirical summaries.
#' @param by optional ID-matched data frame or matrix containing \code{ID} and
#'   one or more of \code{Cohort}, \code{BirthYear}, and \code{Sex}.
#' @param strata optional character vector selecting \code{Cohort} (or its
#'   \code{BirthYear} alias), \code{Sex}, or both. When omitted, all supported
#'   grouping columns in \code{by} are used.
#' @param warn_call_rate warn when the proportion of individuals called at a
#'   group-locus combination is strictly below this value. Set to 0 to disable.
#' @param warn_Ae warn when effective allele number is strictly below this
#'   value. Set to 0 to disable.
#' @return A named list with \code{ALL} first and one data-frame slice per
#'   observed group. Each row is one locus and includes exact frequencies,
#'   likelihood-safe frequencies, call rate, observed allele counts, declared
#'   allele-universe sizes,
#'   \eqn{A_e=1/\sum p_i^2}, and codominant multiallelic PIC.
#' @examples
#' st <- MHStats(pap_example$GenoM, warn_Ae = 0)
#' head(st$ALL[, c("locus", "n_typed", "call_rate", "n_alleles", "Ae", "PIC")])
#' # by cohort
#' lh <- pap_example$LifeHistData
#' by_year <- MHStats(pap_example$GenoM, by = lh[, c("ID", "BirthYear")],
#'                    warn_Ae = 0)
#' names(by_year)
#' @export
MHStats <- function(GenoM, n_alleles = NULL, pseudocount = 0.5, by = NULL,
                    strata = NULL, warn_call_rate = 0.20, warn_Ae = 3) {
  G <- prepare_genotypes(GenoM, n_alleles)
  pseudocount <- .mhstats_pseudocount(pseudocount)
  warn_call_rate <- .mhstats_threshold(warn_call_rate, "warn_call_rate")
  if (warn_call_rate > 1)
    stop("warn_call_rate must not exceed 1", call. = FALSE)
  warn_Ae <- .mhstats_threshold(warn_Ae, "warn_Ae")
  ids <- rownames(G)
  meta <- .mhstats_metadata(by, ids, strata)

  all_counts <- vapply(seq_along(attr(G, "locus_names")), function(l) {
    a <- G[, 2L * l - 1L]
    b <- G[, 2L * l]
    sum(tabulate(c(a[a > 0L], b[b > 0L]),
                 nbins = attr(G, "n_alleles")[l]) > 0L)
  }, integer(1))
  out <- list(ALL = .mhstats_slice(
    G, seq_len(nrow(G)), "ALL", NA_character_, NA_integer_, all_counts,
    pseudocount, warn_call_rate, warn_Ae, meta$strata
  ))

  if (length(meta$strata)) {
    eligible <- which(meta$complete)
    dat <- meta$data
    cohort_codes <- if ("Cohort" %in% meta$strata)
      match(dat$Cohort, unique(dat$Cohort[eligible])) else rep(1L, nrow(dat))
    sex_codes <- if ("Sex" %in% meta$strata)
      match(dat$Sex, unique(dat$Sex[eligible])) else rep(1L, nrow(dat))
    keys <- paste(cohort_codes, sex_codes, sep = ":")
    first_keys <- unique(keys[eligible])
    labels <- vapply(first_keys, function(key) {
      i <- eligible[match(key, keys[eligible])]
      paste(c(if ("Cohort" %in% meta$strata) paste0("Cohort=", dat$Cohort[i]),
              if ("Sex" %in% meta$strata) paste0("Sex=", dat$Sex[i])),
            collapse = ";")
    }, character(1))
    labels <- make.unique(labels)
    for (j in seq_along(first_keys)) {
      rows <- eligible[keys[eligible] == first_keys[j]]
      i <- rows[1L]
      out[[labels[j]]] <- .mhstats_slice(
        G, rows, labels[j],
        if ("Cohort" %in% meta$strata) dat$Cohort[i] else NA_character_,
        if ("Sex" %in% meta$strata) dat$Sex[i] else NA_integer_,
        all_counts, pseudocount, warn_call_rate, warn_Ae, meta$strata
      )
    }
  }
  class(out) <- c("PaPeR_MHStats", "list")

  if (meta$excluded > 0L)
    warning(meta$excluded, " individual(s) with missing or unknown requested strata ",
            "were excluded from grouped slices but retained in ALL", call. = FALSE)
  summary_rows <- do.call(rbind, unclass(out))
  .mhstats_warning(summary_rows[summary_rows$low_call_rate, , drop = FALSE],
                   "call rate", warn_call_rate)
  .mhstats_warning(summary_rows[summary_rows$low_Ae, , drop = FALSE],
                   "Ae", warn_Ae)
  out
}
