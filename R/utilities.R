# Diagnostic / utility functions (Sequoia analogues), all built on the multiallelic
# relationship-likelihood and OH machinery.

#' Report potential duplicate samples
#'
#' Screens every pair of genotyped individuals for the same-individual (MZ)
#' model against unrelated. Reports evidence for review; never removes or
#' merges samples. Run it before \code{\link{assign_parents}} and remove
#' confirmed duplicates. A reported pair need not favor MZ over every
#' biological relationship.
#'
#' @inheritParams assign_parents
#' @param thresh minimum base-10 MZ-versus-unrelated LLR to report.
#' @param control named list; only the error-model and frequency settings of
#'   \code{\link{pap_control}} are used (\code{error_model}, \code{Err},
#'   \code{frequencies}, \code{n_alleles}, \code{ForeignQ}).
#' @return A data frame ordered by decreasing \code{LLR} (MZ versus U), with
#'   \code{LLR_MZ_vs_PO}, FS, HS, GP, FA, HA and U, \code{LLR_MZ_vs_best},
#'   \code{best_other_relationship}, jointly called and missing locus counts,
#'   genotype mismatches, mismatch rate and opposing homozygotes. All LLRs are
#'   base 10. A mismatch compares unordered genotypes at jointly called loci;
#'   mismatch rate is \code{NA} when none are jointly called.
#' @examples
#' G <- pap_example$GenoM[1:60, ]
#' G <- rbind(G, resampled = G[10, ])   # a second sample of one individual
#' dups <- find_duplicates(G, Dropout = 0.005, Foreign = 5e-4, thresh = 5)
#' dups[, c("ID1", "ID2", "LLR", "LLR_MZ_vs_best", "mismatch_rate")]
#' @export
find_duplicates <- function(GenoM, Dropout = NULL, Foreign = NULL, thresh = 0,
                            control = list()) {
  out <- .find_duplicates(GenoM, .pap_utility_tables(GenoM, Dropout, Foreign, control), thresh)
  .pap_check_removed_path()
  out
}

# The same screen on prebuilt locus tables.
.find_duplicates <- function(GenoM, tables, thresh = 0,
                            Complex = c("simp", "full")) {
  Complex <- match.arg(Complex)
  G <- .prepare_for_tables(GenoM, tables)
  e <- rel_edges(G, tables, "MZ", thresh = thresh, oh_max = -1L,
                 Complex = Complex, operation = "duplicate")
  ctx <- .new_inference_context(G, tables, Complex, include_duplicate = TRUE)
  raw <- if (nrow(e)) .score_unconditional_pairs(ctx, e$id1, e$id2) else
    matrix(numeric(), 0L, nrow(ctx$Delta),
           dimnames = list(NULL, rownames(ctx$Delta)))
  other <- c("PO", "FS", "HS", "GP", "FA", "HA", "U")
  contrasts <- lapply(other, function(rel) .llr_log10(
    .log_likelihood_contrast_ln(raw[, "MZ"], raw[, rel], "duplicate report")))
  names(contrasts) <- paste0("LLR_MZ_vs_", other)
  best <- if (nrow(raw)) apply(raw[, other, drop = FALSE], 1L, max) else numeric()
  labels <- if (nrow(raw)) vapply(seq_len(nrow(raw)), function(i)
    paste(other[raw[i, other] == best[i]], collapse = "|"), character(1)) else character()
  n <- nrow(e)
  d <- data.frame(ID1 = e$id1, ID2 = e$id2, LLR = e$llr,
    winning_configuration = e$focal_winner,
    winning_configuration_set = e$focal_winner_set,
    alternative_winner = e$alternative_winner,
    alternative_winner_set = e$alternative_winner_set,
    duplicate_guard_llr = e$llr,
    duplicate_guard_status = rep("candidate", n),
    best_other_relationship = labels,
    LLR_MZ_vs_best = .llr_log10(.log_likelihood_contrast_ln(
      raw[, "MZ"], best, "duplicate report best alternative")),
    stringsAsFactors = FALSE)
  for (name in names(contrasts)) d[[name]] <- contrasts[[name]]
  for (name in c("n_loci", "jointly_called", "missing_id1", "missing_id2",
                 "genotype_mismatches", "opposing_homozygotes")) d[[name]] <- integer(n)
  d$mismatch_rate <- rep(NA_real_, n)
  first <- seq.int(1L, ncol(G), by = 2L); second <- first + 1L
  for (i in seq_len(n)) {
    a <- G[match(d$ID1[i], rownames(G)), ]; b <- G[match(d$ID2[i], rownames(G)), ]
    called_a <- a[first] > 0L; called_b <- b[first] > 0L
    joint <- called_a & called_b
    mismatch <- joint & (a[first] != b[first] | a[second] != b[second])
    d$n_loci[i] <- length(first); d$jointly_called[i] <- sum(joint)
    d$missing_id1[i] <- sum(!called_a); d$missing_id2[i] <- sum(!called_b)
    d$genotype_mismatches[i] <- sum(mismatch)
    d$opposing_homozygotes[i] <- sum(joint & a[first] == a[second] &
      b[first] == b[second] & a[first] != b[first])
    if (any(joint)) d$mismatch_rate[i] <- sum(mismatch) / sum(joint)
  }
  d <- d[order(-d$LLR, d$ID1, d$ID2), , drop = FALSE]
  rownames(d) <- NULL
  attr(d, "relationship_universe") <- attr(e, "relationship_universe")
  attr(d, "log_base") <- .public_log_base
  attr(d, "threshold_log10") <- thresh
  d
}

#' Pairs already related (parent-offspring or share a parent) in a pedigree
#' @keywords internal
#' @noRd
.known_related_pairs <- function(ped) {
  keyp <- function(a, b) paste(pmin(a, b), pmax(a, b))
  k <- character(0)
  po <- !is.na(ped$dam); k <- c(k, keyp(ped$id[po], ped$dam[po]))
  po <- !is.na(ped$sire); k <- c(k, keyp(ped$id[po], ped$sire[po]))
  for (par in c("dam", "sire")) {          # sibs: share a parent
    grp <- split(ped$id, ped[[par]])
    for (g in grp) if (length(g) > 1) {
      cb <- utils::combn(g, 2); k <- c(k, keyp(cb[1, ], cb[2, ]))
    }
  }
  unique(k)
}

#' Find likely relatives that are not in a pedigree
#'
#' Reports pairs whose most likely relationship is a close relative (base-10
#' LLR against unrelated above \code{thresh}) and that are not already
#' parent-offspring or siblings in \code{Pedigree}. Useful for spotting missed
#' relationships (the analogue of Sequoia's \code{GetMaybeRel}).
#'
#' @inheritParams find_duplicates
#' @param Pedigree optional \code{data.frame(id, dam, sire)}; pairs already
#'   related in it are left out.
#' @param LifeHistData unused; kept for call compatibility. All pairs are
#'   scanned.
#' @param thresh minimum base-10 LLR of the top relationship against unrelated.
#' @param exclude_rels relationship labels to drop (none by default).
#' @return \code{data.frame(ID1, ID2, TopRel, LLR, ...)}, strongest first.
#'   \code{LLR} is base 10.
#' @examples
#' G <- pap_example$GenoM[1:60, ]
#' rel <- get_maybe_rel(G, Dropout = 0.005, Foreign = 5e-4, thresh = 5)
#' head(rel[, c("ID1", "ID2", "TopRel", "LLR")])
#' @export
get_maybe_rel <- function(GenoM, Dropout = NULL, Foreign = NULL, Pedigree = NULL,
                          LifeHistData = NULL, thresh = 2, exclude_rels = NULL,
                          control = list()) {
  out <- .get_maybe_rel(GenoM, .pap_utility_tables(GenoM, Dropout, Foreign, control),
                        Pedigree = Pedigree, thresh = thresh, exclude_rels = exclude_rels)
  .pap_check_removed_path()
  out
}

# The same scan on prebuilt locus tables.
.get_maybe_rel <- function(GenoM, tables, Pedigree = NULL, LifeHistData = NULL,
                          thresh = 2, exclude_rels = NULL,
                          Complex = c("simp", "full")) {
  Complex <- match.arg(Complex)
  excluded <- .resolve_relationship_names(
    exclude_rels, Complex, "excluded relationship")
  e <- rel_edges(GenoM, tables, "RELATED", alts = "U",
                 thresh = thresh, oh_max = -1L, Complex = Complex,
                 operation = "diagnostic")
  if (!nrow(e)) {
    res <- data.frame(
      ID1 = character(0), ID2 = character(0), TopRel = character(0),
      TopRelSet = character(0), LLR = numeric(0),
      alternative_winner = character(0), alternative_winner_set = character(0),
      duplicate_guard_llr = numeric(0), duplicate_guard_status = character(0),
      stringsAsFactors = FALSE)
    attr(res, "relationship_universe") <- attr(e, "relationship_universe")
    attr(res, "log_base") <- .public_log_base
    return(res)
  }
  res <- data.frame(ID1 = e$id1, ID2 = e$id2,
                    TopRel = e$focal_winner, TopRelSet = e$focal_winner_set,
                    LLR = e$llr,
                    alternative_winner = e$alternative_winner,
                    alternative_winner_set = e$alternative_winner_set,
                    duplicate_guard_llr = e$duplicate_guard_llr,
                    duplicate_guard_status = e$duplicate_guard_status,
                    stringsAsFactors = FALSE)
  res <- res[!res$TopRel %in% excluded, , drop = FALSE]
  if (!is.null(Pedigree)) {
    known <- .known_related_pairs(Pedigree)
    key <- paste(pmin(res$ID1, res$ID2), pmax(res$ID1, res$ID2))
    res <- res[!key %in% known, , drop = FALSE]
  }
  res <- res[order(-res$LLR), ]
  attr(res, "relationship_universe") <- attr(e, "relationship_universe")
  attr(res, "log_base") <- .public_log_base
  res
}
