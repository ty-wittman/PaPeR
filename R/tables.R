#' Build the per-locus likelihood tables
#'
#' Precomputes compact per-locus data for the pairwise and trio likelihoods:
#' HWE prior, genotype-index metadata, an implicit observed-given-true error
#' kernel, and the unrelated observed marginal. The returned production objects
#' never retain a dense genotype-by-genotype relationship, Jacquard, emission,
#' or error-convolved matrix.
#'
#' @param freqs one explicitly selected group slice returned by \code{MHStats},
#'   or a list of likelihood-safe per-locus allele-frequency vectors. A named
#'   panel is matched to genotype locus names in the same order at inference.
#' @param Err legacy called-genotype discordance probability: one panel-wide scalar.
#'   Default \code{1e-4}, matching Sequoia. Valid
#'   only for \code{ErrorModel = "legacy_sequoia"}; must be strictly positive.
#' @param ErrProfile optional legacy panel-wide transition profile; see
#'   \code{expand_errors}. Valid only for \code{ErrorModel =
#'   "legacy_sequoia"}.
#' @param ErrorModel explicit observation model. \code{"legacy_sequoia"}
#'   preserves the existing \code{Err}/\code{ErrProfile} model and remains the
#'   default. \code{"dropout_foreign"} is opt-in.
#' @param Dropout for \code{"dropout_foreign"}, the conditional probability
#'   that a called true heterozygote collapses to one of its constituent
#'   homozygotes. Supply one strictly positive panel-wide scalar.
#' @param Foreign for \code{"dropout_foreign"}, the per-copy probability of a
#'   foreign-allele draw. A draw may reproduce the precursor allele, so this is
#'   not the probability that a copy necessarily changes. Supply one strictly
#'   positive panel-wide scalar.
#' @param ForeignQ optional foreign-allele distribution. \code{NULL} uses the
#'   locus population frequencies. Otherwise supply one vector for one locus or
#'   an exact per-locus list, with each vector normalized and aligned to the
#'   declared allele universe.
#' @details \code{Dropout} and \code{Foreign} are required explicitly for the
#'   opt-in model; this interface does not assert universal biological defaults.
#'   \code{Dropout} is conditional called-genotype collapse, not raw per-copy
#'   PCR dropout. \code{Foreign} is a hard-call foreign-draw opportunity, not a
#'   contaminating-read fraction. With \code{ForeignQ = NULL}, \code{q = p}
#'   represents the locuswise marginal allele of a random same-population
#'   contaminant. Missing observations remain separate all-ones emissions.
#' @return list of per-locus tables, one per locus, to pass to the LL wrappers.
#' @keywords internal
#' @noRd
build_tables <- function(freqs, Err = 1e-4, ErrProfile = NULL,
                         ErrorModel = c("legacy_sequoia", "dropout_foreign"),
                         Dropout = NULL, Foreign = NULL, ForeignQ = NULL) {
  ErrorModel <- match.arg(ErrorModel)
  if (inherits(freqs, "PaPeR_MHStats"))
    stop("select one MHStats group slice, for example stats[[\"ALL\"]]",
         call. = FALSE)
  if (inherits(freqs, "PaPeR_MHStats_slice")) {
    if (!all(c("locus", "n_alleles_universe", "likelihood_frequencies") %in%
             names(freqs)))
      stop("invalid MHStats group slice", call. = FALSE)
    loci <- freqs$locus
    if (!is.character(loci) || length(loci) != nrow(freqs) || anyNA(loci) ||
        any(!nzchar(loci)) || anyDuplicated(loci))
      stop("MHStats group slice must have unique, non-missing locus names",
           call. = FALSE)
    n_alleles <- freqs$n_alleles_universe
    if (!is.numeric(n_alleles) || length(n_alleles) != nrow(freqs) ||
        anyNA(n_alleles) || any(!is.finite(n_alleles)) ||
        any(n_alleles < 1L | n_alleles != floor(n_alleles)))
      stop("MHStats group slice has invalid declared allele universes",
           call. = FALSE)
    freqs <- freqs$likelihood_frequencies
    if (!is.list(freqs) || length(freqs) != length(loci))
      stop("invalid MHStats likelihood_frequencies", call. = FALSE)
    freqs <- Map(.validate_likelihood_frequencies, freqs, loci, n_alleles)
    names(freqs) <- loci
  } else {
    if (!is.list(freqs) || !length(freqs))
      stop("freqs must be a nonempty list of likelihood frequency vectors",
           call. = FALSE)
    loci <- names(freqs)
    if (!is.null(loci)) {
      if (!is.character(loci) || anyNA(loci) || any(!nzchar(loci)) ||
          anyDuplicated(loci))
        stop("named frequency panels must have unique, non-missing locus names",
             call. = FALSE)
    }
    labels <- if (is.null(loci)) rep("<unnamed>", length(freqs)) else loci
    freqs <- lapply(seq_along(freqs), function(i)
      .validate_likelihood_frequencies(freqs[[i]], labels[i], n_alleles = NULL))
    if (!is.null(loci)) names(freqs) <- loci
  }
  n_loci <- length(freqs)
  if (identical(ErrorModel, "legacy_sequoia")) {
    if (!missing(Dropout) || !missing(Foreign) || !missing(ForeignQ))
      stop("Dropout, Foreign, and ForeignQ require ErrorModel='dropout_foreign'",
           call. = FALSE)
    e <- if (missing(Err)) expand_errors(n_loci = n_loci,
                                         ErrProfile = ErrProfile) else
      expand_errors(Err, n_loci, ErrProfile)
    ans <- lapply(seq_len(n_loci), function(l) {
      profile <- if (e$profile_default) NULL else e$ErrProfile
      build_locus_table(freqs[[l]], e$Err[l], profile)
    })
  } else {
    if (!missing(Err) || !missing(ErrProfile))
      stop("Err and ErrProfile require ErrorModel='legacy_sequoia'",
           call. = FALSE)
    if (is.null(Dropout) || is.null(Foreign))
      stop("Dropout and Foreign must be supplied for ErrorModel='dropout_foreign'",
           call. = FALSE)
    spec <- .validate_dropout_foreign_spec(Dropout, Foreign, ForeignQ,
                                           freqs, n_loci)
    ans <- lapply(seq_len(n_loci), function(l) {
      out <- build_locus_table(freqs[[l]], ErrorModel = ErrorModel,
                               Dropout = spec$Dropout[l],
                               Foreign = spec$Foreign[l],
                               ForeignQ = spec$ForeignQ[[l]])
      out$ForeignQMode <- spec$ForeignQMode[l]
      out
    })
  }
  names(ans) <- names(freqs)
  .require_positive_error_rates(ans)
  ans
}

.validate_likelihood_frequencies <- function(p, locus, n_alleles = NULL) {
  if (!is.numeric(p) || is.object(p) || !length(p) || anyNA(p) ||
      any(!is.finite(p)) || any(p <= 0))
    stop("unavailable or zero likelihood frequencies for locus '", locus,
         "'; supply positive pseudocount and called genotypes", call. = FALSE)
  if (!is.null(n_alleles) && length(p) != n_alleles)
    stop("likelihood frequency length for locus '", locus,
         "' does not match its declared allele universe", call. = FALSE)
  allele_names <- names(p)
  if (!is.null(allele_names) &&
      !identical(as.character(allele_names), as.character(seq_along(p))))
    stop("likelihood frequency allele names for locus '", locus,
         "' must be the integer allele universe in order", call. = FALSE)
  total <- sum(p)
  if (!is.finite(total) || total <= 0)
    stop("likelihood frequencies for locus '", locus,
         "' must have a positive finite total", call. = FALSE)
  p <- as.numeric(p) / total
  if (!is.null(allele_names)) names(p) <- allele_names
  p
}

#' Build one internal microhaplotype likelihood table
#'
#' @param p allele-frequency vector.
#' @param Err called-genotype error probability. Default \code{1e-4}.
#' @param ErrProfile optional three-weight profile or symmetric 3 x 3 adapter.
#' @return A locus likelihood table.
#' @keywords internal
#' @noRd
build_locus_table <- function(p, Err = 1e-4, ErrProfile = NULL,
                              ErrorModel = c("legacy_sequoia", "dropout_foreign"),
                              Dropout = NULL, Foreign = NULL, ForeignQ = NULL) {
  ErrorModel <- match.arg(ErrorModel)
  if (!is.numeric(p) || length(p) < 1L || anyNA(p) || any(!is.finite(p)) ||
      any(p < 0) || sum(p) <= 0)
    stop("p must be a finite, nonnegative allele-frequency vector with positive sum")
  if (identical(ErrorModel, "dropout_foreign")) {
    if (!missing(Err) || !missing(ErrProfile))
      stop("Err and ErrProfile require ErrorModel='legacy_sequoia'",
           call. = FALSE)
    if (is.null(Dropout) || is.null(Foreign))
      stop("Dropout and Foreign must be supplied for ErrorModel='dropout_foreign'",
           call. = FALSE)
    p_norm <- as.numeric(p) / sum(p)
    if (!is.null(names(p))) names(p_norm) <- names(p)
    spec <- .validate_dropout_foreign_spec(Dropout, Foreign, ForeignQ,
                                           list(p_norm), 1L)
    q <- spec$ForeignQ[[1L]]
    ans <- build_locus_table_dropout_foreign_cpp(
      p_norm, spec$Dropout[1L], spec$Foreign[1L], unname(q))
    ans$ForeignQ <- q
    ans$ForeignQMode <- spec$ForeignQMode[1L]
    return(ans)
  }
  if (!missing(Dropout) || !missing(Foreign) || !missing(ForeignQ))
    stop("Dropout, Foreign, and ForeignQ require ErrorModel='dropout_foreign'",
         call. = FALSE)
  e <- if (missing(Err)) expand_errors(n_loci = 1L, ErrProfile = ErrProfile) else
    expand_errors(Err, 1L, ErrProfile)
  if (length(p) == 1L && e$Err[1L] > 0)
    stop("a monomorphic locus cannot represent positive Err")
  w <- if (e$profile_default) {
    E <- e$Err[1L]
    c(w_Het = E - (E / 2)^2, w_Hom = E / 2, w_Far = (E / 2)^2)
  } else e$ErrProfile
  k <- length(p)
  hom_raw <- (k - 1) * w["w_Het"] + (k - 1) * k / 2 * w["w_Far"]
  het_raw <- 2 * w["w_Hom"] + 2 * (k - 2) * w["w_Het"] +
    (k - 2) * (k - 1) / 2 * w["w_Far"]
  if (e$Err[1L] > 0 && (hom_raw <= 0 || (k > 1L && het_raw <= 0)))
    stop("ErrProfile is non-normalizable for at least one true genotype row")
  ans <- build_locus_table_cpp(p, e$Err[1L], unname(w))
  ans$Err <- e$Err[1L]
  ans$ErrProfile <- w
  ans$ErrProfileDefault <- e$profile_default
  ans
}

#' Resolve individual IDs (character) or 1-based rows to 0-based GenoM rows
#' @keywords internal
#' @noRd
.row0 <- function(GenoM, id) {
  if (is.character(id)) {
    r <- match(id, rownames(GenoM))
    if (anyNA(r)) stop("unknown individual id(s): ",
                       paste(id[is.na(r)], collapse = ", "))
    return(r - 1L)
  }
  as.integer(id) - 1L
}

#' Relationship log-likelihoods for one pair of individuals
#'
#' @inheritParams find_duplicates
#' @param i,j individual IDs (character) or 1-based row indices of \code{GenoM}.
#' @return Named natural-log likelihoods of the pair: \code{U} (unrelated),
#'   \code{PO}, \code{FS}, \code{DEG2} (second degree: half siblings,
#'   grandparent or full avuncular, which a pair's genotypes cannot tell
#'   apart), \code{HA} (half avuncular) and \code{MZ} (the same individual).
#'   Missing calls are marginalized; fully missing loci contribute zero. These
#'   are raw log-likelihoods, not LLRs (attributes \code{log_base = "e"} and
#'   \code{value_type = "raw_log_likelihood"}). The attribute
#'   \code{configuration_likelihoods} gives the same values labelled by
#'   relationship (U, PO, FS, HS, GP, FA, HA), and
#'   \code{winning_configuration} names the most likely one.
#' @examples
#' ex <- pap_example
#' kid <- "ID0151"
#' mum <- ex$Pedigree$dam[ex$Pedigree$id == kid]
#' ll <- calc_pair_LL(ex$GenoM, kid, mum, Dropout = 0.005, Foreign = 5e-4)
#' attr(ll, "winning_configuration")
#' round(c(ll) - ll[["U"]], 1)   # natural-log ratios against unrelated
#' @export
calc_pair_LL <- function(GenoM, i, j, Dropout = NULL, Foreign = NULL, control = list()) {
  out <- .calc_pair_LL(GenoM, i, j, .pap_utility_tables(GenoM, Dropout, Foreign, control))
  .pap_check_removed_path()
  out
}

# The same likelihoods on prebuilt locus tables.
.calc_pair_LL <- function(GenoM, i, j, tables, Complex = c("simp", "full")) {
  Complex <- match.arg(Complex)
  universe <- .relationship_universe(Complex, "RELATED", operation = "pair")
  ctx <- .new_inference_context(GenoM, tables, Complex)
  id <- rownames(ctx$G)
  ii <- if (is.character(i)) i else id[as.integer(i)]
  jj <- if (is.character(j)) j else id[as.integer(j)]
  score <- .score_unconditional_pairs(ctx, ii, jj)[1L, ]
  legacy <- calc_pair_LL_cpp(ctx$markers, ctx$G, .row0(ctx$G, i), .row0(ctx$G, j))
  attr(legacy, "Complex") <- Complex
  attr(legacy, "configuration_likelihoods") <- score
  attr(legacy, "winning_configuration") <- names(score)[which.max(score)]
  attr(legacy, "relationship_registry") <- ctx$registry
  attr(legacy, "relationship_universe") <- universe
  attr(legacy, "log_base") <- "e"
  attr(legacy, "value_type") <- "raw_log_likelihood"
  legacy
}

#' Parent-pair (trio) log-likelihoods
#'
#' @inheritParams find_duplicates
#' @param o,d,s offspring, candidate dam and candidate sire (IDs or 1-based
#'   rows of \code{GenoM}).
#' @return Named natural-log likelihoods: H_PP (both are parents), H_D0 (the
#'   dam only), H_0S (the sire only) and H_U (unrelated). Missing calls are
#'   marginalized exactly. These are raw log-likelihoods, not LLRs.
#' @examples
#' ex <- pap_example
#' kid <- "ID0151"
#' mum <- ex$Pedigree$dam[ex$Pedigree$id == kid]
#' dad <- ex$Pedigree$sire[ex$Pedigree$id == kid]
#' ll <- calc_trio_LL(ex$GenoM, kid, mum, dad, Dropout = 0.005, Foreign = 5e-4)
#' round(c(ll) - ll[["H_U"]], 1)   # natural-log ratios against no parents
#' @export
calc_trio_LL <- function(GenoM, o, d, s, Dropout = NULL, Foreign = NULL, control = list()) {
  out <- .calc_trio_LL(GenoM, o, d, s, .pap_utility_tables(GenoM, Dropout, Foreign, control))
  .pap_check_removed_path()
  out
}

# The same likelihoods on prebuilt locus tables.
.calc_trio_LL <- function(GenoM, o, d, s, tables) {
  G <- .prepare_for_tables(GenoM, tables)
  out <- calc_trio_LL_cpp(tables, G, .row0(G, o), .row0(G, d), .row0(G, s))
  attr(out, "log_base") <- "e"
  attr(out, "value_type") <- "raw_log_likelihood"
  out
}
