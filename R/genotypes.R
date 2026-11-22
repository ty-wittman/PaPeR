#' Prepare public microhaplotype genotypes for PaPeR
#'
#' Public PaPeR genotypes have one character column per locus. Called
#' genotypes are unordered allele-identity pairs such as \code{"1/3"}; \code{"0"}
#' is a completely missing locus. This function validates that contract and
#' converts it once to the private two-copy integer matrix consumed by Rcpp.
#'
#' @param GenoM character matrix or data frame with unique, nonempty individual
#'   row names and unique, nonempty locus column names.
#' @param n_alleles optional positive integer vector with one value per locus.
#'   When omitted, each value is inferred as the largest observed allele code.
#' @return A private, class-marked integer matrix. Its attributes retain the
#'   validated row IDs, locus names, declared allele counts, and two-copy column
#'   map. Public functions accept this result so nested package calls can reuse
#'   it without parsing again.
#' @keywords internal
#' @noRd
prepare_genotypes <- function(GenoM, n_alleles = NULL) {
  if (inherits(GenoM, "PaPeR_prepared_genotypes")) {
    .validate_prepared_genotypes(GenoM, n_alleles)
    return(GenoM)
  }
  if (is.matrix(GenoM) && (is.integer(GenoM) || is.numeric(GenoM)))
    stop("unmarked two-column numeric/integer GenoM is a legacy internal format; ",
         "supply one character 'allele1/allele2' column per locus")
  if (!is.matrix(GenoM) && !is.data.frame(GenoM))
    stop("GenoM must be a character matrix or data frame")

  if (is.data.frame(GenoM) && .row_names_info(GenoM, type = 1L) < 0L)
    stop("GenoM must have explicit individual row names")
  ids <- rownames(GenoM)
  loci <- colnames(GenoM)
  .validate_genotype_names(ids, "individual row names")
  .validate_genotype_names(loci, "locus column names")
  if (nrow(GenoM) < 1L) stop("GenoM must contain at least one individual")
  if (ncol(GenoM) < 1L) stop("GenoM must contain at least one locus")

  if (is.data.frame(GenoM)) {
    ok <- vapply(GenoM, function(x) is.character(x) || is.factor(x), logical(1))
    if (!all(ok))
      stop("every GenoM data-frame column must be character or factor")
    x <- as.matrix(data.frame(lapply(GenoM, as.character),
                              check.names = FALSE, stringsAsFactors = FALSE))
    dimnames(x) <- list(ids, loci)
  } else {
    if (!is.character(GenoM)) stop("GenoM matrix must be character")
    x <- GenoM
  }
  if (anyNA(x)) stop("GenoM must not contain NA genotypes")

  k_supplied <- .validate_n_alleles(n_alleles, ncol(x), loci)
  out <- matrix(0L, nrow(x), 2L * ncol(x), dimnames = list(ids, NULL))
  inferred <- integer(ncol(x))
  for (l in seq_len(ncol(x))) {
    z <- x[, l]
    if (any(z == ""))
      stop("blank genotype at locus '", loci[l], "'")
    called <- z != "0"
    if (!any(called)) {
      if (is.null(k_supplied))
        stop("cannot infer n_alleles for all-missing locus '", loci[l],
             "'; supply n_alleles")
      inferred[l] <- k_supplied[l]
      next
    }
    pat <- "^(0|[1-9][0-9]*)/(0|[1-9][0-9]*)$"
    bad <- called & !grepl(pat, z, perl = TRUE)
    if (any(bad))
      stop("invalid genotype '", z[which(bad)[1L]], "' at individual '",
           ids[which(bad)[1L]], "', locus '", loci[l],
           "': use exactly 'a/b' with canonical nonnegative integers, or '0'")
    pieces <- strsplit(z[called], "/", fixed = TRUE)
    av <- vapply(pieces, `[[`, character(1), 1L)
    bv <- vapply(pieces, `[[`, character(1), 2L)
    an <- suppressWarnings(as.numeric(av))
    bn <- suppressWarnings(as.numeric(bv))
    if (any(!is.finite(an) | !is.finite(bn) |
            an > .Machine$integer.max | bn > .Machine$integer.max))
      stop("allele code exceeds the supported integer range at locus '", loci[l], "'")
    miss <- an == 0 | bn == 0
    if (!is.null(k_supplied) && any(c(an, bn) > k_supplied[l]))
      stop("allele code exceeds n_alleles at locus '", loci[l], "'")
    observed_max <- max(c(an[!miss], bn[!miss], 0))
    if (is.null(k_supplied) && observed_max == 0)
      stop("cannot infer n_alleles for all-missing locus '", loci[l],
           "'; supply n_alleles")
    inferred[l] <- if (is.null(k_supplied)) as.integer(observed_max) else k_supplied[l]
    rows <- which(called)[!miss]
    if (length(rows)) {
      aa <- as.integer(an[!miss])
      bb <- as.integer(bn[!miss])
      out[rows, 2L * l - 1L] <- pmin(aa, bb)
      out[rows, 2L * l] <- pmax(aa, bb)
    }
  }
  .mark_prepared_genotypes(out, loci, inferred)
}

.validate_genotype_names <- function(x, what) {
  if (is.null(x) || anyNA(x) || any(x == "") || anyDuplicated(x))
    stop("GenoM must have unique, nonempty ", what)
}

.validate_n_alleles <- function(n_alleles, n_loci, loci) {
  if (is.null(n_alleles)) return(NULL)
  if (!is.numeric(n_alleles) || length(n_alleles) != n_loci ||
      anyNA(n_alleles) || any(!is.finite(n_alleles)) ||
      any(n_alleles < 1 | n_alleles != floor(n_alleles)) ||
      any(n_alleles > .Machine$integer.max))
    stop("n_alleles must contain one positive integer per locus")
  if (!is.null(names(n_alleles)) &&
      !identical(as.character(names(n_alleles)), as.character(loci)))
    stop("named n_alleles must match GenoM locus names and ordering")
  as.integer(n_alleles)
}

.mark_prepared_genotypes <- function(x, loci, n_alleles) {
  storage.mode(x) <- "integer"
  colnames(x) <- as.vector(rbind(paste0(loci, "::copy1"),
                                paste0(loci, "::copy2")))
  attr(x, "row_ids") <- rownames(x)
  attr(x, "locus_names") <- loci
  attr(x, "n_alleles") <- as.integer(n_alleles)
  attr(x, "locus_columns") <- matrix(seq_len(ncol(x)), ncol = 2L, byrow = TRUE,
                                     dimnames = list(loci, c("copy1", "copy2")))
  class(x) <- c("PaPeR_prepared_genotypes", "matrix", "array")
  x
}

.validate_prepared_genotypes <- function(GenoM, n_alleles = NULL) {
  if (!is.matrix(GenoM) || !is.integer(GenoM))
    stop("invalid prepared GenoM object")
  ids <- attr(GenoM, "row_ids", exact = TRUE)
  loci <- attr(GenoM, "locus_names", exact = TRUE)
  k <- attr(GenoM, "n_alleles", exact = TRUE)
  mapping <- attr(GenoM, "locus_columns", exact = TRUE)
  if (!is.character(ids) || !identical(ids, rownames(GenoM)) ||
      !is.character(loci) || ncol(GenoM) != 2L * length(loci) ||
      !is.integer(k) || length(k) != length(loci) ||
      !is.matrix(mapping) || !is.integer(mapping) ||
      !identical(dim(mapping), c(length(loci), 2L)))
    stop("invalid prepared GenoM object")
  .validate_genotype_names(ids, "individual row names")
  .validate_genotype_names(loci, "locus names")
  expected_mapping <- matrix(seq_len(ncol(GenoM)), ncol = 2L, byrow = TRUE,
                             dimnames = list(loci, c("copy1", "copy2")))
  if (!identical(mapping, expected_mapping))
    stop("invalid prepared GenoM copy-column mapping")
  if (!identical(.validate_n_alleles(k, length(loci), loci), k))
    stop("invalid prepared GenoM allele-count metadata")
  if (any(GenoM < 0L) || anyNA(GenoM))
    stop("invalid allele codes in prepared GenoM")
  for (l in seq_along(k)) {
    a <- GenoM[, 2L * l - 1L]
    b <- GenoM[, 2L * l]
    if (any(xor(a == 0L, b == 0L)) || any(a > b) || any(a > k[l] | b > k[l]))
      stop("invalid prepared genotypes at locus '", loci[l], "'")
  }
  if (!is.null(n_alleles)) {
    supplied <- .validate_n_alleles(n_alleles, length(loci), loci)
    if (!identical(supplied, as.integer(k)))
      stop("prepared GenoM n_alleles do not match the supplied values")
  }
  invisible(TRUE)
}

.table_n_alleles <- function(tables) {
  if (!is.list(tables) || !length(tables))
    stop("tables must be a nonempty list from build_tables")
  k <- vapply(tables, function(x) {
    if (!is.list(x) || !is.numeric(x$k) || length(x$k) != 1L ||
        is.na(x$k) || !is.finite(x$k) || x$k < 1L || x$k != floor(x$k))
      return(NA_integer_)
    as.integer(x$k)
  }, integer(1))
  if (anyNA(k)) stop("tables must contain valid per-locus allele counts")
  k
}

.table_locus_names <- function(tables) {
  loci <- names(tables)
  if (is.null(loci)) return(NULL)
  if (!is.character(loci) || length(loci) != length(tables) || anyNA(loci) ||
      any(!nzchar(loci)) || anyDuplicated(loci))
    stop("named tables must have unique, non-missing locus names")
  loci
}

.prepare_for_tables <- function(GenoM, tables) {
  .require_positive_error_rates(tables)
  G <- prepare_genotypes(GenoM, .table_n_alleles(tables))
  table_loci <- .table_locus_names(tables)
  if (!is.null(table_loci) &&
      !identical(table_loci, attr(G, "locus_names", exact = TRUE)))
    stop("named tables must match GenoM locus names and order", call. = FALSE)
  G
}
