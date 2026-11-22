.compact_locus_schema <- "compact-locus-v1"
.compact_locus_schema_v2 <- "compact-locus-v2"

.validate_compact_locus_tables <- function(tables) {
  if (!is.list(tables) || !length(tables))
    stop("tables must be a nonempty compact locus-table list", call. = FALSE)
  obsolete <- c("S", "YS", "X", "Delta", "relnames")
  required_common <- c("schema_version", "schema_id", "k", "G", "geno", "p", "pi",
                       "mObs", "p_oh_true_po", "p_oh_unrelated",
                       "p_me_true_trio")
  for (l in seq_along(tables)) {
    tb <- tables[[l]]
    schema <- if (is.list(tb)) tb$schema_version else NULL
    if (!is.character(schema) || length(schema) != 1L || is.na(schema) ||
        !schema %in% c(.compact_locus_schema, .compact_locus_schema_v2))
      stop("locus table ", l, " must use ", .compact_locus_schema,
           "; rebuild with build_tables()", call. = FALSE)
    bad <- intersect(obsolete, names(tb))
    if (length(bad))
      stop("locus table ", l, " retains obsolete dense field(s): ",
           paste(bad, collapse = ", "), "; rebuild with build_tables()", call. = FALSE)
    required_model <- if (identical(schema, .compact_locus_schema))
      c("Err", "ErrProfile", "Z_hom", "Z_het") else
      c("ErrorModel", "Dropout", "Foreign", "ForeignQ", "ForeignQMode")
    missing <- setdiff(c(required_common, required_model), names(tb))
    if (length(missing))
      stop("compact locus table ", l, " is missing field(s): ",
           paste(missing, collapse = ", "), call. = FALSE)
    if (!is.numeric(tb$schema_id) || length(tb$schema_id) != 1L ||
        tb$schema_id != if (identical(schema, .compact_locus_schema)) 1L else 2L)
      stop("compact locus table ", l, " has inconsistent schema_id", call. = FALSE)
    if (!is.numeric(tb$k) || length(tb$k) != 1L || is.na(tb$k) || tb$k < 1L ||
        tb$k != as.integer(tb$k) || !is.numeric(tb$G) || length(tb$G) != 1L ||
        tb$G != tb$k * (tb$k + 1L) / 2L || !is.integer(tb$geno) ||
        !identical(dim(tb$geno), c(as.integer(tb$G), 2L)) ||
        !is.numeric(tb$p) || length(tb$p) != tb$k || !is.numeric(tb$pi) ||
        length(tb$pi) != tb$G || !is.numeric(tb$mObs) || length(tb$mObs) != tb$G ||
        !is.numeric(tb$p_oh_true_po) || length(tb$p_oh_true_po) != 1L ||
        !is.numeric(tb$p_oh_unrelated) || length(tb$p_oh_unrelated) != 1L ||
        !is.numeric(tb$p_me_true_trio) || length(tb$p_me_true_trio) != 1L ||
        !is.finite(tb$p_oh_true_po) || !is.finite(tb$p_oh_unrelated) ||
        !is.finite(tb$p_me_true_trio) ||
        tb$p_oh_true_po < 0 || tb$p_oh_true_po > 1 ||
        tb$p_oh_unrelated < 0 || tb$p_oh_unrelated > 1 ||
        tb$p_me_true_trio < 0 || tb$p_me_true_trio > 1)
      stop("compact locus table ", l, " has malformed dimensions", call. = FALSE)
    if (identical(schema, .compact_locus_schema)) {
      if (any(c("ErrorModel", "Dropout", "Foreign", "ForeignQ",
                "ForeignQMode") %in% names(tb)))
        stop("legacy compact locus table ", l,
             " contains new-model fields", call. = FALSE)
    } else {
      if (!identical(tb$ErrorModel, "dropout_foreign") ||
          !is.numeric(tb$Dropout) || length(tb$Dropout) != 1L ||
          !is.finite(tb$Dropout) || tb$Dropout < 0 || tb$Dropout > 1 ||
          !is.numeric(tb$Foreign) || length(tb$Foreign) != 1L ||
          !is.finite(tb$Foreign) || tb$Foreign < 0 || tb$Foreign > 1 ||
          !is.numeric(tb$ForeignQ) || length(tb$ForeignQ) != tb$k ||
          anyNA(tb$ForeignQ) || any(!is.finite(tb$ForeignQ)) ||
          any(tb$ForeignQ < 0) || sum(tb$ForeignQ) <= 0 ||
          abs(sum(tb$ForeignQ) - 1) > 1e-12 ||
          !is.character(tb$ForeignQMode) || length(tb$ForeignQMode) != 1L ||
          !tb$ForeignQMode %in% c("p", "supplied"))
        stop("compact locus table ", l, " has malformed dropout_foreign parameters",
             call. = FALSE)
      if (!is.null(names(tb$ForeignQ))) {
        expected <- if (is.null(names(tb$p))) as.character(seq_len(tb$k)) else
          names(tb$p)
        if (!identical(as.character(names(tb$ForeignQ)), as.character(expected)))
          stop("compact locus table ", l,
               " ForeignQ names do not match the allele universe", call. = FALSE)
      }
      if (any(c("Err", "ErrProfile", "Z_hom", "Z_het") %in% names(tb)))
        stop("dropout_foreign compact locus table ", l,
             " contains legacy fields", call. = FALSE)
    }
    numeric_matrices <- Filter(function(x) is.matrix(x) && is.numeric(x), tb)
    if (any(vapply(numeric_matrices, function(x)
      identical(dim(x), rep(as.integer(tb$G), 2L)), logical(1))))
      stop("compact locus table ", l, " retains a forbidden numeric G x G matrix",
           call. = FALSE)
  }
  invisible(tables)
}

.prepare_observed_indices <- function(G, tables) {
  if (!is.matrix(G) || ncol(G) != 2L * length(tables))
    stop("prepared genotypes and tables have incompatible dimensions",
         call. = FALSE)
  out <- matrix(-1L, nrow(G), length(tables))
  for (l in seq_along(tables)) {
    k <- tables[[l]]$k
    a <- G[, 2L * l - 1L]
    b <- G[, 2L * l]
    typed <- a > 0L & b > 0L
    lo <- pmin(a[typed], b[typed])
    hi <- pmax(a[typed], b[typed])
    out[typed, l] <- as.integer(
      (lo - 1L) * k - (lo - 1L) * (lo - 2L) / 2L + (hi - lo)
    )
  }
  dimnames(out) <- list(rownames(G), names(tables))
  out
}

.new_likelihood_cache <- function(enabled = TRUE, max_entries = 4096L,
                                  max_bytes = 4 * 1024^2) {
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled))
    stop("cache_enabled must be one non-missing logical value", call. = FALSE)
  if (!is.numeric(max_entries) || length(max_entries) != 1L || is.na(max_entries) ||
      !is.finite(max_entries) || max_entries < 0 ||
      max_entries != floor(max_entries) || max_entries > .Machine$integer.max)
    stop("cache_entries must be one nonnegative whole number", call. = FALSE)
  if (!is.numeric(max_bytes) || length(max_bytes) != 1L || is.na(max_bytes) ||
      !is.finite(max_bytes) || max_bytes < 0)
    stop("cache_bytes must be one nonnegative finite number", call. = FALSE)
  cache <- new.env(parent = emptyenv())
  cache$enabled <- enabled
  cache$max_entries <- as.integer(max_entries)
  cache$max_bytes <- as.double(max_bytes)
  cache$descriptor <- new.env(parent = emptyenv())
  cache$pair <- new.env(parent = emptyenv())
  cache$order <- character()
  cache$bytes <- 0
  cache$stats <- list(hits = 0L, misses = 0L, insertions = 0L, evictions = 0L,
                      bypasses = 0L, current_size = 0L, current_bytes = 0)
  cache
}

.new_inference_context <- function(GenoM, tables,
                                   Complex = c("simp", "full"),
                                   cache_enabled = TRUE,
                                   cache_entries = 4096L,
                                   cache_bytes = 4 * 1024^2,
                                   include_duplicate = FALSE) {
  Complex <- match.arg(Complex)
  .validate_compact_locus_tables(tables)
  G <- .prepare_for_tables(GenoM, tables)
  registry <- .catalog_relationship_registry(Complex, include_duplicate = TRUE)
  static_registry <- registry[registry[["ordinary_fallback"]] &
    (isTRUE(include_duplicate) | registry$config_id != "MZ"), , drop = FALSE]
  delta <- do.call(rbind, unclass(static_registry$Delta))
  rownames(delta) <- static_registry$config_id
  colnames(delta) <- paste0("D", 1:9)

  ctx <- new.env(parent = emptyenv())
  ctx$schema_version <- "inference-context-v1"
  ctx$Complex <- Complex
  ctx$G <- G
  ctx$obsidx <- .prepare_observed_indices(G, tables)
  ctx$markers <- tables
  ctx$registry <- registry
  ctx$static_registry <- static_registry
  ctx$Delta <- delta
  ctx$ids <- rownames(G)
  ctx$unique_observed_states <- lapply(seq_along(tables), function(l)
    sort(unique(ctx$obsidx[, l])))
  ctx$likelihood_cache <- .new_likelihood_cache(
    cache_enabled, cache_entries, cache_bytes)
  class(ctx) <- "PaPeR_inference_context"
  for (name in c("schema_version", "Complex", "G", "obsidx", "markers", "registry",
                 "static_registry", "Delta", "ids", "unique_observed_states",
                 "likelihood_cache"))
    lockBinding(name, ctx)
  ctx
}

.context_rows <- function(ctx, x) {
  if (!inherits(ctx, "PaPeR_inference_context"))
    stop("ctx must be a PaPeR inference context", call. = FALSE)
  if (is.character(x)) {
    rows <- match(x, rownames(ctx$G))
    if (anyNA(rows))
      stop("unknown individual ID: ", paste(x[is.na(rows)], collapse = ", "),
           call. = FALSE)
    return(rows - 1L)
  }
  rows <- as.integer(x)
  if (anyNA(rows) || any(rows < 1L) || any(rows > nrow(ctx$G)))
    stop("individual row is out of range", call. = FALSE)
  rows - 1L
}

.score_unconditional_pairs <- function(ctx, i, j) {
  ii <- .context_rows(ctx, i)
  jj <- .context_rows(ctx, j)
  if (length(ii) != length(jj))
    stop("i and j must have equal length", call. = FALSE)
  # The ordinary production route intentionally bypasses the R-level cache:
  # one compiled call owns a whole pair batch and deduplicates locally.
  out <- score_basis_pairs_cpp(ctx$markers, ctx$obsidx, ii, jj, ctx$Delta)
  # The compiled entry point retains its test-only counters.  The higher-level
  # scorer remains a plain likelihood matrix for compatibility with callers.
  attr(out, "compact_batch_counters") <- NULL
  colnames(out) <- ctx$static_registry$config_id
  rownames(out) <- paste0(rownames(ctx$G)[ii + 1L], "|",
                          rownames(ctx$G)[jj + 1L])
  out
}

