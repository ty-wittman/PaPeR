# Independent dense test-only oracles.  Production constructors must never
# return these matrices.

compact_legacy_delta <- function() {
  out <- rbind(
    U = c(0, 0, 0, 0, 0, 0, 0, 0, 1),
    PO = c(0, 0, 0, 0, 0, 0, 0, 1, 0),
    FS = c(0, 0, 0, 0, 0, 0, .25, .5, .25),
    DEG2 = c(0, 0, 0, 0, 0, 0, 0, .5, .5),
    HA = c(0, 0, 0, 0, 0, 0, 0, .25, .75),
    MZ = c(0, 0, 0, 0, 0, 0, 1, 0, 0))
  colnames(out) <- paste0("D", 1:9)
  out
}

dense_emission_oracle_compact <- function(tab) {
  out <- matrix(0, tab$G, tab$G)
  for (truth in seq_len(tab$G)) for (observed in seq_len(tab$G)) {
    gt <- tab$geno[truth, ]; go <- tab$geno[observed, ]
    if (truth == observed) {
      out[truth, observed] <- 1 - tab$Err
    } else if (tab$Err > 0) {
      overlap <- sum(pmin(tabulate(gt, nbins = tab$k),
                          tabulate(go, nbins = tab$k)))
      cls <- if (overlap == 0L) "w_Far" else if (go[1L] == go[2L]) "w_Hom" else "w_Het"
      z <- if (gt[1L] == gt[2L]) tab$Z_hom else tab$Z_het
      out[truth, observed] <- tab$Err * tab$ErrProfile[[cls]] / z
    }
  }
  out
}

dense_jacquard_basis_oracle <- function(tab) {
  groups <- list(
    c(1L, 1L, 1L, 1L), c(1L, 1L, 2L, 2L), c(1L, 1L, 1L, 2L),
    c(1L, 1L, 2L, 3L), c(1L, 2L, 1L, 1L), c(1L, 2L, 3L, 3L),
    c(1L, 2L, 1L, 2L), c(1L, 2L, 1L, 3L), c(1L, 2L, 3L, 4L))
  out <- lapply(groups, function(group) {
    ans <- matrix(0, tab$G, tab$G)
    for (gx in seq_len(tab$G)) for (gy in seq_len(tab$G)) {
      x <- tab$geno[gx, ]; y <- tab$geno[gy, ]
      x_orders <- if (x[1L] == x[2L]) list(x) else list(x, rev(x))
      y_orders <- if (y[1L] == y[2L]) list(y) else list(y, rev(y))
      for (xo in x_orders) for (yo in y_orders) {
        allele <- c(xo, yo)
        first <- match(group, group)
        valid <- all(allele == allele[first])
        if (valid) ans[gx, gy] <- ans[gx, gy] + prod(tab$p[allele[which(first == seq_along(group))]])
      }
    }
    ans
  })
  names(out) <- paste0("D", 1:9)
  out
}

dense_relationship_oracle <- function(tab, delta) {
  basis <- dense_jacquard_basis_oracle(tab)
  Reduce(`+`, Map(`*`, as.numeric(delta), basis))
}

dense_observed_relationship_oracle <- function(tab, delta) {
  e <- dense_emission_oracle_compact(tab)
  crossprod(e, dense_relationship_oracle(tab, delta) %*% e)
}

compact_pair_state_founder_oracle <- function(tab, left, right) {
  emission <- function(observed, a, b) {
    if (is.na(observed) || observed < 0L) return(1)
    go <- tab$geno[observed + 1L, ]
    gt <- sort(c(a, b))
    if (identical(gt, go)) return(1 - tab$Err)
    if (tab$Err == 0) return(0)
    overlap <- sum(pmin(tabulate(gt, nbins = tab$k),
                        tabulate(go, nbins = tab$k)))
    cls <- if (overlap == 0L) "w_Far" else if (go[1L] == go[2L]) "w_Hom" else "w_Het"
    tab$Err * tab$ErrProfile[[cls]] / if (gt[1L] == gt[2L]) tab$Z_hom else tab$Z_het
  }
  d <- function(observed) vapply(seq_len(tab$k), function(a) emission(observed, a, a), numeric(1))
  r <- function(observed) vapply(seq_len(tab$k), function(a)
    sum(tab$p * vapply(seq_len(tab$k), function(b) emission(observed, a, b), numeric(1))), numeric(1))
  di <- d(left); dj <- d(right); ri <- r(left); rj <- r(right)
  hi <- sum(tab$p * di); hj <- sum(tab$p * dj)
  mi <- sum(tab$p * ri); mj <- sum(tab$p * rj)
  a_left <- outer(seq_len(tab$k), seq_len(tab$k),
                  Vectorize(function(a, b) emission(left, a, b)))
  a_right <- outer(seq_len(tab$k), seq_len(tab$k),
                   Vectorize(function(a, b) emission(right, a, b)))
  d7 <- sum(outer(tab$p, tab$p) * a_left * a_right)
  c(D1 = sum(tab$p * di * dj), D2 = hi * hj,
    D3 = sum(tab$p * di * rj), D4 = hi * mj,
    D5 = sum(tab$p * ri * dj), D6 = mi * hj,
    D7 = d7, D8 = sum(tab$p * ri * rj), D9 = mi * mj)
}

compact_recursive_dense_matrix <- function(x, G) {
  if (is.matrix(x) && is.numeric(x) && identical(dim(x), rep(as.integer(G), 2L)))
    return(TRUE)
  if (!is.list(x)) return(FALSE)
  any(vapply(x, compact_recursive_dense_matrix, logical(1), G = G))
}
