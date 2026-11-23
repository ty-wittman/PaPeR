
test_that("simp omits the duplicate guard while retaining typed classes", {
  simp <- PaPeR:::.catalog_relationship_registry("simp")
  expect_setequal(simp$config_id, c("U", "PO", "FS", "HS", "GP", "FA", "HA"))
  expect_false("MZ" %in% simp$config_id)
  expect_identical(PaPeR:::.relationship_aliases("simp")$DEG2,
                   c("HS", "GP", "FA"))
})

.exact_pair_from_pedigree_complex <- function(topology, p) {
  checked <- PaPeR:::.validate_complex_topology(topology)
  k <- length(p)
  states <- list(list(genotypes = list(), probability = 1))
  ordered <- as.matrix(expand.grid(seq_len(k), seq_len(k)))
  for (founder in topology$founders) {
    out <- vector("list", length(states) * nrow(ordered))
    z <- 0L
    for (state in states) for (i in seq_len(nrow(ordered))) {
      z <- z + 1L
      next_state <- state
      next_state$genotypes[[founder]] <- ordered[i, ]
      next_state$probability <- state$probability *
        p[ordered[i, 1L]] * p[ordered[i, 2L]]
      out[[z]] <- next_state
    }
    states <- out
  }
  inheritance <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
  for (child in checked$order) {
    par <- topology$parents[[child]]
    out <- vector("list", length(states) * 4L)
    z <- 0L
    for (state in states) for (i in seq_len(4L)) {
      z <- z + 1L
      next_state <- state
      next_state$genotypes[[child]] <- c(
        state$genotypes[[par[1L]]][inheritance[i, 1L]],
        state$genotypes[[par[2L]]][inheritance[i, 2L]])
      next_state$probability <- state$probability / 4
      out[[z]] <- next_state
    }
    states <- out
  }
  G <- k * (k + 1L) / 2L
  index <- function(a, b) {
    if (a > b) { z <- a; a <- b; b <- z }
    (a - 1L) * k - (a - 1L) * (a - 2L) / 2L + (b - a) + 1L
  }
  ans <- matrix(0, G, G)
  for (state in states) {
    gx <- state$genotypes[[topology$focal[1L]]]
    gy <- state$genotypes[[topology$focal[2L]]]
    ans[index(gx[1L], gx[2L]), index(gy[1L], gy[2L])] <-
      ans[index(gx[1L], gx[2L]), index(gy[1L], gy[2L])] + state$probability
  }
  ans
}

test_that("representative full topologies equal independent genotype enumeration", {
  catalog <- PaPeR:::.catalog_topology_catalog()
  reg <- PaPeR:::.catalog_relationship_registry("full", include_duplicate = TRUE)
  cases <- list(
    HS_PO = c(.55, .30, .15),
    HS_HA = c(.65, .35),
    FS_PARENTS_FS = c(.55, .30, .15))
  for (nm in names(cases)) {
    p <- cases[[nm]]
    exact <- .exact_pair_from_pedigree_complex(catalog[[nm]], p)
    delta <- reg$Delta[[match(nm, reg$config_id)]]
    tab <- PaPeR:::build_locus_table(p, Err = 1e-12)
    expect_equal(dense_relationship_oracle(tab, delta), exact,
                 tolerance = 8e-14, info = nm)
  }
})

.emission_ref_complex <- function(tab, true, observed) {
  true <- sort(as.integer(true)); observed <- sort(as.integer(observed))
  if (identical(true, observed)) return(1 - tab$Err)
  if (tab$Err == 0) return(0)
  overlap <- sum(pmin(tabulate(true, nbins = tab$k),
                      tabulate(observed, nbins = tab$k)))
  cls <- if (overlap == 0L) "w_Far" else if (observed[1L] == observed[2L])
    "w_Hom" else "w_Het"
  z <- if (true[1L] == true[2L]) tab$Z_hom else tab$Z_het
  tab$Err * tab$ErrProfile[[cls]] / z
}

.exact_observation_probability_complex <- function(topology, p, observed, tab) {
  checked <- PaPeR:::.validate_complex_topology(topology)
  k <- length(p)
  states <- list(list(genotypes = list(), probability = 1))
  ordered <- as.matrix(expand.grid(seq_len(k), seq_len(k)))
  for (founder in topology$founders) {
    out <- vector("list", length(states) * nrow(ordered)); z <- 0L
    for (state in states) for (i in seq_len(nrow(ordered))) {
      z <- z + 1L; next_state <- state
      next_state$genotypes[[founder]] <- ordered[i, ]
      next_state$probability <- state$probability *
        p[ordered[i, 1L]] * p[ordered[i, 2L]]
      out[[z]] <- next_state
    }
    states <- out
  }
  inheritance <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
  for (child in checked$order) {
    par <- topology$parents[[child]]
    out <- vector("list", length(states) * 4L); z <- 0L
    for (state in states) for (i in seq_len(4L)) {
      z <- z + 1L; next_state <- state
      next_state$genotypes[[child]] <- c(
        state$genotypes[[par[1L]]][inheritance[i, 1L]],
        state$genotypes[[par[2L]]][inheritance[i, 2L]])
      next_state$probability <- state$probability / 4
      out[[z]] <- next_state
    }
    states <- out
  }
  sum(vapply(states, function(state) {
    value <- state$probability
    for (node in names(observed))
      value <- value * .emission_ref_complex(
        tab, state$genotypes[[node]], observed[[node]])
    value
  }, numeric(1)))
}

.complex_validation_fixture <- function(topology, p, n_loci = 1L) {
  checked <- PaPeR:::.validate_complex_topology(topology)
  ids <- checked$nodes
  calls <- vapply(seq_along(ids), function(i) {
    a <- (i - 1L) %% length(p) + 1L
    b <- i %% length(p) + 1L
    paste(sort(c(a, b)), collapse = "/")
  }, character(1))
  G <- matrix(rep(calls, n_loci), nrow = length(ids), ncol = n_loci,
              dimnames = list(ids, paste0("L", seq_len(n_loci))))
  tables <- PaPeR:::build_tables(
    setNames(rep(list(p), n_loci), colnames(G)), Err = 0.01)
  state <- data.frame(id = ids, dam = NA_character_, sire = NA_character_,
                      stringsAsFactors = FALSE)
  for (child in names(topology$parents)) {
    state$dam[match(child, ids)] <- topology$parents[[child]][1L]
    state$sire[match(child, ids)] <- topology$parents[[child]][2L]
  }
  list(ids = ids, G = G, tables = tables, state = state,
       node_map = setNames(ids, checked$nodes), checked = checked)
}

.complex_validation_observations <- function(G, locus = 1L) {
  setNames(lapply(seq_len(nrow(G)), function(i) {
    as.integer(strsplit(G[i, locus], "/", fixed = TRUE)[[1L]])
  }), rownames(G))
}

# Independent vectorized oracle used only for the k=3 high-risk cases.  It
# enumerates founder allele copies and Mendelian copy choices directly, rather
# than calling the production inheritance-map builder.
.emission_probability_complex_vectorized <- function(tab, true, observed) {
  true <- as.matrix(true)
  observed <- sort(as.integer(observed))
  true_a <- pmin(true[, 1L], true[, 2L])
  true_b <- pmax(true[, 1L], true[, 2L])
  same <- true_a == observed[1L] & true_b == observed[2L]
  if (tab$Err == 0) return(as.numeric(same))
  overlap <- ifelse(true_a == true_b,
                    as.integer(observed[1L] == true_a) +
                      as.integer(observed[2L] == true_a),
                    as.integer(observed[1L] == true_a |
                               observed[2L] == true_a) +
                      as.integer(observed[1L] == true_b |
                                 observed[2L] == true_b))
  weight <- ifelse(overlap == 0L, tab$ErrProfile[["w_Far"]],
                   if (observed[1L] == observed[2L])
                     tab$ErrProfile[["w_Hom"]] else tab$ErrProfile[["w_Het"]])
  out <- tab$Err * weight / ifelse(true_a == true_b, tab$Z_hom, tab$Z_het)
  out[same] <- 1 - tab$Err
  out
}

.exact_observation_probability_complex_vectorized <- function(topology, p,
                                                              observed, tab) {
  founders <- topology$founders
  pending <- names(topology$parents)
  done <- founders
  order <- character()
  while (length(pending)) {
    ready <- pending[vapply(topology$parents[pending], function(x)
      all(x %in% done), logical(1))]
    if (!length(ready)) stop("oracle topology cycle")
    order <- c(order, ready)
    done <- c(done, ready)
    pending <- setdiff(pending, ready)
  }
  source_maps <- list(setNames(lapply(seq_along(founders), function(i)
    c(2L * i - 1L, 2L * i)), founders))
  choices <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
  for (child in order) {
    parents <- topology$parents[[child]]
    next_maps <- vector("list", length(source_maps) * 4L)
    z <- 0L
    for (source_map in source_maps) for (choice in seq_len(nrow(choices))) {
      z <- z + 1L
      next_map <- source_map
      next_map[[child]] <- c(source_map[[parents[1L]]][choices[choice, 1L]],
                             source_map[[parents[2L]]][choices[choice, 2L]])
      next_maps[[z]] <- next_map
    }
    source_maps <- next_maps
  }
  founder_assignments <- as.matrix(expand.grid(rep(list(seq_len(length(p))),
                                                   2L * length(founders))))
  founder_probability <- apply(founder_assignments, 1L, function(x)
    prod(p[x]))
  log_weight <- -length(order) * log(4)
  sum(vapply(source_maps, function(source_map) {
    probability <- founder_probability * exp(log_weight)
    for (node in names(observed)) {
      true <- founder_assignments[, source_map[[node]], drop = FALSE]
      probability <- probability * .emission_probability_complex_vectorized(
        tab, true, observed[[node]])
    }
    sum(probability)
  }, numeric(1)))
}





