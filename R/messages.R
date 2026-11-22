# Reusable production messages for typed-parent Stage-1 search.
#
# Current Sequoia publishes individual genotype messages (LindX) in
# UpdateAllProbs(), then uses bounded CalcLind/ParProb/setParTmp calculations
# while SelectParent() temporarily changes local parent slots.  The production
# route below follows that boundary: one immutable ancestry-message snapshot is
# published per message_revision, while focal and nuclear-family contractions
# are derived lazily from those frozen boundary messages.

.message_schema <- "simp-message-snapshot-v2"

.message_log_normalize <- function(x) {
  z <- .logsumexp(x)
  if (identical(z, -Inf)) return(x)
  x - z
}

.message_log_transmitted <- function(message_ln, transmission_log_weight) {
  k <- ncol(transmission_log_weight)
  vapply(seq_len(k), function(a)
    .logsumexp(message_ln + transmission_log_weight[, a]), numeric(1))
}

.message_child_from_alleles <- function(dam_allele_ln, sire_allele_ln,
                                     child_index) {
  k <- nrow(child_index)
  out <- rep(-Inf, k * (k + 1L) / 2L)
  for (a in seq_len(k)) for (b in a:k) {
    g <- child_index[a, b]
    if (a == b) out[g] <- dam_allele_ln[a] + sire_allele_ln[b] else
      out[g] <- .logsumexp(c(dam_allele_ln[a] + sire_allele_ln[b],
                             dam_allele_ln[b] + sire_allele_ln[a]))
  }
  out
}

.message_child_from_parent_messages <- function(locus, dam_ln = NULL,
                                             sire_ln = NULL) {
  dam <- if (is.null(dam_ln)) locus$log_p else
    .message_log_transmitted(dam_ln, locus$transmission_log_weight)
  sire <- if (is.null(sire_ln)) locus$log_p else
    .message_log_transmitted(sire_ln, locus$transmission_log_weight)
  .message_child_from_alleles(dam, sire, locus$child_index)
}

.message_locus_data <- function(ctx, locus, ids,
                             child_observation_cache_enabled = TRUE) {
  tab <- ctx$markers[[locus]]
  k <- as.integer(tab$k); G <- as.integer(tab$G)
  geno <- tab$geno
  child_index <- matrix(0L, k, k)
  for (a in seq_len(k)) for (b in a:k)
    child_index[a, b] <- (a - 1L) * k - (a - 1L) * (a - 2L) / 2L +
      (b - a) + 1L
  log_p <- log(as.numeric(tab$p))
  prior <- vapply(seq_len(G), function(g) {
    a <- geno[g, 1L]; b <- geno[g, 2L]
    if (a == b) 2 * log_p[a] else log(2) + log_p[a] + log_p[b]
  }, numeric(1))
  transmission_log_weight <- matrix(-Inf, G, k)
  for (g in seq_len(G)) {
    a <- geno[g, 1L]; b <- geno[g, 2L]
    transmission_log_weight[g, a] <- log(if (a == b) 1 else .5)
    transmission_log_weight[g, b] <- log(if (a == b) 1 else .5)
  }
  observation <- lapply(seq_len(nrow(ctx$G)), function(row) {
    observed <- ctx$obsidx[row, locus]
    if (observed < 0L) return(rep(0, G))
    emission <- emission_probabilities_cpp(
      tab, as.integer(seq_len(G) - 1L), rep.int(as.integer(observed), G))
    .log_probability_ln(emission, "simp observed-genotype emission")
  })
  names(observation) <- ids
  list(tab = tab, k = k, G = G, geno = geno, child_index = child_index,
       log_p = log_p, prior_ln = prior,
       transmission_log_weight = transmission_log_weight,
       observation_ln = observation,
       parent_joint_cache = new.env(parent = emptyenv()),
       child_observation_cache = new.env(parent = emptyenv()),
       child_observation_cache_enabled = isTRUE(
         child_observation_cache_enabled))
}

.message_topological_order <- function(parents) {
  n <- nrow(parents)
  indegree <- rowSums(parents > 0L)
  reverse <- vector("list", n)
  for (child in seq_len(n)) for (parent in parents[child, ])
    if (parent > 0L) reverse[[parent]] <- c(reverse[[parent]], child)
  queue <- which(indegree == 0L); out <- integer()
  while (length(queue)) {
    node <- queue[1L]; queue <- queue[-1L]
    out <- c(out, node)
    for (child in reverse[[node]]) {
      indegree[child] <- indegree[child] - 1L
      if (indegree[child] == 0L) queue <- c(queue, child)
    }
  }
  if (length(out) != n) stop("simp snapshot requires an acyclic parent graph",
                             call. = FALSE)
  out
}

.message_ancestor_distances <- function(parents, node) {
  distance <- setNames(0L, as.character(node))
  queue <- as.integer(node)
  while (length(queue)) {
    current <- queue[1L]; queue <- queue[-1L]
    d <- distance[[as.character(current)]]
    for (parent in parents[current, ]) if (parent > 0L) {
      key <- as.character(parent)
      old <- unname(distance[match(key, names(distance))])
      if (!length(old) || is.na(old) || old > d + 1L) {
        distance[[key]] <- d + 1L
        queue <- c(queue, parent)
      }
    }
  }
  distance
}

# Source-defined ordinary relationship correction for connected typed parents.
# Direct PO, FS and HS paths are exact IBD classes.  Repeated deeper ancestral
# paths use the same bounded kinship approximation represented by simp ordinary
# alternatives; it prevents an invalid product of connected marginals without
# introducing a general pedigree-inference engine.
.message_parent_ibd <- function(parents, first, second) {
  if (first == second) return(list(k0 = 0, k1 = 0, k2 = 1,
                                  relation = "duplicate"))
  p1 <- parents[first, parents[first, ] > 0L]
  p2 <- parents[second, parents[second, ] > 0L]
  shared <- intersect(p1, p2)
  if (length(p1) == 2L && length(p2) == 2L && length(shared) == 2L)
    return(list(k0 = .25, k1 = .5, k2 = .25, relation = "full_siblings"))
  d1 <- .message_ancestor_distances(parents, first)
  d2 <- .message_ancestor_distances(parents, second)
  direct12 <- unname(d1[match(as.character(second), names(d1))])
  direct21 <- unname(d2[match(as.character(first), names(d2))])
  direct <- c(if (!length(direct12) || is.na(direct12)) Inf else direct12,
              if (!length(direct21) || is.na(direct21)) Inf else direct21)
  if (any(is.finite(direct))) {
    d <- min(direct)
    k1 <- min(1, 2^(1 - d))
    return(list(k0 = 1 - k1, k1 = k1, k2 = 0,
                relation = if (d == 1L) "parent_offspring" else
                  paste0("direct_ancestor_degree_", d)))
  }
  common <- intersect(names(d1), names(d2))
  common <- setdiff(common, c(as.character(first), as.character(second)))
  if (!length(common))
    return(list(k0 = 1, k1 = 0, k2 = 0, relation = "unconnected"))
  r <- sum(vapply(common, function(a)
    2^(-(d1[[a]] + d2[[a]])), numeric(1)))
  k1 <- min(1, 2 * r)
  list(k0 = 1 - k1, k1 = k1, k2 = 0,
       relation = if (length(shared) == 1L) "half_siblings" else
         "shared_ancestry_search_approximation")
}

.message_parent_joint_ln <- function(locus, first_ln, second_ln, ibd) {
  # The compact C++ primitive performs the same independent or IBD-weighted
  # natural-log construction and normalization without repeated R matrix
  # allocation inside final family-boundary recursion.
  parent_joint_log_cpp(
    locus$prior_ln, locus$log_p, locus$child_index,
    as.numeric(first_ln), as.numeric(second_ln),
    ibd$k0, ibd$k1, ibd$k2)
}

.message_child_from_parent_joint <- function(locus, joint_ln) {
  dam_allele <- locus$transmission_log_weight
  sire_allele <- locus$transmission_log_weight
  allele_joint <- matrix(-Inf, locus$k, locus$k)
  for (a in seq_len(locus$k)) for (b in seq_len(locus$k)) {
    terms <- numeric()
    for (gd in seq_len(locus$G)) for (gs in seq_len(locus$G))
      terms <- c(terms, joint_ln[gd, gs] + dam_allele[gd, a] +
                   sire_allele[gs, b])
    allele_joint[a, b] <- .logsumexp(terms)
  }
  out <- rep(-Inf, locus$G)
  for (a in seq_len(locus$k)) for (b in a:locus$k) {
    g <- locus$child_index[a, b]
    out[g] <- if (a == b) allele_joint[a, b] else
      .logsumexp(c(allele_joint[a, b], allele_joint[b, a]))
  }
  out
}

.message_build_ancestry_messages <- function(locus, parents, order, parent_ibd) {
  n <- nrow(parents); out <- vector("list", n); likelihood <- numeric(n)
  founder <- NULL  # the same inherited message for every parentless individual
  for (child in order) {
    dam <- parents[child, 1L]; sire <- parents[child, 2L]
    if (dam == 0L && sire == 0L) {
      if (is.null(founder)) founder <- .message_child_from_parent_messages(locus)
      inherited <- founder
    } else if (dam > 0L && sire > 0L) {
      ibd <- parent_ibd[[child]]
      if (ibd$k1 > 0 || ibd$k2 > 0) {
        joint <- .message_parent_joint_ln(locus, out[[dam]], out[[sire]], ibd)
        inherited <- .message_child_from_parent_joint(locus, joint)
      } else inherited <- .message_child_from_parent_messages(
        locus, out[[dam]], out[[sire]])
    } else inherited <- .message_child_from_parent_messages(
      locus, if (dam > 0L) out[[dam]] else NULL,
      if (sire > 0L) out[[sire]] else NULL)
    unscaled <- inherited + locus$observation_ln[[child]]
    likelihood[child] <- .logsumexp(unscaled)
    out[[child]] <- .message_log_normalize(unscaled)
  }
  list(genotype_ln = out, individual_likelihood_ln = likelihood)
}

.message_snapshot_families <- function(ids, parents) {
  assigned <- which(rowSums(parents > 0L) > 0L)
  if (!length(assigned)) return(list(families = list(), by_parent =
    setNames(vector("list", length(ids)), ids)))
  key <- paste(parents[assigned, 1L], parents[assigned, 2L], sep = ":")
  groups <- split(assigned, key)
  families <- lapply(names(groups), function(k) {
    child <- groups[[k]]
    list(key = k, dam = parents[child[1L], 1L],
         sire = parents[child[1L], 2L], children = child)
  })
  names(families) <- names(groups)
  by_parent <- vector("list", length(ids)); names(by_parent) <- ids
  for (i in seq_along(families)) for (parent in
    unique(c(families[[i]]$dam, families[[i]]$sire)))
    if (parent > 0L) by_parent[[parent]] <- c(by_parent[[parent]], i)
  list(families = families, by_parent = by_parent)
}

.message_snapshot_lock <- function(snapshot) {
  for (nm in ls(snapshot, all.names = TRUE)) lockBinding(nm, snapshot)
  snapshot
}

.message_compact_snapshot_inputs <- function(snapshot) {
  ids <- snapshot$genotype_state_metadata$ids
  list(
    n_individuals = as.integer(length(ids)),
    message_revision = as.integer(snapshot$message_revision),
    prior_ln = lapply(snapshot$loci, `[[`, "prior_ln"),
    log_p = lapply(snapshot$loci, `[[`, "log_p"),
    transmission_log_weight = lapply(
      snapshot$loci, `[[`, "transmission_log_weight"),
    child_index = lapply(snapshot$loci, `[[`, "child_index"),
    observations_ln = lapply(snapshot$loci, function(locus)
      do.call(cbind, locus$observation_ln)),
    messages_ln = lapply(snapshot$messages, function(messages)
      do.call(cbind, messages)),
    genotype_alleles = snapshot$context_identity$G,
    observed_genotypes = snapshot$context_identity$obsidx,
    id_index = setNames(as.integer(seq_along(ids)), ids)
  )
}


#' Build the initial message snapshot.
#' @noRd
.build_message_snapshot <- function(ctx, state,
                                         child_observation_cache_enabled = TRUE,
                                         purpose = c("search", "final-output-scoring")) {
  purpose <- match.arg(purpose)
  if (!inherits(ctx, "PaPeR_inference_context") ||
      !ctx$Complex %in% c("simp", "full"))
    stop("the Stage-1 message snapshot requires a Complex='simp' or 'full' context",
         call. = FALSE)
  if (!inherits(state, "PaPeR_PedigreeState"))
    stop("the simp message snapshot requires a ParentageState", call. = FALSE)
  ids <- .parentage_state_ids(state)
  if (!identical(ids, rownames(ctx$G)))
    stop("state IDs must match the inference context", call. = FALSE)
  snapshot <- new.env(parent = emptyenv())
  snapshot$schema_version <- .message_schema
  snapshot$context_identity <- ctx
  snapshot$purpose <- purpose
  parents <- state$parents
  order <- .message_topological_order(parents)
  parent_ibd <- lapply(seq_len(nrow(parents)), function(child) {
    if (all(parents[child, ] > 0L))
      .message_parent_ibd(parents, parents[child, 1L], parents[child, 2L]) else
        list(k0 = 1, k1 = 0, k2 = 0, relation = "unconnected")
  })
  loci <- lapply(seq_along(ctx$markers), function(l)
    .message_locus_data(
      ctx, l, ids,
      child_observation_cache_enabled = child_observation_cache_enabled))
  ancestry <- lapply(seq_along(loci), function(l)
    .message_build_ancestry_messages(
      loci[[l]], parents = parents, order = order, parent_ibd = parent_ibd))
  names(loci) <- names(ctx$markers); names(ancestry) <- names(ctx$markers)
  routing <- .message_snapshot_families(ids, parents)
  snapshot$context_token <- paste0(ctx$schema_version, ":",
    length(ids), ":", length(loci), ":", state$message_revision)
  snapshot$purpose <- purpose
  snapshot$prepared_locus_metadata <- lapply(loci, function(x)
    list(k = x$k, G = x$G,
         maximum_temporary_joint_cells = x$G^2,
         maximum_temporary_joint_bytes = 8 * x$G^2))
  snapshot$genotype_state_metadata <- list(
    ids = ids, id_index = setNames(as.integer(seq_along(ids)), ids),
    G_dim = dim(ctx$G), obs_dim = dim(ctx$obsidx),
    locus_names = attr(ctx$G, "locus_names", exact = TRUE),
    n_alleles = attr(ctx$G, "n_alleles", exact = TRUE))
  snapshot$message_revision <- as.integer(state$message_revision)
  snapshot$graph_revision <- as.integer(state$graph_revision)
  snapshot$message_signature <- paste0("message-revision=", state$message_revision)
  snapshot$parents_at_construction <- parents
  snapshot$topological_order <- order
  snapshot$loci <- loci
  snapshot$messages <- lapply(ancestry, `[[`, "genotype_ln")
  snapshot$individual_likelihood_ln <- lapply(
    ancestry, `[[`, "individual_likelihood_ln")
  snapshot$families <- routing$families
  snapshot$families_by_parent <- routing$by_parent
  snapshot$graph_provenance <- list(
    message_revision = as.integer(state$message_revision),
    graph_revision = as.integer(state$graph_revision),
    parents_at_construction = parents,
    purpose = purpose)
  snapshot$bounded_cache_metadata <- list(
    persistent_dense_parent_pair_tables = 0L,
    cached_loci = length(loci),
    cached_individual_locus_messages = length(ids) * length(loci),
    child_observation_cache_enabled = isTRUE(child_observation_cache_enabled),
    maximum_temporary_joint_cells = max(vapply(loci, function(x) x$G^2,
                                                numeric(1))),
    maximum_temporary_joint_bytes = max(vapply(loci, function(x) 8 * x$G^2,
                                                numeric(1))))
  snapshot$diagnostic_metadata <- list(
    topology_status = "relationship_source_mapped",
    supported_motifs = c("repeated_full_siblings", "maternal_half_siblings",
      "paternal_half_siblings", "repeated_breeders", "multigenerational_typed",
      "direct_ancestral_mating", "deeper_shared_ancestry_search_approximation"),
    exact_equality_boundary = c("relationship_acyclic_parentage",
      "repeated_full_sibling_families", "half_sibling_families",
      "direct_parent_offspring_co_parents"),
    approximation_boundary =
      "multiply connected and deeper shared-ancestry simp routing",
    publication = if (identical(purpose, "final-output-scoring"))
      "final-output-only" else "frozen within one parentage pass",
    snapshot_purpose = purpose,
    local_scratch = "bounded current-family calculations",
    log_base = "e")
  if (purpose %in% c("search", "final-output-scoring")) {
    compact_inputs <- .message_compact_snapshot_inputs(snapshot)
    snapshot$compiled <- build_snapshot_cpp(compact_inputs)
  } else snapshot$compiled <- NULL
  class(snapshot) <- c("PaPeR_message_snapshot",
                       "PaPeR_message_provider")
  .message_snapshot_lock(snapshot)
  snapshot
}
