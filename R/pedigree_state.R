# Authoritative typed-only parentage state and transactional configuration patches.

.empty_edges <- function() {
  data.frame(child = character(), side = character(), parent = character(),
             stringsAsFactors = FALSE)
}

.parentage_state_ids <- function(state) state$nodes$id

.parentage_parent_ids <- function(state) {
  ids <- .parentage_state_ids(state)
  out <- matrix(NA_character_, nrow(state$parents), 2L,
                dimnames = list(ids, c("dam", "sire")))
  used <- state$parents > 0L
  out[used] <- ids[state$parents[used]]
  out
}

.current_parent_configuration <- function(state, child) {
  child <- match(as.character(child)[1L], .parentage_state_ids(state))
  if (is.na(child)) stop("parent configuration requested for an unknown ID",
                         call. = FALSE)
  .parentage_parent_ids(state)[child, ]
}

# Deterministic read-only view. The returned data frame is never stored as a
# second mutable pedigree representation.

.parentage_state_data_frame <- function(state) {
  p <- .parentage_parent_ids(state)
  data.frame(id = rownames(p), dam = unname(p[, "dam"]),
             sire = unname(p[, "sire"]), stringsAsFactors = FALSE)
}

.parent_indices <- function(ids, values, field) {
  values <- as.character(values)
  values[is.na(values) | !nzchar(values) | values == "0"] <- NA_character_
  bad <- unique(values[!is.na(values) & !values %in% ids])
  if (length(bad))
    stop(field, " names IDs absent from the typed inference context: ",
         paste(sort(bad), collapse = ", "), call. = FALSE)
  out <- match(values, ids)
  out[is.na(out)] <- 0L
  as.integer(out)
}


.derive_parentage_indexes <- function(ids, parents, sex) {
  reverse <- setNames(vector("list", length(ids)), ids)
  for (child in seq_len(nrow(parents))) for (side in seq_len(2L)) {
    parent <- parents[child, side]
    if (parent > 0L) reverse[[parent]] <- sort(unique(c(reverse[[parent]], child)))
  }
  role_counts <- matrix(0L, length(ids), 2L,
                        dimnames = list(ids, c("dam", "sire")))
  role <- setNames(rep.int(0L, length(ids)), ids)
  for (side in seq_len(2L)) {
    used <- parents[, side]; used <- used[used > 0L]
    if (length(used)) role_counts[, side] <- tabulate(used, nbins = length(ids))
    used <- unique(used)
    unknown <- used[is.na(sex[used])]
    if (length(unknown)) role[unknown] <- side
  }
  list(reverse_offspring = reverse, inferred_role = role,
       role_counts = role_counts)
}

.new_parentage_state_empty <- function(ids, genotype_rows, life_history,
                                       min_age_parent, max_age_parent,
                                       ctx = NULL) {
  state <- new.env(parent = emptyenv())
  state$nodes <- data.frame(
    id = ids, node_type = "individual", genotype_row = as.integer(genotype_rows),
    active = TRUE, stringsAsFactors = FALSE)
  state$id_to_index <- setNames(seq_along(ids), ids)
  state$parents <- matrix(0L, length(ids), 2L,
                          dimnames = list(ids, c("dam", "sire")))
  state$sex <- setNames(as.integer(life_history$Sex), ids)
  state$life_history <- life_history
  state$birth_year_masses <- .life_history_masses(life_history)
  state$birth_year_support <- setNames(lapply(ids, function(id)
    .birth_year_support(life_history, id)), ids)
  expected_year <- vapply(state$birth_year_masses, function(x)
    if (is.null(x)) Inf else sum(as.numeric(names(x)) * x), numeric(1))
  state$focal_order <- ids[order(expected_year, ids)]
  state$year_last <- setNames(as.integer(life_history$Year.last), ids)
  state$min_age_parent <- setNames(as.integer(min_age_parent), c("M", "P"))
  state$max_age_parent <- setNames(as.integer(max_age_parent), c("M", "P"))
  state$graph_revision <- 0L
  state$node_revision <- setNames(integer(length(ids)), ids)
  state$component_revision <- setNames(integer(length(ids)), ids)
  state$message_revision <- 0L
  state$current_messages <- NULL
  state$message_snapshot <- NULL
  state$node_likelihood_ln <- setNames(rep(NA_real_, length(ids)), ids)
  state$total_likelihood_ln <- NA_real_
  state$candidate_pool <- NULL
  state$active_shortlist_diagnostics <- list()
  state$unresolved_diagnostics <- list()
  state$accepted_move_history <- list()
  state$round_history <- list()
  state$state_overlay_cache <- new.env(parent = emptyenv())
  state$configuration_score_cache <- new.env(parent = emptyenv())
  state$message_score_cache <- new.env(parent = emptyenv())
  state$ctx <- ctx
  indexes <- .derive_parentage_indexes(ids, state$parents, state$sex)
  state$reverse_offspring <- indexes$reverse_offspring
  state$inferred_role <- indexes$inferred_role
  state$role_counts <- indexes$role_counts
  class(state) <- c("PaPeR_ParentageState", "PaPeR_PedigreeState")
  state
}

.normalize_supplied_pedigree <- function(pedigree, ids) {
  out <- matrix(0L, length(ids), 2L, dimnames = list(ids, c("dam", "sire")))
  if (is.null(pedigree)) return(out)
  if (!is.data.frame(pedigree) || !all(c("id", "dam", "sire") %in% names(pedigree)))
    stop("supplied pedigree must be a data.frame with id, dam, and sire",
         call. = FALSE)
  if (!nrow(pedigree)) return(out)
  child <- as.character(pedigree$id)
  if (anyNA(child) || any(!nzchar(child)) || anyDuplicated(child))
    stop("supplied pedigree IDs must be unique nonempty strings", call. = FALSE)
  if (any(!child %in% ids))
    stop("supplied pedigree contains a child absent from the typed inference context",
         call. = FALSE)
  rows <- match(child, ids)
  out[rows, "dam"] <- .parent_indices(ids, pedigree$dam, "supplied dam")
  out[rows, "sire"] <- .parent_indices(ids, pedigree$sire, "supplied sire")
  out
}

.parentage_cycle <- function(parents) {
  colour <- integer(nrow(parents))
  visit <- function(node) {
    if (colour[node] == 1L) return(TRUE)
    if (colour[node] == 2L) return(FALSE)
    colour[node] <<- 1L
    for (parent in parents[node, ]) if (parent > 0L && visit(parent)) return(TRUE)
    colour[node] <<- 2L
    FALSE
  }
  any(vapply(seq_len(nrow(parents)), visit, logical(1)))
}

.consistency_result <- function(valid, code = "ok", reason = "",
                                affected_ids = character()) {
  list(valid = isTRUE(valid), code = as.character(code),
       reason = as.character(reason),
       affected_ids = sort(unique(as.character(affected_ids))))
}

.check_parent_matrix <- function(state, parents) {
  invalid <- function(code, reason, rows = integer())
    .consistency_result(FALSE, code, reason, .parentage_state_ids(state)[rows])
  ids <- .parentage_state_ids(state)
  if (!is.matrix(parents) || !identical(dim(parents), c(length(ids), 2L)) ||
      anyNA(parents) || any(parents < 0L) || any(parents > length(ids)))
    return(invalid("invalid_parent_matrix", "authoritative parent matrix is invalid"))
  child <- seq_len(nrow(parents))
  if (any(parents[, "dam"] == child & parents[, "dam"] > 0L) ||
      any(parents[, "sire"] == child & parents[, "sire"] > 0L))
    return(invalid("self_parentage", "self-parentage is not permitted",
                   which(parents[, "dam"] == child | parents[, "sire"] == child)))
  same <- which(parents[, "dam"] > 0L & parents[, "dam"] == parents[, "sire"])
  if (length(same)) return(invalid(
    "same_parent_both_sides", "a parent cannot occupy both sides for one child", same))
  dam_used <- unique(parents[, "dam"]); dam_used <- dam_used[dam_used > 0L]
  sire_used <- unique(parents[, "sire"]); sire_used <- sire_used[sire_used > 0L]
  mixed <- intersect(dam_used, sire_used)
  if (length(mixed)) return(invalid(
    "inferred_role_conflict", "an individual cannot serve globally as both dam and sire",
    mixed))
  sex <- state$sex
  bad_dam <- dam_used[!is.na(sex[dam_used]) & sex[dam_used] == 2L]
  bad_sire <- sire_used[!is.na(sex[sire_used]) & sex[sire_used] == 1L]
  if (length(bad_dam) || length(bad_sire)) return(invalid(
    "known_sex_role_conflict", "known sex conflicts with parent role",
    c(bad_dam, bad_sire)))
  if (.parentage_cycle(parents))
    return(invalid("ancestry_cycle", "ancestry cycle detected", seq_along(ids)))
  for (i in seq_len(nrow(parents))) for (j in seq_len(2L)) {
    parent <- parents[i, j]
    if (parent == 0L) next
    side <- if (j == 1L) "M" else "P"
    child_id <- ids[i]; parent_id <- ids[parent]
    if (!.parent_age_feasible(state$life_history, child_id, parent_id, side,
                              state$min_age_parent, state$max_age_parent))
      return(invalid("impossible_age_direction",
                     "parent and offspring birth-year supports violate age bounds",
                     c(i, parent)))
    if (!.parent_year_last_feasible(state$life_history, child_id, parent_id))
      return(invalid("year_last_conflict",
                     "offspring birth-year support occurs after parent Year.last",
                     c(i, parent)))
    if (!.parent_age_year_feasible(state$life_history, child_id, parent_id, side,
                                   state$min_age_parent, state$max_age_parent))
      return(invalid("impossible_age_year_support",
                     "no joint birth-year support satisfies age and Year.last",
                     c(i, parent)))
  }
  lh <- state$life_history
  age <- typed_age_query_cpp(typed_age_state_cpp(parents, as.integer(lh$BY.min),
    as.integer(lh$BY.max), as.integer(lh$Year.last),
    as.numeric(state$min_age_parent), as.numeric(state$max_age_parent)), 0L, integer())
  if (!age$feasible) return(invalid("impossible_age_component",
    "no jointly feasible typed birth years under the assigned links", seq_along(ids)))
  .consistency_result(TRUE)
}

.validate_candidate_backed_parents <- function(state, parents = state$parents,
                                               candidate_pool = state$candidate_pool,
                                               ctx = state$ctx,
                                               edges_to_validate = NULL,
                                               pair_children = NULL,
                                               scored_proposal = NULL) {
  if (is.null(candidate_pool)) return(.consistency_result(TRUE))
  invalid <- function(code, reason, ids)
    .consistency_result(FALSE, code, reason, ids)
  if (!inherits(candidate_pool, "PaPeR_candidate_pool"))
    return(invalid("invalid_candidate_pool", "candidate pool has the wrong type", character()))
  if (!inherits(ctx, "PaPeR_inference_context"))
    return(invalid("missing_candidate_context",
                   "candidate-backed validation requires an inference context", character()))
  ids <- .parentage_state_ids(state)
  edges <- .parentage_edges_from_matrix(ids, parents)
  if (!is.null(edges_to_validate)) {
    wanted <- if (!nrow(edges_to_validate)) character() else paste(
      edges_to_validate$child, edges_to_validate$side, edges_to_validate$parent,
      sep = "\r")
    keys <- if (!nrow(edges)) character() else paste(
      edges$child, edges$side, edges$parent, sep = "\r")
    edges <- edges[keys %in% wanted, , drop = FALSE]
  }
  for (i in seq_len(nrow(edges))) {
    row <- which(candidate_pool$rows$offspring == edges$child[i] &
                 candidate_pool$rows$side == edges$side[i] &
                 candidate_pool$rows$candidate == edges$parent[i])
    pair_ids <- c(edges$child[i], edges$parent[i])
    if (length(row) != 1L)
      return(invalid("candidate_absent", "supplied parent is absent from the retained universe",
                     pair_ids))
    row <- candidate_pool$rows[row, , drop = FALSE]
    if (!isTRUE(row$oh_pass))
      return(invalid("oh_threshold", "supplied parent fails the permanent OH threshold",
                     pair_ids))
    if (!isTRUE(row$permanent_pass))
      return(invalid("permanent_candidate_ineligible",
                     paste0("supplied parent fails permanent eligibility: ",
                            row$permanent_failure_reason), pair_ids))
  }
  pair_rows <- which(parents[, "dam"] > 0L & parents[, "sire"] > 0L)
  if (!is.null(pair_children)) pair_rows <- intersect(
    pair_rows, match(as.character(pair_children), ids, nomatch = 0L))
  for (child in pair_rows) {
    me <- count_mendelian_errors_batch_cpp(
      ctx$G, state$nodes$genotype_row[child] - 1L,
      state$nodes$genotype_row[parents[child, "dam"]] - 1L,
      state$nodes$genotype_row[parents[child, "sire"]] - 1L)
    if (as.integer(me[1L]) > candidate_pool$MaxMendelE)
      return(invalid("mendelian_error_threshold",
                     "supplied parent pair fails the permanent ME threshold",
                     ids[c(child, parents[child, ])]))
  }
  .consistency_result(TRUE)
}

.parentage_edges_from_matrix <- function(ids, parents) {
  out <- lapply(seq_len(2L), function(j) {
    keep <- parents[, j] > 0L
    if (!any(keep)) return(NULL)
    data.frame(child = ids[keep], side = c("dam", "sire")[j],
               parent = ids[parents[keep, j]], stringsAsFactors = FALSE)
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(.empty_edges())
  out <- do.call(rbind, out)
  out <- out[order(out$child, match(out$side, c("dam", "sire")), out$parent),
             , drop = FALSE]
  rownames(out) <- NULL
  out
}

.initialize_parentage_state <- function(ids, genotype_rows = NULL,
                                        life_history = NULL,
                                        min_age_parent = 1L,
                                        max_age_parent = NULL,
                                        pedigree = NULL, ctx = NULL,
                                        candidate_pool = NULL) {
  ids <- as.character(ids)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("state IDs must be unique nonempty strings", call. = FALSE)
  if (is.null(genotype_rows)) genotype_rows <- seq_along(ids)
  genotype_rows <- as.integer(genotype_rows)
  if (length(genotype_rows) != length(ids) || anyNA(genotype_rows) ||
      any(genotype_rows < 1L) || anyDuplicated(genotype_rows))
    stop("genotype_rows must be unique positive rows, one per ID", call. = FALSE)
  if (is.null(life_history)) life_history <- data.frame(
    ID = ids, Sex = NA_integer_, BirthYear = NA_integer_, stringsAsFactors = FALSE)
  lh <- .normalize_life_history(life_history, ids)
  bounds <- .parent_bounds(lh, min_age_parent, max_age_parent)
  state <- .new_parentage_state_empty(
    ids, genotype_rows, lh, bounds$min, bounds$max, ctx)
  state$parents[,] <- .normalize_supplied_pedigree(pedigree, ids)
  structural <- .check_parent_matrix(state, state$parents)
  if (!structural$valid)
    stop("supplied pedigree structural validation failed [", structural$code,
         "]: ", structural$reason, call. = FALSE)
  indexes <- .derive_parentage_indexes(ids, state$parents, state$sex)
  state$reverse_offspring <- indexes$reverse_offspring
  state$inferred_role <- indexes$inferred_role
  state$role_counts <- indexes$role_counts
  if (!is.null(candidate_pool)) .attach_parentage_candidate_pool(state, candidate_pool, ctx)
  state
}

# Compatibility constructor used by finalized scorers and older focused tests.
# It creates the same authoritative ParentageState and does not store edges.

.attach_parentage_candidate_pool <- function(state, candidate_pool, ctx = state$ctx) {
  check <- .validate_candidate_backed_parents(state, state$parents,
                                              candidate_pool, ctx)
  if (!check$valid)
    stop("supplied pedigree candidate validation failed [", check$code,
         "]: ", check$reason, call. = FALSE)
  state$candidate_pool <- candidate_pool
  state$ctx <- ctx
  invisible(state)
}


.patch_staleness <- function(state, patch) {
  if (!inherits(patch, "PaPeR_ParentagePatch"))
    return(.consistency_result(FALSE, "invalid_patch",
                               "patch must be a ParentagePatch"))
  if (!identical(as.integer(patch$expected_graph_revision),
                 as.integer(state$graph_revision)))
    return(.consistency_result(FALSE, "stale_graph_revision",
                               "patch was created under a different graph revision",
                               patch$affected_nodes))
  if (!identical(as.integer(patch$expected_message_revision),
                 as.integer(state$message_revision)))
    return(.consistency_result(FALSE, "stale_message_revision",
                               "patch was scored under a different message revision",
                               patch$affected_nodes))
  expected <- patch$expected_node_revisions
  if (is.null(names(expected)) || any(!names(expected) %in% names(state$node_revision)) ||
      any(state$node_revision[names(expected)] != expected))
    return(.consistency_result(FALSE, "stale_node_revision",
                               "patch has stale affected-node revisions",
                               names(expected)))
  .consistency_result(TRUE)
}

.check_pedigree_consistency <- function(state, patch = NULL) {
  if (!inherits(state, "PaPeR_PedigreeState"))
    return(.consistency_result(FALSE, "invalid_state",
                               "state must be a typed ParentageState"))
  if (is.null(patch)) return(.check_parent_matrix(state, state$parents))
  stale <- .patch_staleness(state, patch)
  if (!stale$valid) return(stale)
  if (identical(patch$move_type, "no_change"))
    return(.consistency_result(FALSE, "no_change", "patch does not change parentage",
                               patch$child_id))
  current <- state$parents[patch$child, ]
  if (!identical(as.integer(current),
                 as.integer(c(patch$old_dam, patch$old_sire))))
    return(.consistency_result(FALSE, "stale_parent_configuration",
                               "patch old configuration no longer matches state",
                               patch$child_id))
  preview <- state$parents
  preview[patch$child, ] <- c(patch$proposed_dam, patch$proposed_sire)
  structural <- .check_parent_matrix(state, preview)
  if (!structural$valid) return(structural)
  .validate_candidate_backed_parents(
    state, preview, edges_to_validate = patch$edges_added,
    pair_children = patch$child_id, scored_proposal = patch$scored_proposal)
}

.validate_state_graph <- function(state) {
  result <- .check_pedigree_consistency(state)
  if (!result$valid) stop(result$reason, call. = FALSE)
  invisible(TRUE)
}
