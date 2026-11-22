# Compact relationship coefficients and genotype-state enumeration.
# Coefficients are small fixed rows; topology enumeration checks catalog construction.

.gcd_int <- function(a, b) {
  a <- abs(as.integer(a)); b <- abs(as.integer(b))
  while (b != 0L) { z <- a %% b; a <- b; b <- z }
  a
}

.rational_delta_key <- function(numerator, denominator) {
  if (length(numerator) != 9L || length(denominator) != 1L ||
      is.na(denominator) || denominator <= 0L)
    stop("Jacquard numerator/denominator has invalid dimensions", call. = FALSE)
  paste(vapply(as.integer(numerator), function(x) {
    if (x == 0L) return("0/1")
    g <- .gcd_int(x, denominator)
    paste0(x %/% g, "/", denominator %/% g)
  }, character(1)), collapse = ",")
}

.catalog_topology <- function(founders, parents, focal)
  list(founders = as.character(founders), parents = parents,
       focal = as.character(focal))

.catalog_structural_node_colours <- function(topology) {
  checked <- .validate_complex_topology(topology)
  nodes <- checked$nodes
  initial <- setNames(vapply(nodes, function(node) paste0(
    if (node %in% topology$founders) "F" else "N", ":",
    paste(which(topology$focal == node), collapse = ",")), character(1)), nodes)
  rank_values <- function(x) match(x, sort(unique(x), method = "radix"))
  colours <- setNames(rank_values(initial), nodes)
  descriptors <- initial
  for (iteration in seq_len(length(nodes) + 2L)) {
    descriptors <- setNames(vapply(nodes, function(node) {
      parent_colours <- if (node %in% names(topology$parents))
        colours[topology$parents[[node]]] else integer()
      child_roles <- character()
      for (child in names(topology$parents)) for (side in seq_len(2L))
        if (identical(topology$parents[[child]][side], node))
          child_roles <- c(child_roles, paste(colours[[child]], side, sep = ":"))
      paste(initial[[node]], colours[[node]],
            paste(parent_colours, collapse = ","),
            paste(sort(child_roles, method = "radix"), collapse = ","), sep = "|")
    }, character(1)), nodes)
    next_colours <- setNames(rank_values(descriptors), nodes)
    if (identical(unname(next_colours), unname(colours))) break
    colours <- next_colours
  }
  list(colours = colours, descriptors = descriptors)
}

.catalog_registry_topology_signature <- function(topology) {
  checked <- .validate_complex_topology(topology)
  refined <- .catalog_structural_node_colours(topology)
  edge_key <- if (length(topology$parents)) paste(sort(vapply(
    names(topology$parents), function(child) paste(
      refined$colours[[child]],
      paste(refined$colours[topology$parents[[child]]], collapse = ","),
      sep = "<-"), character(1)), method = "radix"), collapse = ";") else ""
  paste(paste(sort(refined$colours[topology$founders]), collapse = ","),
        edge_key, paste(refined$colours[topology$focal], collapse = ","),
        sep = "|")
}

# These are the finite typed-parent configurations represented by the active
# Complex = "full" registry (see docs/MODEL.md, section 6).  No dummy node, sibship cluster, or reconstruction
# topology is included here.
.catalog_topology_catalog <- function() {
  T <- .catalog_topology
  list(
    U = T(c("A", "B"), list(), c("A", "B")),
    PO = T(c("B", "Q"), list(A = c("B", "Q")), c("A", "B")),
    FS = T(c("P", "Q"), list(A = c("P", "Q"), B = c("P", "Q")),
           c("A", "B")),
    HS = T(c("P", "Q", "R"),
           list(A = c("P", "Q"), B = c("P", "R")), c("A", "B")),
    GP = T(c("B", "Q", "R"),
           list(P = c("B", "Q"), A = c("P", "R")), c("A", "B")),
    FA = T(c("P", "Q", "R"),
           list(B = c("P", "Q"), S = c("P", "Q"), A = c("S", "R")),
           c("A", "B")),
    HA = T(c("P", "Q", "R", "S"),
           list(B = c("P", "Q"), T = c("P", "R"), A = c("T", "S")),
           c("A", "B")),
    MZ = T("A", list(), c("A", "A")),

    # Typed-parent double/inbred routes used by current Sequoia's Complex=2
    # branches.  The topology, rather than a label alone, is the mathematical
    # definition of each row.
    HS_PO = T(c("P", "Q"),
              list(B = c("P", "Q"), A = c("B", "P")), c("A", "B")),
    PO_HA = T(c("R", "S", "T"),
              list(P = c("R", "S"), Q = c("R", "T"),
                   A = c("P", "Q")), c("A", "P")),
    HS_HA = T(c("P", "R", "S"),
              list(B = c("P", "R"), Q = c("R", "S"),
                   A = c("P", "Q")), c("A", "B")),
    HS_F1C = T(c("P", "G", "H"),
               list(Q = c("G", "H"), R = c("G", "H"),
                    A = c("P", "Q"), B = c("P", "R")), c("A", "B")),
    HS_FA = T(c("P", "R"),
              list(B = c("P", "R"), Q = c("P", "R"),
                   A = c("P", "Q")), c("A", "B")),
    HS_GP = T(c("P", "S", "R"),
              list(B = c("P", "S"), Q = c("B", "R"),
                   A = c("P", "Q")), c("A", "B")),
    HS_GGG = T(c("P", "S", "R", "T"),
               list(B = c("P", "S"), C = c("B", "R"),
                    Q = c("C", "T"), A = c("P", "Q")), c("A", "B")),
    FA_HA = T(c("G", "H", "J"),
              list(B = c("G", "H"), P = c("G", "H"),
                   Q = c("G", "J"), A = c("P", "Q")), c("A", "B")),
    FS_HA = T(c("P", "R"),
              list(Q = c("P", "R"), A = c("P", "Q"),
                   B = c("P", "Q")), c("A", "B")),
    FS_H1C = T(c("G", "R", "S"),
               list(P = c("G", "R"), Q = c("G", "S"),
                    A = c("P", "Q"), B = c("P", "Q")), c("A", "B")),
    FS_PARENTS_FS = T(c("G", "H"),
                      list(P = c("G", "H"), Q = c("G", "H"),
                           A = c("P", "Q"), B = c("P", "Q")), c("A", "B"))
  )
}

.validate_complex_topology <- function(topology) {
  if (!is.list(topology) ||
      !all(c("founders", "parents", "focal") %in% names(topology)))
    stop("topology must contain founders, parents, and focal", call. = FALSE)
  founders <- topology$founders
  parents <- topology$parents
  focal <- topology$focal
  if (!is.character(founders) || !length(founders) || anyNA(founders) ||
      any(!nzchar(founders)) || anyDuplicated(founders))
    stop("topology founders must be unique nonempty IDs", call. = FALSE)
  if (!is.list(parents) || (length(parents) && is.null(names(parents))))
    stop("topology parents must be a named list", call. = FALSE)
  if (length(parents) && (anyNA(names(parents)) || any(!nzchar(names(parents))) ||
                          anyDuplicated(names(parents))))
    stop("topology nonfounder IDs must be unique nonempty names", call. = FALSE)
  if (length(intersect(founders, names(parents))))
    stop("a topology node cannot be both founder and nonfounder", call. = FALSE)
  if (!is.character(focal) || length(focal) != 2L || anyNA(focal) ||
      any(!nzchar(focal)))
    stop("topology focal must contain two node IDs", call. = FALSE)
  nodes <- c(founders, names(parents))
  if (any(!focal %in% nodes))
    stop("topology focal node is missing from topology", call. = FALSE)
  if (length(parents)) {
    bad_shape <- vapply(parents, function(x)
      !is.character(x) || length(x) != 2L || anyNA(x) || any(!nzchar(x)),
      logical(1))
    if (any(bad_shape))
      stop("every topology nonfounder must have exactly two named parents",
           call. = FALSE)
    missing_parent <- setdiff(unique(unlist(parents, use.names = FALSE)), nodes)
    if (length(missing_parent))
      stop("topology parent references missing parent: ",
           paste(missing_parent, collapse = ", "), call. = FALSE)
    if (any(vapply(names(parents), function(x) x %in% parents[[x]], logical(1))))
      stop("topology self-parentage is not supported; this is an ancestry cycle",
           call. = FALSE)
    if (any(vapply(parents, function(x) x[1L] == x[2L], logical(1))))
      stop("topology selfing is not supported", call. = FALSE)
  }

  done <- founders
  pending <- names(parents)
  order <- character(0)
  while (length(pending)) {
    ready <- pending[vapply(parents[pending], function(x) all(x %in% done),
                            logical(1))]
    if (!length(ready))
      stop("topology contains an ancestry cycle", call. = FALSE)
    order <- c(order, ready)
    done <- c(done, ready)
    pending <- setdiff(pending, ready)
  }
  list(nodes = nodes, order = order)
}

.catalog_registry_rows <- list(
  list("U", 1L, c(0,0,0,0,0,0,0,0,1), "U", "CalcU;CheckPair", "none",
       "none", FALSE, "supported", ""),
  list("PO", 1L, c(0,0,0,0,0,0,0,1,0), "PO", "CalcP2;PairPO;CheckPair",
       "forward", "parent-offspring", FALSE, "supported", ""),
  list("FS", 4L, c(0,0,0,0,0,0,1,2,1), "FS",
       "PairFullSib;CheckPair", "symmetric", "typed full siblings", FALSE,
       "supported", ""),
  list("HS", 2L, c(0,0,0,0,0,0,0,1,1), "HS",
       "PairHalfSib;CheckPair", "symmetric", "typed half siblings", FALSE,
       "supported", ""),
  list("GP", 2L, c(0,0,0,0,0,0,0,1,1), "GP",
       "PairGP;CheckPair", "forward", "typed grandparent path", FALSE,
       "supported", ""),
  list("FA", 2L, c(0,0,0,0,0,0,0,1,1), "FA",
       "PairUA;CheckPair", "signed", "typed full-avuncular path", TRUE,
       "supported", ""),
  list("HA", 4L, c(0,0,0,0,0,0,0,1,3), "HA",
       "PairUA;PairCC;CheckPair", "signed", "typed distant path", TRUE,
       "supported", ""),
  list("MZ", 1L, c(0,0,0,0,0,0,1,0,0), NA_character_,
       "PairSelf;CheckPair", "symmetric", "duplicate guard", FALSE,
       "duplicate_guard", ""),

  list("HS_PO", 4L, c(0,0,1,0,0,0,1,2,0), "PO",
       "PairHSPO;PairPO;CheckPair", "contextual",
       "typed parent-offspring plus half-sib connection", TRUE,
       "supported", "PO_GP|PO_FA"),
  list("PO_HA", 8L, c(0,0,1,0,0,0,1,6,0), "PO",
       "PairPOHA;CheckPair", "forward", "typed parent plus connected HA", FALSE,
       "supported", ""),
  list("HS_HA", 8L, c(0,0,0,0,0,0,1,4,3), "HS",
       "PairHSHA;CheckPair", "symmetric", "typed HS plus HA connection", FALSE,
       "supported", ""),
  list("HS_F1C", 8L, c(0,0,0,0,0,0,1,4,3), "HS",
       "PairHSCC;CheckPair", "symmetric", "typed HS plus full-cousin connection",
       FALSE, "supported", "TQS"),
  list("HS_FA", 8L, c(0,0,1,1,0,0,1,4,1), "HS",
       "pairHSHA;CheckPair", "contextual",
       "typed HS plus FA connection", TRUE, "supported", ""),
  list("HS_GP", 8L, c(0,0,1,0,0,0,1,4,2), "HS",
       "pairHSGP;CheckPair", "contextual",
       "typed HS plus GP connection", TRUE, "supported", ""),
  list("HS_GGG", 16L, c(0,0,1,0,0,0,1,8,6), "HA",
       "PairGGP;PairGA;CheckPair", "signed",
       "typed HS plus great-grandparent connection", TRUE,
       "supported", ""),
  list("FA_HA", 16L, c(0,0,1,1,0,0,1,8,5), "FA",
       "pairFAHA;CheckPair", "signed", "typed FA plus HA connection", TRUE,
       "supported", ""),
  list("FS_HA", 16L, c(1,0,2,1,2,1,3,6,0), "FS",
       "pairFSHA;CheckPair", "symmetric", "typed FS plus HA connection", FALSE,
       "supported", ""),
  list("FS_H1C", 32L, c(1,0,2,1,2,1,7,14,4), "FS",
       "pairFSFC;CheckPair", "symmetric", "typed FS plus HS/FA parents", FALSE,
       "supported", "FS_PARENTS_FA"),
  list("FS_PARENTS_FS", 32L, c(2,1,4,1,4,1,7,10,2), "FS",
       "pairFSFC;CheckPair", "symmetric", "typed FS with FS parents", FALSE,
       "supported", "FSI")
)

.catalog_topology_edges <- function(topology) {
  checked <- .validate_complex_topology(topology)
  if (!length(topology$parents))
    return(data.frame(child = character(), parent = character(),
                      stringsAsFactors = FALSE))
  out <- do.call(rbind, lapply(names(topology$parents), function(child) {
    data.frame(child = rep.int(child, 2L),
               parent = as.character(topology$parents[[child]]),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out[order(match(out$child, checked$order), out$parent), , drop = FALSE]
}

.catalog_registry_applicability <- function(config_id, topology,
                                             age_direction,
                                             both_orientations) {
  ordinary <- config_id %in% c("U", "PO", "FS", "HS", "GP", "FA", "HA", "MZ")
  checked <- .validate_complex_topology(topology)
  edges <- .catalog_topology_edges(topology)
  focal_edge <- function(x) {
    setequal(c(x$child, x$parent), topology$focal)
  }
  contextual_edges <- if (!ordinary && nrow(edges)) {
    if (config_id %in% c("PO", "HS_PO", "PO_HA")) {
      # Every edge into the focal child belongs to the tested parent
      # hypothesis. Required state context is the surrounding HS/HA structure.
      edges[edges$child != topology$focal[1L], , drop = FALSE]
    } else edges[!vapply(seq_len(nrow(edges)), function(i)
      focal_edge(edges[i, , drop = FALSE]), logical(1)), , drop = FALSE]
  }
  else data.frame(child = character(), parent = character(),
                  stringsAsFactors = FALSE)
  required_nodes <- if (ordinary) character() else
    setdiff(checked$nodes, topology$focal)
  allowed <- if (isTRUE(both_orientations)) c("auto", "forward", "reverse")
    else if (identical(age_direction, "symmetric")) c("auto", "symmetric")
    else if (identical(age_direction, "none")) c("auto", "none")
    else c("auto", "forward")
  list(required_nodes = required_nodes,
       required_edges = contextual_edges,
       forbidden_edges = data.frame(child = character(), parent = character(),
                                     stringsAsFactors = FALSE),
       allowed_orientations = allowed,
       missing_context_status = "not_applicable",
       "ordinary_fallback" = ordinary,
       context_required = !ordinary && length(required_nodes) > 0L)
}

.catalog_relationship_registry <- local({
  rows <- .catalog_registry_rows
  ids <- vapply(rows, `[[`, character(1), 1L)
  den <- vapply(rows, `[[`, integer(1), 2L)
  num <- lapply(rows, `[[`, 3L)
  topologies <- .catalog_topology_catalog()
  if (!identical(ids, names(topologies)))
    stop("complex topology catalog and registry are out of sync")
  data <- data.frame(
    config_id = ids,
    delta_key = vapply(seq_along(rows), function(i)
      .rational_delta_key(num[[i]], den[i]), character(1)),
    major_bin = vapply(rows, `[[`, character(1), 4L),
    source_routines = vapply(rows, `[[`, character(1), 5L),
    age_direction = vapply(rows, `[[`, character(1), 6L),
    required_evidence = vapply(rows, `[[`, character(1), 7L),
    both_orientations = vapply(rows, `[[`, logical(1), 8L),
    support_status = vapply(rows, `[[`, character(1), 9L),
    aliases = vapply(rows, `[[`, character(1), 10L),
    stringsAsFactors = FALSE)
  data$delta_id <- paste0("J", match(data$delta_key, unique(data$delta_key)))
  data$Delta <- I(lapply(seq_along(rows), function(i) {
    z <- num[[i]] / den[i]; names(z) <- paste0("D", 1:9); z
  }))
  data$topology <- I(unname(topologies[data$config_id]))
  # Keep the registry row bound to the exact executable topology.  A compact
  # coefficient row is not sufficient identity for a full-search member:
  # callers must be able to validate the topology and its focal slots before
  # scoring.  MZ is intentionally a one-true-individual kernel; its search
  # request supplies the two observed IDs separately (see the Stage-1 full
  # hypothesis layer), rather than treating the repeated focal slot as one
  # observed row.
  data$topology_signature <- vapply(data$topology, .catalog_registry_topology_signature,
                                    character(1))
  data$topology_binding <- I(lapply(seq_len(nrow(data)), function(i) {
    topology <- data$topology[[i]]
    list(configuration = data$config_id[i],
         topology_signature = data$topology_signature[i],
         nodes = .validate_complex_topology(topology)$nodes,
         focal = topology$focal,
         parents = topology$parents,
         duplicate_observation_kernel = identical(data$config_id[i], "MZ"))
  }))
  applicability <- lapply(seq_along(rows), function(i)
    .catalog_registry_applicability(
      data$config_id[i], topologies[[data$config_id[i]]],
      data$age_direction[i], data$both_orientations[i]))
  data$required_nodes <- I(lapply(applicability, `[[`, "required_nodes"))
  data$required_edges <- I(lapply(applicability, `[[`, "required_edges"))
  data$forbidden_edges <- I(lapply(applicability, `[[`, "forbidden_edges"))
  data$allowed_orientations <- I(lapply(applicability, `[[`, "allowed_orientations"))
  data$missing_context_status <- vapply(
    applicability, `[[`, character(1), "missing_context_status")
  data[["ordinary_fallback"]] <- vapply(
    applicability, `[[`, logical(1), "ordinary_fallback")
  data$applicability <- I(applicability)

  function(Complex = c("simp", "full"), include_duplicate = FALSE) {
    Complex <- match.arg(Complex)
    keep <- if (Complex == "simp")
      data$config_id %in% c("U", "PO", "FS", "HS", "GP", "FA", "HA", "MZ")
    else rep(TRUE, nrow(data))
    if (!include_duplicate) keep <- keep & data$support_status != "duplicate_guard"
    out <- data[keep, , drop = FALSE]
    rownames(out) <- NULL
    out
  }
})

.registry_delta_matrix <- function(Complex = c("simp", "full"),
                                   include_duplicate = FALSE) {
  reg <- .catalog_relationship_registry(match.arg(Complex), include_duplicate)
  out <- do.call(rbind, unclass(reg$Delta))
  rownames(out) <- reg$config_id
  colnames(out) <- paste0("D", 1:9)
  out
}

# Immutable, once-per-native-run representation of the reviewed full registry.
# Production C++ receives compact integer rows and Delta coefficients only;
# relationship names remain at this boundary for diagnostics and tests.
.full_native_registry_payload <- function() {
  registry <- .catalog_relationship_registry("full", include_duplicate = FALSE)
  ids <- as.character(registry$config_id)
  orientation_mask <- vapply(registry$allowed_orientations, function(x) {
    # Bit 0 is canonical/symmetric-or-forward; bit 1 is reverse.  The native
    # hot path receives no relationship strings and permutes only reverse rows.
    bitwOr(if ("symmetric" %in% x || "forward" %in% x) 1L else 0L,
           if ("reverse" %in% x) 2L else 0L)
  }, integer(1))
  list(
    schema_version = "stage1-full-native-registry-v1",
    relationship_code = as.integer(seq_along(ids)),
    relationship_id = ids,
    delta = unname(.registry_delta_matrix("full", include_duplicate = FALSE)),
    orientation_mask = as.integer(orientation_mask),
    parent_numerator = ids %in% .full_parent_numerator_owners,
    denominator_alternative = ids %in% .full_relationship_alternatives)
}

.relationship_aliases <- function(Complex = c("simp", "full")) {
  Complex <- match.arg(Complex)
  reg <- .catalog_relationship_registry(Complex, include_duplicate = TRUE)
  out <- setNames(lapply(reg$config_id, function(x) x), reg$config_id)
  for (bin in unique(reg$major_bin[!is.na(reg$major_bin)]))
    out[[bin]] <- unique(reg$config_id[!is.na(reg$major_bin) & reg$major_bin == bin])
  for (i in seq_len(nrow(reg))) {
    aliases <- strsplit(reg$aliases[i], "|", fixed = TRUE)[[1L]]
    aliases <- aliases[nzchar(aliases)]
    for (alias in aliases) out[[alias]] <- unique(c(out[[alias]], reg$config_id[i]))
  }
  out$RELATED <- reg$config_id[
    !is.na(reg$major_bin) & reg$major_bin != "U" &
      reg$support_status != "duplicate_guard"]
  out$DEG2 <- c("HS", "GP", "FA")
  if (Complex == "simp") {
    out$TQS <- "HS"
    out$FSI <- "FS"
  } else {
    out$TQS <- "HS_F1C"
    out$FSI <- "FS_PARENTS_FS"
  }
  out
}


.resolve_relationship_names <- function(names, Complex = c("simp", "full"),
                                        what = "relationship") {
  Complex <- match.arg(Complex)
  if (is.null(names)) return(character(0))
  if (!is.character(names) || !length(names) || anyNA(names) ||
      any(!nzchar(names)))
    stop(what, " names must be nonempty character values", call. = FALSE)
  aliases <- .relationship_aliases(Complex)
  bad <- setdiff(names, names(aliases))
  if (length(bad))
    stop("unknown ", what, " name(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  unique(unlist(aliases[names], use.names = FALSE))
}

.registry_operation_routing <- function(reg, operation) {
  if (operation %in% c("pair", "diagnostic")) return(reg$major_bin)
  if (operation %in% c("parent", "final_llr"))
    return(ifelse(reg$major_bin == "PO", "PO", NA_character_))
  if (operation == "duplicate")
    return(ifelse(reg$config_id == "MZ", "MZ", NA_character_))
  stop("unsupported relationship operation: ", operation, call. = FALSE)
}

.relationship_universe <- function(
    Complex = c("simp", "full"), focal, alternatives = NULL,
    include_duplicate = TRUE,
    applicability = c("unconditional", "conditional"),
    operation = c("pair", "parent", "duplicate", "diagnostic", "final_llr",
                  "trio"), applicable_configurations = NULL) {
  Complex <- match.arg(Complex)
  applicability <- match.arg(applicability)
  operation <- match.arg(operation)

  if (operation == "trio") {
    trio <- c("H_PP", "H_D0", "H_0S", "H_U")
    if (!is.character(focal) || length(focal) != 1L || is.na(focal) ||
        !focal %in% trio)
      stop("trio focal relationship must be one of: ", paste(trio, collapse = ", "),
           call. = FALSE)
    if (is.null(alternatives)) alternatives <- setdiff(trio, focal)
    if (!is.character(alternatives) || !length(alternatives) || anyNA(alternatives) ||
        any(!alternatives %in% trio))
      stop("trio alternatives must be H_PP, H_D0, H_0S, or H_U", call. = FALSE)
    alternatives <- unique(setdiff(alternatives, focal))
    if (!length(alternatives)) stop("trio relationship universe has no alternatives",
                                    call. = FALSE)
    return(structure(list(
      Complex = Complex, operation = operation, applicability = applicability,
      score_kind = "trio", registry = data.frame(config_id = trio),
      routing = setNames(trio, trio), orientation = setNames(rep("none", 4L), trio),
      focal_configurations = focal, alternative_configurations = alternatives,
      duplicate_guards = character(), alias_expansion = setNames(as.list(trio), trio),
      duplicate_guard_indices = integer(), duplicate_guard_threshold_ln = NA_real_,
      duplicate_guard_model = NA_character_,
      static_configurations = trio, contextual_configurations = character(),
      score_delta = NULL, score_configurations = trio,
      focal_indices = match(focal, trio), alternative_indices = match(alternatives, trio)),
      class = "PaPeR_relationship_universe"))
  }

  reg <- .catalog_relationship_registry(Complex, include_duplicate = TRUE)
  routing <- .registry_operation_routing(reg, operation)
  if (!is.character(focal) || !length(focal) || anyNA(focal) || any(!nzchar(focal)))
    stop("focal relationship names must be nonempty character values", call. = FALSE)

  eligible <- if (applicability == "conditional")
    reg$support_status == "supported" | reg$config_id == "U"
  else
    reg[["ordinary_fallback"]] | reg$support_status == "duplicate_guard"
  if (!is.null(applicable_configurations)) {
    selected <- .resolve_relationship_names(applicable_configurations, Complex,
                                             "applicable relationship")
    eligible <- eligible & reg$config_id %in% selected
  }

  focal_named <- if (identical(focal, "RELATED")) character() else
    .resolve_relationship_names(focal, Complex, "focal relationship")
  unavailable_focal <- focal_named[!eligible[match(focal_named, reg$config_id)]]
  relationship_labels <- c("RELATED", "U", "PO", "FS", "HS", "GP", "FA", "HA")
  if (applicability == "unconditional" && length(unavailable_focal) &&
      any(unavailable_focal %in% reg$config_id[!reg[["ordinary_fallback"]]]) &&
      any(!focal %in% relationship_labels))
    stop("contextual relationship requires typed evidence and is not applicable " ,
         "to unconditional static screening: ", paste(unavailable_focal, collapse = ", "),
         call. = FALSE)

  if (identical(focal, "RELATED")) {
    focal_configs <- reg$config_id[eligible & !is.na(routing) & routing != "U" &
                                    reg$support_status != "duplicate_guard"]
  } else {
    direct_allowed <- if (operation %in% c("parent", "final_llr", "duplicate"))
      !is.na(routing) else rep(TRUE, nrow(reg))
    focal_configs <- reg$config_id[eligible &
      ((reg$config_id %in% focal_named & direct_allowed) | routing %in% focal_named |
       routing %in% focal[focal %in% c("U", "PO", "FS", "HS", "GP", "FA", "HA")])]
  }
  if (!include_duplicate) focal_configs <- setdiff(focal_configs, "MZ")
  if (!length(focal_configs)) {
    if (operation %in% c("parent", "final_llr"))
      stop("focal relationship is not a parent relationship for this operation",
           call. = FALSE)
    stop("focal relationship is unavailable", call. = FALSE)
  }

  if (operation == "duplicate") {
    if (!identical(focal_configs, "MZ"))
      stop("duplicate operation requires focal relationship MZ", call. = FALSE)
    alternative_configs <- "U"
  } else if (is.null(alternatives)) {
    alternative_configs <- reg$config_id[eligible & !reg$config_id %in% focal_configs]
  } else {
    selected <- .resolve_relationship_names(alternatives, Complex,
                                             "alternative relationship")
    unavailable <- selected[!eligible[match(selected, reg$config_id)]]
    if (applicability == "unconditional" && length(unavailable) &&
        any(unavailable %in% reg$config_id[!reg[["ordinary_fallback"]]]) &&
        any(!alternatives %in% relationship_labels))
      stop("contextual relationship requires typed evidence and is not applicable " ,
           "to unconditional static screening: ", paste(unavailable, collapse = ", "),
           call. = FALSE)
    alternative_configs <- reg$config_id[eligible &
      (reg$config_id %in% selected | routing %in% selected |
       routing %in% alternatives[alternatives %in% c("U", "PO", "FS", "HS", "GP", "FA", "HA")])]
    alternative_configs <- setdiff(alternative_configs, focal_configs)
  }
  alternative_configs <- unique(c(alternative_configs,
    if (operation != "duplicate" && !"U" %in% focal_configs &&
        eligible[match("U", reg$config_id)]) "U"))
  # MZ is available for explicit duplicate reports, not ordinary inference.
  alternative_configs <- setdiff(alternative_configs, "MZ")
  alternative_configs <- intersect(alternative_configs, reg$config_id[eligible])
  if (!length(alternative_configs)) stop("relationship universe has no alternatives",
                                         call. = FALSE)

  used <- unique(c(focal_configs, alternative_configs))
  used_rows <- reg[match(used, reg$config_id), , drop = FALSE]
  delta_ids <- unique(used_rows$delta_id)
  score_rows <- used_rows[match(delta_ids, used_rows$delta_id), , drop = FALSE]
  delta <- do.call(rbind, unclass(score_rows$Delta))
  rownames(delta) <- score_rows$delta_id
  colnames(delta) <- paste0("D", 1:9)
  score_index <- setNames(match(reg$delta_id, score_rows$delta_id), reg$config_id)
  structure(list(
    Complex = Complex, operation = operation, applicability = applicability,
    score_kind = "pair", registry = reg, routing = setNames(routing, reg$config_id),
    orientation = setNames(reg$age_direction, reg$config_id),
    focal_configurations = focal_configs,
    alternative_configurations = alternative_configs,
    duplicate_guards = character(), duplicate_guard_indices = integer(),
    duplicate_guard_threshold_ln = NA_real_, duplicate_guard_model = NA_character_,
    static_configurations = reg$config_id[reg[["ordinary_fallback"]]],
    contextual_configurations = reg$config_id[!reg[["ordinary_fallback"]] &
                                                reg$support_status == "supported"],
    alias_expansion = .relationship_aliases(Complex),
    score_delta = delta, score_configurations = score_rows$config_id,
    focal_indices = unname(unique(score_index[focal_configs])),
    alternative_indices = unname(unique(score_index[alternative_configs]))),
    class = "PaPeR_relationship_universe")
}

# Finalized parentage uses only the routes that own the requested output in
# the source operation.  This table is deliberately separate from the global
# registry: a registered relationship that is not in an operation's family is
# excluded, not an unsupported competitor.  MZ is handled by the duplicate
# operation and never enters a biological parent contrast.


.relationship_winner_labels <- function(indices, universe, side) {
  if (!inherits(universe, "PaPeR_relationship_universe") ||
      !side %in% c("focal", "alternative"))
    stop("side must be focal or alternative for relationship winner mapping",
         call. = FALSE)
  configurations <- if (side == "focal") universe$focal_configurations else
    universe$alternative_configurations
  score_rows <- universe$score_configurations
  row_keys <- universe$registry$delta_id[match(score_rows, universe$registry$config_id)]
  config_keys <- universe$registry$delta_id[
    match(configurations, universe$registry$config_id)]
  label_sets <- lapply(indices, function(index) {
    if (is.na(index) || index < 1L || index > length(row_keys)) return(character())
    configurations[config_keys == row_keys[index]]
  })
  winners <- vapply(label_sets, function(x) if (length(x)) x[1L] else NA_character_, character(1))
  sets <- vapply(label_sets, function(x) paste(x, collapse = "|"), character(1))
  list(winner = winners, set = sets)
}

# -------------------------------------------------------------------------
# Stage-1 Complex="full" search request contract
# -------------------------------------------------------------------------
#
# The finalized scorer above predates the explicit full-search request.  Keep
# it intact for the frozen output-only paths and expose a small, pure request
# layer for Stage-1 search.  These helpers intentionally contain no chooser or
# native-engine calls.  A request is a value object: callers may retain it as
# provenance, but scoring never writes to it or to the supplied state.

.full_parent_numerator_owners <- c("PO", "HS_PO", "PO_HA")
.full_relationship_alternatives <- c(
  "U", "FS", "HS", "GP", "FA", "HA", "HS_HA", "HS_F1C",
  "HS_FA", "HS_GP", "HS_GGG", "FA_HA", "FS_HA", "FS_H1C",
  "FS_PARENTS_FS")
