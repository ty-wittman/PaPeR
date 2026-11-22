test_that("simp has one canonical ordinary relationship universe", {
  parent <- PaPeR:::.relationship_universe("simp", "PO", operation = "parent")
  pair <- PaPeR:::.relationship_universe("simp", "RELATED", operation = "pair")
  expect_identical(parent$focal_configurations, "PO")
  expect_setequal(parent$alternative_configurations,
                  c("U", "FS", "HS", "GP", "FA", "HA"))
  expect_setequal(pair$focal_configurations,
                  c("PO", "FS", "HS", "GP", "FA", "HA"))
  expect_identical(parent$duplicate_guards, character())
  expect_identical(unname(parent$routing["PO"]), "PO")
  expect_identical(unname(parent$orientation["PO"]), "forward")
  full <- PaPeR:::.relationship_universe(
    "full", "PO", alternatives = "U", operation = "parent",
    applicability = "conditional")
  expect_identical(full$focal_configurations, c("PO", "HS_PO", "PO_HA"))
  expect_identical(full$alternative_configurations, "U")
})

test_that("aliases and operation-specific routing resolve centrally", {
  parent <- PaPeR:::.relationship_universe(
    "simp", "PO", alternatives = "DEG2", operation = "parent")
  pair <- PaPeR:::.relationship_universe("simp", "HS", operation = "pair")
  aliases <- PaPeR:::.relationship_aliases("simp")
  expect_identical(aliases$DEG2, c("HS", "GP", "FA"))
  expect_setequal(aliases$RELATED, c("PO", "FS", "HS", "GP", "FA", "HA"))
  expect_setequal(parent$alternative_configurations, c("HS", "GP", "FA", "U"))
  expect_identical(unname(pair$orientation["HS"]), "symmetric")
  expect_error(PaPeR:::.relationship_universe(
    "simp", "FS", operation = "parent"), "not a parent")
  expect_error(PaPeR:::.relationship_universe(
    "simp", "PO", alternatives = "invented", operation = "parent"),
    "unknown alternative relationship")
})

test_that("HS GP and FA retain separate labels with their shared coefficient row", {
  universe <- PaPeR:::.relationship_universe("simp", "RELATED", operation = "pair")
  reg <- universe$registry
  expect_true(all(c("HS", "GP", "FA") %in% universe$focal_configurations))
  expect_identical(reg$delta_id[reg$config_id == "HS"],
                   reg$delta_id[reg$config_id == "GP"])
  expect_identical(reg$delta_id[reg$config_id == "HS"],
                   reg$delta_id[reg$config_id == "FA"])
})

test_that("compiled screening uses side-specific winner mapping", {
  x <- rbind(a = rep("1/2", 4L), b = rep("1/1", 4L), c = rep("2/2", 4L))
  colnames(x) <- paste0("L", seq_len(ncol(x)))
  tables <- build_tables(rep(list(c(0.6, 0.4)), ncol(x)), Err = 0.01)
  names(tables) <- colnames(x)
  universe <- PaPeR:::.relationship_universe(
    "simp", "GP", alternatives = c("HS", "FA"), operation = "pair")
  edges <- rel_edges(x, tables, "GP", alts = c("HS", "FA"),
                     thresh = -Inf, oh_max = -1L)
  expect_true(all(edges$focal_winner == "GP"))
  expect_true(all(edges$focal_winner %in% universe$focal_configurations))
  expect_true(all(edges$alternative_winner %in% universe$alternative_configurations))
  expect_true(all(edges$focal_winner_set == "GP"))
  expect_true(any(edges$alternative_winner_set == "HS|FA"))
  expect_true(all(edges$alternative_winner_set %in% c("U", "HS|FA")))

  perm <- c("c", "a", "b")
  edges_perm <- rel_edges(x[perm, , drop = FALSE], tables, "GP",
                          alts = c("HS", "FA"), thresh = -Inf, oh_max = -1L)
  key <- function(z) paste(sort(c(z$id1, z$id2)), collapse = "|")
  expect_identical(
    edges$focal_winner[order(vapply(seq_len(nrow(edges)), function(i) key(edges[i, ]), character(1)))],
    edges_perm$focal_winner[order(vapply(seq_len(nrow(edges_perm)), function(i) key(edges_perm[i, ]), character(1)))])
})

test_that("full static screening excludes contextual rows", {
  universe <- PaPeR:::.relationship_universe("full", "RELATED", operation = "pair")
  expect_true(all(universe$score_configurations %in% universe$static_configurations))
  expect_false(any(c("HS_HA", "HS_PO", "PO_HA") %in% universe$score_configurations))
  expect_error(PaPeR:::.relationship_universe(
    "full", "HS_HA", operation = "pair"), "contextual relationship")
  x <- matrix(c("1/1", "1/2"), nrow = 2L,
              dimnames = list(c("a", "b"), "L1"))
  tables <- build_tables(list(L1 = c(.6, .4)), Err = .01)
  expect_error(rel_edges(x, tables, "HS_HA", Complex = "full"),
    "contextual relationship")
})


