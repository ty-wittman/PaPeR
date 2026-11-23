.consistency_state <- function(ids = c("p", "q", "r", "s", "a", "b", "c"),
                               sex = rep(NA_integer_, length(ids)), dam, sire) {
  PaPeR:::.initialize_parentage_state(ids,
    life_history = data.frame(ID = ids, Sex = sex, BirthYear = NA_integer_),
    pedigree = data.frame(id = ids, dam = dam, sire = sire))
}

test_that("one integer matrix is the authoritative typed pedigree", {
  state <- .consistency_state(dam = c(NA, NA, NA, NA, "p", NA, NA),
                              sire = c(NA, NA, NA, NA, "r", NA, NA))
  expect_s3_class(state, "PaPeR_ParentageState")
  expect_type(state$parents, "integer")
  expect_identical(colnames(state$parents), c("dam", "sire"))
  expect_null(state$edges)
  derived <- PaPeR:::.parentage_state_data_frame(state)
  changed <- derived
  changed$dam[] <- "q"
  expect_identical(PaPeR:::.parentage_state_data_frame(state), derived)
  expect_true(all(state$nodes$node_type == "individual"))
  expect_true(all(c("graph_revision", "node_revision", "component_revision",
                    "message_revision", "current_messages", "node_likelihood_ln",
                    "total_likelihood_ln", "candidate_pool",
                    "unresolved_diagnostics", "accepted_move_history",
                    "round_history") %in% ls(state)))
})
