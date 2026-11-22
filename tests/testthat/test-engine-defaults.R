
test_that('width telemetry uses an actual domain and permits a reference override', {
  expect_true(.engine_defaults(list(pedigree_likelihood=TRUE))$likelihood_actual_width_domain)
  expect_false(.engine_defaults(list(pedigree_likelihood=TRUE,likelihood_actual_width_domain=FALSE))$likelihood_actual_width_domain)
})
