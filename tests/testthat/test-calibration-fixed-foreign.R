test_that('Foreign is profiled but never adjusted under any calibration policy', {
  for (foreign_bounds in list(c(.004,.008), c(.00001,.0001))) {
    fit <- list(status='fitted', estimated_rates=c(Dropout=.04,Foreign=.006),
      interval=rbind(Dropout=c(lower=.03,upper=.05),
        Foreign=setNames(foreign_bounds,c('lower','upper'))))
    current <- c(Dropout=.02,Foreign=.001)
    for (direction in c('raise_only','lower_only','two_sided')) {
      x <- .parentage_error_bound_decision(fit,current,direction)
      expect_identical(unname(x$used_rate['Foreign']),.001)
      expect_false(x$parameters$Foreign$adjust)
      expect_equal(x$parameters$Foreign$lower_95,foreign_bounds[1])
    }
    x <- .parentage_error_raise_guarded_step(fit,current,100,100,1)
    expect_equal(x$next_rates,c(Dropout=.03,Foreign=.001))
    fit$interval['Dropout',] <- c(.003,.005)
    x <- .parentage_error_raise_guarded_step(fit,current,100,100,1)
    expect_equal(x$next_rates,c(Dropout=.005,Foreign=.001))
    x <- .parentage_error_raise_guarded_step(fit,current,100,100,2)
    expect_false(x$continue)
    expect_identical(x$next_rates,current)
  }
})
