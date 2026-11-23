
test_that('iteration tracks directions by parameter and stops before a reversal', {
 fit<-function(lo,hi)list(status='fitted',interval=rbind(Dropout=c(lower=lo[1],upper=hi[1]),Foreign=c(lower=lo[2],upper=hi[2])))
 current<-c(Dropout=.01,Foreign=.01);previous<-c(Dropout=0,Foreign=0)
 a<-.parentage_error_iteration_step(fit(c(.02,.001),c(.04,.005)),current,previous,1L)
 expect_true(a$continue);expect_equal(a$next_rates,c(Dropout=.02,Foreign=.01))
 expect_equal(a$previous_direction,c(Dropout=1,Foreign=0))
 b<-.parentage_error_iteration_step(fit(c(.001,.001),c(.01,.003)),a$next_rates,a$previous_direction,2L)
 expect_false(b$continue);expect_identical(b$stop_reason,'direction_reversal')
 expect_identical(b$reversal_parameters,'Dropout');expect_identical(b$next_rates,a$next_rates)
 retained<-.parentage_error_iteration_step(fit(c(.001,.001),c(.04,.008)),a$next_rates,a$previous_direction,2L)
 expect_identical(retained$stop_reason,'supplied_rate_retained')
 expect_equal(retained$previous_direction,a$previous_direction)
})

