test_that("truncated count probabilities preserve missingness and selection mass", {
  p <- c(.1,.3,.4)
  called <- rbind(c(TRUE,TRUE,TRUE),c(TRUE,FALSE,TRUE))
  pmf <- .parentage_error_count_pmf(p,called,1L)
  reference <- function(q)c(prod(1-q),sum(vapply(seq_along(q),function(j)
    q[j]*prod(1-q[-j]),0.0)))
  expect_equal(pmf[1,],reference(p))
  expect_equal(pmf[2,],reference(p[c(1,3)]))
  expect_true(all(rowSums(pmf)<1))
  expect_equal(rowSums(pmf/rowSums(pmf)),c(1,1))
})

test_that("ceiling correction recovers a separately constructed truncated population", {
  L<-30L;N<-2000L;ceiling<-2L;rate<-.04
  tables<-build_tables(rep(list(c(.5,.5)),L),Err=.001)
  true_table<-build_tables(list(c(.5,.5)),Err=rate)[[1]]
  p<-true_table$p_me_true_trio
  weights<-stats::dbinom(0:ceiling,L,p);weights<-weights/sum(weights)
  counts<-findInterval((seq_len(N)-.5)/N,cumsum(weights))
  g<-matrix('1/1',N+2L,L,dimnames=list(paste0('i',seq_len(N+2L)),paste0('l',seq_len(L))))
  for(i in seq_len(N))if(counts[i])g[i+2L,seq_len(counts[i])]<-'1/2'
  parents<-matrix(0L,N+2L,2);parents[3L:(N+2L),]<-rep(c(1L,2L),each=N)
  fit<-estimate_parentage_error(g,tables,parents,ceiling)
  expect_equal(fit$status,'fitted')
  expect_lt(abs(fit$estimate-rate),1e-4)
  expect_lt(fit$naive_estimate,fit$estimate-.005)
  expect_true(fit$profile_excludes_initial)
  expect_equal(fit$observations$count,counts)
  expect_equal(sum(fit$histogram$observed),N)
  expect_equal(sum(fit$histogram$fitted_given_ceiling),N)
})

test_that("calibration preserves profiles and refuses unidentifiable cohorts", {
  tables<-build_tables(list(c(.7,.3)),ErrorModel='dropout_foreign',Dropout=.005,Foreign=.0001)
  # MH calibration profiles two rates; the historical common-scale adapter
  # is no longer in the estimator or assignment path.
  rebuilt <- .parentage_mh_rebuild(tables,c(Dropout=.02,Foreign=.0004))
  expect_equal(rebuilt[[1]]$Dropout,.02)
  expect_equal(rebuilt[[1]]$Foreign,.0004)
  expect_equal(rebuilt[[1]]$ForeignQ,tables[[1]]$ForeignQ)
  expect_equal(rebuilt[[1]]$ForeignQMode,tables[[1]]$ForeignQMode)
  g<-matrix('1/1',3L,1L,dimnames=list(c('d','s','c'),'l'))
  parents<-rbind(c(0L,0L),c(0L,0L),c(1L,2L))
  expect_equal(estimate_parentage_error(g,tables,parents,1)$status,'too_few_trios')
  expect_equal(estimate_parentage_error(g,tables,parents,0)$status,'unidentifiable_zero_ceiling')
  expect_error(estimate_parentage_error(g,tables,parents,1,confident=c(TRUE,NA,TRUE)),'confident')
  expect_error(estimate_parentage_error(g,tables,parents,-1),'MaxMendelE')
})

