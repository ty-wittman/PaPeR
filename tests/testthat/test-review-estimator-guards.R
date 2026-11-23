test_that('OH ceiling audit uses each called parental side independently', {
  g<-rbind(d=c('1/1','1/1'),s=c('2/2','0'),c=c('1/2','2/2'));colnames(g)<-c('L1','L2')
  t<-build_tables(rep(list(c(.5,.5)),2),Err=.01)
  p<-data.frame(id=c('d','s','c'),dam=c(NA,NA,'d'),sire=c(NA,NA,'s'))
  obs<-.parentage_error_observations(g,t,p)
  expect_equal(obs$n_called,1);expect_equal(unname(obs$dam_oh),1);expect_equal(unname(obs$sire_oh),0)
  fit<-.parentage_error_ceiling_audit(list(observations=obs),1)
  audit<-fit$opposite_homozygote_ceiling
  expect_true(audit$available);expect_false(audit$modeled)
  expect_equal(c(audit$dam_at_T,audit$sire_at_T,audit$any_at_T,audit$any_above_T),c(1,0,1,0))
})

test_that('clean cohorts retain strictly positive fitted rates and interval endpoints', {
  g<-matrix('1/1',12,1,dimnames=list(paste0('i',1:12),'L1'))
  p<-matrix(0L,12,2);p[3:12,]<-rep(c(1L,2L),each=10)
  for(t in list(build_tables(list(c(.5,.5)),Err=.001),
    build_tables(list(c(.5,.5)),ErrorModel='dropout_foreign',Dropout=.005,Foreign=.0001))){
    fit<-estimate_parentage_error(g,t,p,1L)
    expect_identical(fit$status,'fitted')
    expect_true(all(fit$estimated_rates>0));expect_true(all(fit$interval>0))
    expect_equal(fit$positive_rate_floor,1e-12)
    expect_true(any(fit$at_positive_rate_floor))
    if(!is.null(fit$boundary_fit))expect_true(fit$boundary_fit)
  }
})

