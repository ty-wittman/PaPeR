test_that('ME categories match an independent full observed-trio enumeration', {
 p<-c(.3,.7);t<-build_locus_table(p,ErrorModel='dropout_foreign',Dropout=.07,Foreign=.013)
 g<-t$geno;G<-nrow(g);prior<-p[g[,1]]*p[g[,2]]*ifelse(g[,1]==g[,2],1,2)
 # Independent generative contraction: enumerate true parents, four gamete
 # choices, and every observed trio. No containment identities are used.
 E<-matrix(emission_probabilities_cpp(t,rep(0:(G-1),G),rep(0:(G-1),each=G)),G,G)
 mass<-array(0,c(G,G,G))
 index<-function(a,b)which(g[,1]==min(a,b)&g[,2]==max(a,b))
 for(d in seq_len(G))for(s in seq_len(G))for(a in 1:2)for(b in 1:2){
  c<-index(g[d,a],g[s,b])
  for(od in seq_len(G))for(os in seq_len(G))
   mass[,od,os]<-mass[,od,os]+prior[d]*prior[s]/4*E[d,od]*E[s,os]*E[c,]
 }
 expected<-c(0,0)
 for(c in seq_len(G))for(d in seq_len(G))for(s in seq_len(G)){
  compatible<-(g[c,1]%in%g[d,]&g[c,2]%in%g[s,])|(g[c,2]%in%g[d,]&g[c,1]%in%g[s,])
  if(!compatible){j<-if(g[c,1]==g[c,2])1L else 2L;expected[j]<-expected[j]+mass[c,d,s]}
 }
 expect_equal(sum(mass),1,tolerance=1e-12)
 expect_equal(unname(.parentage_mh_me_categories(t)),expected,tolerance=1e-12)
 expect_equal(sum(expected),t$p_me_true_trio,tolerance=1e-12)
})

test_that('compiled two-rate probabilities reproduce independent direct points', {
 for(p in list(c(.3,.7),c(0,.2,.3,.5),rep(1/8,8))){
  t<-build_locus_table(p,ErrorModel='dropout_foreign',Dropout=.005,Foreign=.0001)
  model<-.parentage_mh_category_model(list(t))
  for(rate in list(c(1e-12,1e-12),c(.04,.002),c(.002,.08),c(.25,.25),c(.30,.015))){
   rebuilt<-.parentage_mh_rebuild(list(t),setNames(rate,c('Dropout','Foreign')))[[1]]
   expect_equal(as.vector(model(rate)),unname(.parentage_mh_me_categories(rebuilt)),tolerance=1e-10)
  }
 }
})

test_that('Dropout uses its own profile bound while Foreign stays supplied', {
 fit<-list(status='fitted',interval=rbind(Dropout=c(lower=.03,upper=.05),Foreign=c(lower=0,upper=.003)))
 d<-.parentage_error_bound_decision(fit,c(Dropout=.005,Foreign=.0001),'raise_only')
 expect_true(d$adjust);expect_equal(d$used_rate,c(Dropout=.03,Foreign=.0001))
 expect_false(d$lower);expect_true(d$parameters$Dropout$raise);expect_false(d$parameters$Foreign$raise)
 fit$status<-'too_few_trios'
 expect_false(.parentage_error_bound_decision(fit,c(Dropout=.005,Foreign=.0001),'raise_only')$adjust)
})

test_that('two-sided policy changes only to the appropriate profile endpoint', {
 fit<-list(status='fitted',interval=c(lower=.01,upper=.03))
 expect_equal(.parentage_error_bound_decision(fit,.001,'two_sided')$used_rate,.01)
 expect_equal(.parentage_error_bound_decision(fit,.1,'two_sided')$used_rate,.03)
 expect_equal(.parentage_error_bound_decision(fit,.02,'two_sided')$used_rate,.02)
 expect_false(.parentage_error_bound_decision(list(status='too_few_trios'),.1,'two_sided')$adjust)
})

test_that('a valid supplied dropout outside the fit branch still has a diagnostic likelihood', {
 t<-build_tables(list(c(.5,.5)),ErrorModel='dropout_foreign',Dropout=.30,Foreign=.015)
 g<-matrix('1/1',12L,1L,dimnames=list(paste0('i',1:12),'L1'))
 p<-matrix(0L,12L,2L);p[3:12,]<-rep(c(1L,2L),each=10L)
 fit<-estimate_parentage_error(g,t,p,MaxMendelE=1L)
 expect_identical(fit$status,'fitted');expect_true(is.finite(fit$loglik_gain))
 expect_true(all(is.finite(fit$interval)));expect_true(all(fit$estimated_rates<=.25))
})
