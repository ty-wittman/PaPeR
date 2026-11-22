test_that('converged question records retain scored pairs across blank skips', {
 p<-matrix(c(0L,0L,1L,0L,0L,2L),3L,2L)
 i<-list(likelihood_evidence='ancestral',identity_margin_ln=1,Tassign_ln=1)
 r<-data.frame(focal=3L,selected_dam=1L,selected_sire=2L,
  graph_revision_after=2L,complete_pair_skipped=0L,identity_margin_ln=5,full_relationship_ln=6,
  dam_support_ln=10,sire_support_ln=11,pair_support_ln=12)
 blank<-r;blank$complete_pair_skipped<-1L;blank$identity_margin_ln<- -Inf;blank$full_relationship_ln<-NA_real_;blank$graph_revision_after<-2L
 a<-.parentage_confidence_from_records(c('d','s','c'),p,i,list(r,blank),TRUE)
 expect_identical(a$confident,c(FALSE,FALSE,TRUE))
 expect_equal(a$native_audit_calls,0L);expect_equal(a$diagnostics$evaluated_revision,2L)
 expect_null(.parentage_confidence_from_records(c('d','s','c'),p,i,list(r),FALSE))
 i$likelihood_evidence<-'full'
 expect_null(.parentage_confidence_from_records(c('d','s','c'),p,i,list(r),TRUE))
 i$likelihood_evidence<-'ancestral';wrong<-r;wrong$selected_sire<-1L;wrong$graph_revision_after<-2L
 expect_null(.parentage_confidence_from_records(c('d','s','c'),p,i,list(r,wrong),TRUE))
 expect_null(.parentage_confidence_from_records(c('d','s','c'),p,i,list(blank),TRUE))
})

test_that('confidence is certified on a fixed final pedigree without mutations', {
 source(test_path('helper-parentage-fixture.R'),local=TRUE)
 calls<-rbind(A=rep('1/2',60L),D=rep('1/1',60L),S=rep('2/2',60L));colnames(calls)<-paste0('L',1:60)
 f<-.scientific_fixture(calls,dam_candidates='D',sire_candidates='S',Complex='full',birth_year=c(A=4L,D=2L,S=2L))
 i<-f$input;i$pedigree_likelihood<-TRUE;i$likelihood_evidence<-'ancestral';i$pedigree_hypothesis_suite<-'S7';i<- .engine_defaults(i)
 i$local_conditioning<-i$local_core_only<-TRUE
 p<-i$parents;p[1,]<-c(2L,3L)
 z<-.parentage_confident_final_trios(f$compact_inputs,i,p)
 expect_true(z$frozen_state_nonmutation)
 expect_identical(z$confident,c(TRUE,FALSE,FALSE))
 expect_identical(i$parents,f$input$parents)
 expect_gt(z$diagnostics$identity_ln,z$identity_threshold_ln)
 expect_gt(z$diagnostics$relationship_ln,z$relationship_threshold_ln)
})

test_that('the reverse policy uses the upper bound and never raises', {
 fit<-function(point,lower,upper)list(status='fitted',estimate=point,interval=c(lower=lower,upper=upper))
 d<-.parentage_error_bound_decision(fit(.001,.0008,.0012),.05,'lower_only')
 expect_true(d$adjust);expect_true(d$lower);expect_false(d$raise)
 expect_equal(d$used_rate,.0012)
 expect_identical(d$status,'lowered_to_upper_bound')
 expect_false(.parentage_error_bound_decision(fit(.02,.01,.06),.05,'lower_only')$adjust)
 expect_false(.parentage_error_bound_decision(fit(.02,.01,.05),.05,'lower_only')$adjust)
 expect_false(.parentage_error_bound_decision(fit(.02,.01,.03),.001,'lower_only')$adjust)
 expect_false(.parentage_error_bound_decision(list(status='too_few_trios'),.05,'lower_only')$adjust)
})

test_that('bounded frozen batches agree with independent final-trio audits', {
 source(test_path('helper-parentage-fixture.R'),local=TRUE)
 calls<-rbind(A=rep('1/2',60L),D=rep('1/1',60L),S=rep('2/2',60L),O=rep('1/2',60L));colnames(calls)<-paste0('L',1:60)
 f<-.scientific_fixture(calls,dam=c(A='D',O='D'),sire=c(A='S',O='S'),
  dam_candidates='D',sire_candidates='S',Complex='full',birth_year=c(A=4L,D=2L,S=2L,O=4L))
 i<-f$input;i$pedigree_likelihood<-TRUE;i$likelihood_evidence<-'ancestral';i$pedigree_hypothesis_suite<-'S7';i<-.engine_defaults(i)
 i$local_conditioning<-i$local_core_only<-TRUE
 i$diagnostic_start_position<-1L;i$diagnostic_end_position<-1L
 # The helper deliberately declares candidates only for A; the second
 # audit also needs its incumbent pair in its candidate arrays.
 for(side in c('dam','sire')) {
  i[[paste0(side,'_ids')]]<-c(i[[paste0(side,'_ids')]],if(side=='dam')2L else 3L)
  i[[paste0(side,'_offsets')]][5L]<-2L
  i[[paste0(side,'_quick')]]<-c(0,0)
  i[[paste0(side,'_guard')]]<-c(-Inf,-Inf)
  i[[paste0(side,'_permanent')]]<-c(TRUE,TRUE)
 }
 a<-.parentage_confident_final_trios(f$compact_inputs,i,i$parents,batch_size=1L)
 b<-.parentage_confident_final_trios(f$compact_inputs,i,i$parents,batch_size=16L)
 expect_identical(a$confident,b$confident)
 expect_equal(a$diagnostics,b$diagnostics,tolerance=1e-6)
 expect_equal(a$native_audit_calls,2L);expect_equal(b$native_audit_calls,1L)
 expect_true(a$frozen_state_nonmutation);expect_true(b$frozen_state_nonmutation)
})

test_that('the scalar calibration specification is SNP-only; MH keeps the two-rate estimator', {
  snp <- build_tables(list(c(.3,.7)),Err=.01)
  spec <- .parentage_error_scale_spec(snp)
  expect_identical(spec$parameter,'Err')
  expect_equal(spec$initial,.01)
  expect_equal(spec$rebuild(.02)[[1]]$Err,.02)
  expect_equal(spec$rebuild(.02)[[1]]$p,c(.3,.7))
  mh <- build_tables(list(c(.3,.7)),ErrorModel='dropout_foreign',Dropout=.02,Foreign=.001)
  expect_error(.parentage_error_scale_spec(mh),'SNP-only')
})
