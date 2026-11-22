test_that('confidence audits retain decisions and scores without raw contractions', {
 source(test_path('helper-parentage-fixture.R'),local=TRUE)
 calls<-rbind(A=rep('1/2',60L),D=rep('1/1',60L),S=rep('2/2',60L));colnames(calls)<-paste0('L',1:60)
 f<-.scientific_fixture(calls,dam_candidates='D',sire_candidates='S',Complex='full',birth_year=c(A=4L,D=2L,S=2L))
 i<-f$input;profile<-.engine_profile();i[names(profile)]<-profile;i<-.engine_defaults(i)
 p<-i$parents;p[1,]<-c(2L,3L)
 audit<-function(skip){
  e<-new.env(parent=asNamespace('PaPeR'));frames<-list()
  e$run_parentage_cpp<-function(...){z<-PaPeR:::run_parentage_cpp(...)
   frames[[length(frames)+1L]]<<-as.data.frame(z$configuration_diagnostics)
   z}
  run<-PaPeR:::.parentage_confident_final_trios;environment(run)<-e
  j<-i;j$likelihood_confidence_raw_elision<-skip
  r<-run(f$compact_inputs,j,p)
  list(result=r,rows=do.call(rbind,frames))
 }
 old<-audit(FALSE);slim<-audit(TRUE)
 # Raw-score availability replaces the retired per-visit planning receipt.
 expect_true(any(old$rows$raw_comparison_available))
 expect_true(any(is.finite(old$rows$objective_ln)))
 expect_true(all(is.na(slim$rows$objective_ln)))
 expect_true(all(!slim$rows$raw_comparison_available))
 expect_identical(slim$result$confident,old$result$confident)
 expect_identical(slim$result$confident,c(TRUE,FALSE,FALSE))
 expect_equal(slim$result$diagnostics,old$result$diagnostics,tolerance=1e-9)
 expect_true(slim$result$frozen_state_nonmutation)
 expect_identical(i$parents,f$input$parents)
})
