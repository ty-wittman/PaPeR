
test_that('confidence requires a current evaluation and never resurrects an unavailable row', {
  p<-matrix(c(0L,0L,1L,0L,0L,2L),3,2)
  i<-list(likelihood_evidence='ancestral',identity_margin_ln=1,Tassign_ln=1)
  r<-data.frame(focal=3L,selected_dam=1L,selected_sire=2L,graph_revision_after=1L,
    complete_pair_skipped=0L,identity_margin_ln=5,full_relationship_ln=6,
    dam_support_ln=10,sire_support_ln=11,pair_support_ln=12)
  run<-function(rows,revision=2L,current=TRUE).parentage_confidence_from_records(
    c('d','s','c'),p,i,rows,TRUE,final_revision=revision,current_records=current)
  expect_null(run(list(r)))
  latest<-r;latest$graph_revision_after<-2L;latest$identity_margin_ln<-latest$full_relationship_ln<-NA_real_
  expect_null(run(list(r,latest)))
  expect_true(run(list(r,latest),current=FALSE)$confident[3])
  r$graph_revision_after<-2L
  skip<-latest;skip$complete_pair_skipped<-1L
  expect_true(run(list(r,skip))$confident[3])
  expect_null(run(list(r,latest)))
})
