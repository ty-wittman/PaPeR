.confidence_record_fixture <- function() {
  ids<-c('A','B','D','S','G','H')
  p<-matrix(0L,6L,2L);p[1:2,]<-rep(c(3L,4L),each=2L)
  d<-data.frame(focal=1:2,selected_dam=3L,selected_sire=4L,
    graph_revision_after=c(3L,5L),complete_pair_skipped=0L,
    identity_margin_ln=10,full_relationship_ln=10,
    dam_support_ln=10,sire_support_ln=10,pair_support_ln=10)
  state<-list(schema_version='complete-trio-dependency-v1',enabled=TRUE,
    converged=TRUE,final_graph_revision=9L,final_parents=p,
    last_evaluated_revision=c(NA_integer_,3L,5L,rep(NA_integer_,4L)),
    recheck_pending=integer(7L))
  list(ids=ids,parents=p,input=list(likelihood_evidence='ancestral',identity_margin_ln=.5,Tassign_ln=1),
    blocks=list(d),converged=TRUE,final_revision=9L,dependency_state=state,allow_partial=TRUE)
}

test_that('unrelated revisions do not expire dependency-current trio records', {
  f<-.confidence_record_fixture();x<-do.call(.parentage_confidence_from_records,f)
  expect_identical(which(x$confident),1:2)
  expect_identical(x$needs_audit,integer())
  expect_equal(x$diagnostics$evaluated_revision,c(3L,5L))
  f$dependency_state<-NULL;f$allow_partial<-FALSE
  expect_null(do.call(.parentage_confidence_from_records,f))
})

test_that('stale, missing and differently parented records need individual audits', {
  f<-.confidence_record_fixture();f$dependency_state$recheck_pending[3L]<-1L
  x<-do.call(.parentage_confidence_from_records,f)
  expect_identical(x$needs_audit,2L);expect_identical(which(x$confident),1L)
  f<-.confidence_record_fixture();f$blocks[[1L]]<-f$blocks[[1L]][1L,,drop=FALSE]
  expect_identical(do.call(.parentage_confidence_from_records,f)$needs_audit,2L)
  f<-.confidence_record_fixture();f$blocks[[1L]]$selected_sire[2L]<-5L
  expect_identical(do.call(.parentage_confidence_from_records,f)$needs_audit,2L)
  f<-.confidence_record_fixture();f$dependency_state$last_evaluated_revision[3L]<-6L
  expect_identical(do.call(.parentage_confidence_from_records,f)$needs_audit,2L)
  f<-.confidence_record_fixture();f$blocks[[1L]]$full_relationship_ln[2L]<-NA_real_
  expect_identical(do.call(.parentage_confidence_from_records,f)$needs_audit,2L)
})

test_that('a dependency certificate must describe the same final graph', {
  f<-.confidence_record_fixture();f$dependency_state$final_graph_revision<-8L
  expect_identical(do.call(.parentage_confidence_from_records,f)$needs_audit,1:2)
  f<-.confidence_record_fixture();f$dependency_state$final_parents[1L,2L]<-5L
  expect_identical(do.call(.parentage_confidence_from_records,f)$needs_audit,1:2)
  f<-.confidence_record_fixture();f$dependency_state$enabled<-FALSE
  expect_identical(do.call(.parentage_confidence_from_records,f)$needs_audit,1:2)
})

test_that('only the pending trio enters a frozen audit', {
  f<-.confidence_record_fixture();f$dependency_state$recheck_pending[3L]<-1L
  compact<-list(n_individuals=6L,id_index=setNames(seq_along(f$ids),f$ids));seen<-NULL
  local_mocked_bindings(.parentage_confident_final_trios=function(compact,input,parents,
    extra_ln=log(10),batch_size=16L,focals=NULL){
    seen<<-focals;c<-rep(FALSE,6L);c[focals]<-TRUE
    list(confident=c,diagnostics=data.frame(id=f$ids[focals],confident=TRUE,
      identity_ln=10,relationship_ln=10,dam_support_ln=10,sire_support_ln=10,pair_support_ln=10),
      native_audit_calls=1L,native_audit_focals=length(focals),numerical_work=c(questions=1))
  },.package='PaPeR')
  x<-.parentage_final_confidence(compact,f$input,f$parents,f$blocks,TRUE,9L,f$dependency_state)
  expect_identical(seen,2L);expect_identical(which(x$confident),1:2)
  expect_identical(x$native_audit_focals,1L)
  expect_identical(x$reused_focals,1L)
})

test_that('an empty audit performs no native preparation', {
  f<-.confidence_record_fixture();compact<-list(n_individuals=6L,id_index=setNames(seq_along(f$ids),f$ids))
  local_mocked_bindings(build_snapshot_cpp=function(...)stop('unexpected native compilation'),.package='PaPeR')
  x<-.parentage_confident_final_trios(compact,f$input,f$parents,focals=integer())
  expect_identical(x$native_audit_calls,0L)
  expect_identical(x$native_audit_focals,0L)
  expect_false(any(x$confident))
})

