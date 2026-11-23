test_that('parentage requires exact years or explicit inclusive bounds', {
  lh <- PaPeR:::.normalize_life_history(data.frame(
    ID=c('a','b'),Sex=c(1,2),BirthYear=c(1L,NA_integer_)))
  expect_error(PaPeR:::.require_birth_bounds(lh),'missing for: b')
  lh$BY.min[2]<-0L;lh$BY.max[2]<-4L
  expect_silent(PaPeR:::.require_birth_bounds(lh))
})

test_that('typed proposals propagate through offspring and relax removed links', {
  parents<-matrix(c(0L,0L,0L,2L,0L,0L),3L,2L)
  make<-function(p)PaPeR:::typed_age_state_cpp(p,c(2L,0L,2L),c(2L,4L,4L),
    rep(NA_integer_,3),c(1,1),c(5,5))
  pointer<-make(parents)
  q<-PaPeR:::typed_age_query_cpp(pointer,0L,integer())
  expect_true(q$feasible);expect_identical(q$birth_hi,c(2L,1L,4L))
  proposed<-PaPeR:::typed_age_query_cpp(pointer,2L,c(3L,0L))
  expect_false(proposed$feasible)
  # The query is transactional and cannot leave the state partially changed.
  expect_identical(PaPeR:::typed_age_query_cpp(pointer,0L,integer()),q)
  cut<-PaPeR:::typed_age_query_cpp(pointer,1L,c(0L,0L))
  expect_true(cut$feasible);expect_identical(cut$birth_hi,c(2L,4L,4L))
  parents[1,2]<-0L
  expect_true(PaPeR:::typed_age_query_cpp(make(parents),2L,c(3L,0L))$feasible)
})

test_that('public fixed scoring rejects the R09 age chain before numerical scoring', {
  calls <- rbind(A = c('1/2', '1/2'), D = c('1/1', '1/1'), C = c('1/2', '1/2'))
  colnames(calls) <- c('L1', 'L2')
  life <- data.frame(ID = rownames(calls), Sex = c(NA_integer_, 1L, NA_integer_),
    BirthYear = c(NA_integer_, NA_integer_, 2L),
    BY.min = c(0L, 2L, 2L), BY.max = c(4L, 4L, 2L))
  pedigree <- data.frame(id = rownames(calls), dam = c('D', NA_character_, NA_character_),
    sire = c(NA_character_, NA_character_, 'A'))
  control <- list(error_model = 'sequoia', Err = .01,
    frequencies = rep(list(c(.5, .5)), 2), error_calibration = 'none',
    MinAgeParent = 1, MaxAgeParent = 5)
  expect_error(calc_llr(pedigree, calls, life, control = control),
    'impossible_age_component')
})
