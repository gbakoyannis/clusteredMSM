# Helper to construct a small msdata-style example
make_msd <- function() {
  # 4 subjects, illness-death model, transitions 1=H->I, 2=H->D, 3=I->D
  data.frame(
    id     = c(1,1, 2,2, 3,3,3, 4,4),
    Tstart = c(0,0, 0,0, 0,0,1.6, 0,0),
    Tstop  = c(1.9,1.9, 2.2,2.2, 1.6,1.6,2.7, 0.4,0.4),
    from   = c(1,1, 1,1, 1,1,2, 1,1),
    to     = c(2,3, 2,3, 2,3,3, 2,3),
    trans  = c(1,2, 1,2, 1,2,3, 1,2),
    status = c(0,0, 0,1, 1,0,0, 1,0)
  )
}

test_that("state_at returns starting state at s = 0", {
  msd <- make_msd()
  out <- state_at(msd, s = 0.0)
  expect_equal(nrow(out), 4)
  expect_true(all(out$state == 1))
})

test_that("state_at identifies subjects who have transitioned to illness", {
  msd <- make_msd()
  out <- state_at(msd, s = 2.0)
  expect_equal(out$state[out$id == 3], 2)
})

test_that("state_at excludes absorbed subjects", {
  msd <- make_msd()
  out <- state_at(msd, s = 1.0)
  expect_false(4 %in% out$id)
})

test_that("state_at excludes subjects whose follow-up ended", {
  msd <- make_msd()
  out <- state_at(msd, s = 3.0)
  expect_false(2 %in% out$id)
})

test_that("state_at deduplicates competing-transition rows", {
  msd <- make_msd()
  out <- state_at(msd, s = 0.5)
  expect_equal(anyDuplicated(out$id), 0)
})

test_that("state_at validates input", {
  expect_error(state_at(data.frame(x = 1), s = 0), "missing columns")
  expect_error(state_at(make_msd(), s = c(1, 2)), "numeric scalar")
  expect_error(state_at(make_msd(), s = "a"), "numeric scalar")
})
