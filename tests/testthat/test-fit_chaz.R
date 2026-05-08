# Helper: build a minimal long-format dataset with known transitions
make_simple_msd <- function() {
  # 5 subjects in illness-death model
  # transitions: 1 = H->I, 2 = H->D, 3 = I->D
  # subj 1: censored healthy at 2.0
  # subj 2: H->D at 1.5
  # subj 3: H->I at 1.0, then I->D at 2.5
  # subj 4: H->I at 0.5, censored ill at 3.0
  # subj 5: censored healthy at 4.0
  data.frame(
    id     = c(1,1, 2,2, 3,3,3, 4,4,4, 5,5),
    Tstart = c(0,0, 0,0, 0,0,1.0, 0,0,0.5, 0,0),
    Tstop  = c(2.0,2.0, 1.5,1.5, 1.0,1.0,2.5, 0.5,0.5,3.0, 4.0,4.0),
    from   = c(1,1, 1,1, 1,1,2, 1,1,2, 1,1),
    to     = c(2,3, 2,3, 2,3,3, 2,3,3, 2,3),
    trans  = c(1,2, 1,2, 1,2,3, 1,2,3, 1,2),
    status = c(0,0, 0,1, 1,0,1, 1,0,0, 0,0)
  )
}

test_that("fit_chaz returns expected columns and types", {
  msd <- make_simple_msd()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  haz <- fit_chaz(msd, tmat)

  expect_s3_class(haz, "data.frame")
  expect_named(haz, c("time", "Haz", "trans"))
  expect_type(haz$trans, "integer")
  expect_true(all(haz$Haz >= 0))
})

test_that("fit_chaz output is sorted by transition then time", {
  msd <- make_simple_msd()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  haz <- fit_chaz(msd, tmat)

  expect_equal(haz$trans, sort(haz$trans))
  for (tr in unique(haz$trans)) {
    times_tr <- haz$time[haz$trans == tr]
    expect_equal(times_tr, sort(times_tr))
  }
})

test_that("fit_chaz cumulative hazards are monotonically non-decreasing", {
  msd <- make_simple_msd()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  haz <- fit_chaz(msd, tmat)

  for (tr in unique(haz$trans)) {
    h_tr <- haz$Haz[haz$trans == tr]
    expect_true(all(diff(h_tr) >= -1e-12))
  }
})

test_that("fit_chaz contains all expected transition IDs", {
  msd <- make_simple_msd()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  haz <- fit_chaz(msd, tmat)

  expect_setequal(unique(haz$trans), c(1L, 2L, 3L))
})

test_that("fit_chaz weights argument is honored", {
  msd <- make_simple_msd()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  unweighted <- fit_chaz(msd, tmat)
  half_w <- fit_chaz(msd, tmat, weights = rep(0.5, nrow(msd)))

  # With uniform weights of 0.5, hazards should equal unweighted hazards
  # (Breslow estimator is invariant to a global rescaling for stratified
  # Cox without covariates).
  expect_equal(half_w$Haz, unweighted$Haz, tolerance = 1e-10)

  # Non-uniform weights should produce different hazards
  varied_w <- runif(nrow(msd), 0.1, 2)
  varied <- fit_chaz(msd, tmat, weights = varied_w)
  expect_false(isTRUE(all.equal(varied$Haz, unweighted$Haz)))
})

test_that("fit_chaz validates inputs", {
  msd <- make_simple_msd()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(fit_chaz(data.frame(x = 1), tmat), "missing columns")
  expect_error(fit_chaz(msd, tmat, weights = 1:3),
               "length nrow")
  expect_error(fit_chaz(msd, tmat, weights = rep(-1, nrow(msd))),
               "non-negative")

  bad <- msd
  bad$trans[1] <- 99
  expect_error(fit_chaz(bad, tmat), "not in tmat")
})
