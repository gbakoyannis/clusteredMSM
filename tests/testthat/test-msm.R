test_that("msm constructs a valid object with correct columns", {
  obj <- msm(c(0, 1.5), c(1.5, 3.0), c(1, 2), c(2, 3))

  expect_s3_class(obj, "msm")
  expect_true(is.matrix(obj))
  expect_equal(colnames(obj), c("Tstart", "Tstop", "Sstart", "Sstop"))
  expect_equal(nrow(obj), 2)
})

test_that("msm accepts integer-valued states and rejects fractional ones", {
  # Whole-number doubles are accepted (matrix storage is double; values
  # are validated to be integer-valued).
  expect_silent(msm(c(0, 1), c(1, 2), c(1.0, 2.0), c(2.0, 3.0)))
  # Fractional state values are rejected.
  expect_error(msm(c(0, 1), c(1, 2), c(1.5, 2.0), c(2.0, 3.0)),
               "integer-valued")
})

test_that("msm errors on length mismatch", {
  expect_error(msm(c(0, 1), c(1, 2, 3), c(1, 1), c(2, 2)),
               "same length")
})

test_that("msm errors on missing arguments", {
  expect_error(msm(), "four arguments")
  expect_error(msm(0), "four arguments")
})

test_that("msm errors when Tstart >= Tstop", {
  expect_error(msm(c(0, 2), c(1, 1), c(1, 1), c(2, 2)),
               "strictly less than")
  expect_error(msm(c(0), c(0), c(1), c(2)),
               "strictly less than")
})

test_that("msm errors on non-numeric times", {
  expect_error(msm(c("a", "b"), c(1, 2), c(1, 1), c(2, 2)),
               "numeric")
})

test_that("msm.print runs without error", {
  obj <- msm(c(0, 1, 2), c(1, 2, 3), c(1, 1, 2), c(1, 2, 3))
  expect_output(print(obj), "msm object")
})

test_that("msm allows censoring rows (Sstart == Sstop)", {
  # Censoring is permitted at the constructor level; it's the
  # validator's job to enforce that it's only on the last row.
  obj <- msm(0, 1, 1, 1)
  expect_equal(unname(obj[1, "Sstart"]), 1L)
  expect_equal(unname(obj[1, "Sstop"]), 1L)
})
