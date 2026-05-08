test_that("trans_mat builds the standard illness-death matrix correctly", {
  M <- trans_mat(list(c(2, 3), 3, integer(0)),
                 names = c("Healthy", "Ill", "Dead"))

  expect_equal(dim(M), c(3, 3))
  expect_equal(M["Healthy", "Ill"], 1L)
  expect_equal(M["Healthy", "Dead"], 2L)
  expect_equal(M["Ill", "Dead"], 3L)
  expect_true(is.na(M["Ill", "Healthy"]))
  expect_true(is.na(M["Dead", "Healthy"]))
  expect_true(is.na(M["Healthy", "Healthy"]))
})

test_that("trans_mat supports recovery (cyclic transitions)", {
  M <- trans_mat(list(c(2, 3), c(1, 3), integer(0)),
                 names = c("H", "I", "D"))

  expect_equal(M["I", "H"], 3L)              # recovery transition
  expect_equal(max(M, na.rm = TRUE), 4L)
  expect_false(is.na(M["I", "H"]))
})

test_that("trans_mat assigns IDs in row-major order", {
  M <- trans_mat(list(c(2, 3), c(1, 3), integer(0)))
  expect_equal(M[1, 2], 1L)                  # row 1 transitions first
  expect_equal(M[1, 3], 2L)
  expect_equal(M[2, 1], 3L)                  # row 2 transitions next
})

test_that("trans_mat handles default and custom names", {
  M1 <- trans_mat(list(2, integer(0)))
  expect_equal(dimnames(M1), list(from = c("1", "2"), to = c("1", "2")))

  M2 <- trans_mat(list(2, integer(0)), names = c("A", "B"))
  expect_equal(dimnames(M2), list(from = c("A", "B"), to = c("A", "B")))
})

test_that("trans_mat rejects invalid input", {
  expect_error(trans_mat("not a list"), "must be a list")
  expect_error(trans_mat(list(2)), "at least 2")
  expect_error(trans_mat(list(c(2, 5), integer(0))), "outside 1:2")
  expect_error(trans_mat(list(1, integer(0))), "self-transition")
  expect_error(trans_mat(list(c(2, 2), integer(0))), "duplicated")
  expect_error(trans_mat(list(2, integer(0)), names = c("A")), "length 2")
  expect_error(trans_mat(list(2, integer(0)), names = c("A", "A")), "unique")
  expect_error(trans_mat(list(integer(0), integer(0))), "no transitions")
})

test_that("trans_mat output works as input to existing patp pipeline", {
  M <- trans_mat(list(c(2, 3), 3, integer(0)),
                 names = c("Healthy", "Ill", "Dead"))
  expect_true(is.matrix(M))
  expect_true(is.integer(M[!is.na(M)]))
})
