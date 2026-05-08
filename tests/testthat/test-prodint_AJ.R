# Helper: build a small Haz data frame manually
# tmat: 3 states, transitions 1=1->2, 2=1->3, 3=2->3
make_haz <- function() {
  data.frame(
    time  = c(0.5, 1.0, 1.5,    # trans 1 jumps
              0.8, 2.0,         # trans 2 jumps
              1.2, 1.8),        # trans 3 jumps
    Haz   = c(0.1, 0.25, 0.5,
              0.05, 0.2,
              0.15, 0.4),
    trans = c(1, 1, 1,
              2, 2,
              3, 3)
  )
}

test_that("prodint_AJ output has correct shape and columns", {
  haz <- make_haz()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  P <- prodint_AJ(haz, tmat, predt = 0, h = 1)

  expect_s3_class(P, "data.frame")
  expect_named(P, c("time", "pstate1", "pstate2", "pstate3"))
  expect_equal(nrow(P), length(unique(haz$time)) + 1L)
})

test_that("prodint_AJ rows sum to 1 at every time point", {
  haz <- make_haz()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  P <- prodint_AJ(haz, tmat, predt = 0, h = 1)

  row_sums <- rowSums(P[, paste0("pstate", 1:3)])
  expect_equal(row_sums, rep(1, nrow(P)), tolerance = 1e-10)
})

test_that("prodint_AJ starts at indicator vector e_h", {
  haz <- make_haz()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  for (h in 1:3) {
    P <- prodint_AJ(haz, tmat, predt = 0, h = h)
    expected <- numeric(3); expected[h] <- 1
    expect_equal(unname(unlist(P[1, paste0("pstate", 1:3)])),
                 expected, tolerance = 1e-12)
  }
})

test_that("prodint_AJ probabilities are in [0, 1]", {
  haz <- make_haz()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  P <- prodint_AJ(haz, tmat, predt = 0, h = 1)

  vals <- as.matrix(P[, paste0("pstate", 1:3)])
  expect_true(all(vals >= -1e-12 & vals <= 1 + 1e-12))
})

test_that("prodint_AJ absorbing state is monotonically non-decreasing", {
  haz <- make_haz()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  P <- prodint_AJ(haz, tmat, predt = 0, h = 1)

  # State 3 is absorbing: P(state 3) should never decrease
  expect_true(all(diff(P$pstate3) >= -1e-12))
})

test_that("prodint_AJ predt argument truncates correctly", {
  haz <- make_haz()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  P0 <- prodint_AJ(haz, tmat, predt = 0,   h = 1)
  P1 <- prodint_AJ(haz, tmat, predt = 1.0, h = 1)

  expect_true(all(P1$time >= 1.0))
  expect_true(nrow(P1) < nrow(P0))
})

test_that("prodint_AJ supports recovery (cyclic transitions)", {
  # Recovery model: 1<->2, 1->3, 2->3
  haz <- data.frame(
    time  = c(0.5, 1.5,    # trans 1: 1->2
              1.0, 2.0,    # trans 2: 1->3
              0.8, 1.8,    # trans 3: 2->1 (recovery)
              1.2, 2.5),   # trans 4: 2->3
    Haz   = c(0.1, 0.3,
              0.05, 0.15,
              0.08, 0.2,
              0.1, 0.25),
    trans = c(1, 1, 2, 2, 3, 3, 4, 4)
  )
  tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)))

  P <- prodint_AJ(haz, tmat, predt = 0, h = 1)

  expect_equal(rowSums(P[, paste0("pstate", 1:3)]),
               rep(1, nrow(P)), tolerance = 1e-10)

  # State 3 (Dead) is still absorbing -> non-decreasing
  expect_true(all(diff(P$pstate3) >= -1e-12))

  # Probability of being in state 1 (Healthy) need not be monotone --
  # this is the whole point of recovery. Just check it dips below 1.
  expect_true(min(P$pstate1) < 1)
})

test_that("prodint_AJ rejects invalid input", {
  haz <- make_haz()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(prodint_AJ(list(), tmat), "data.frame")
  expect_error(prodint_AJ(haz, "not a matrix"), "matrix")
  expect_error(prodint_AJ(haz, tmat, h = 99), "h %in%")
  expect_error(prodint_AJ(haz, tmat, predt = c(0, 1)), "length\\(predt\\)")

  bad <- haz; bad$trans[1] <- 99
  expect_error(prodint_AJ(bad, tmat), "not present in tmat")
})
