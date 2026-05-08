make_simple_fit <- function(B = 50, two_sample = FALSE) {
  set.seed(1)
  rows <- list()
  k <- 1L
  for (c in 1:10) {
    for (i in 1:5) {
      t_ill <- rexp(1, 0.5)
      t_dth <- rexp(1, 0.3)
      if (t_ill < t_dth) {
        rows[[length(rows)+1L]] <- data.frame(
          pid = k, site = c, arm = c %% 2,
          t0 = 0, t1 = t_ill, s0 = 1L, s1 = 2L
        )
        rows[[length(rows)+1L]] <- data.frame(
          pid = k, site = c, arm = c %% 2,
          t0 = t_ill, t1 = t_ill + 0.5, s0 = 2L, s1 = 3L
        )
      } else {
        rows[[length(rows)+1L]] <- data.frame(
          pid = k, site = c, arm = c %% 2,
          t0 = 0, t1 = t_dth, s0 = 1L, s1 = 3L
        )
      }
      k <- k + 1L
    }
  }
  d <- do.call(rbind, rows)
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fmla <- if (two_sample) {
    msm(t0, t1, s0, s1) ~ arm
  } else {
    msm(t0, t1, s0, s1) ~ 1
  }

  patp(fmla, data = d, tmat = tmat,
       id = "pid", cluster = "site",
       h = 1, j = 2, s = 0, B = B, seed = 1)
}


test_that("print.patp prints expected sections for one-sample fit", {
  fit <- make_simple_fit()
  out <- capture.output(print(fit))
  expect_true(any(grepl("Population-Averaged", out)))
  expect_true(any(grepl("Subjects:", out)))
  expect_true(any(grepl("Bootstrap", out)))
  expect_true(any(grepl("Curve summary", out)))
})

test_that("print.patp shows test info for two-sample fit", {
  fit <- make_simple_fit(two_sample = TRUE)
  out <- capture.output(print(fit))
  expect_true(any(grepl("Two-sample test", out)))
  expect_true(any(grepl("Statistic:", out)))
  expect_true(any(grepl("p-value:", out)))
})

test_that("print.patp signals B = 0 explicitly", {
  fit <- make_simple_fit(B = 0)
  out <- capture.output(print(fit))
  expect_true(any(grepl("B = 0", out)))
  expect_true(any(grepl("not computed", out)))
})

test_that("summary.patp returns a summary.patp object", {
  fit <- make_simple_fit()
  s <- summary(fit)
  expect_s3_class(s, "summary.patp")
  expect_true("curves" %in% names(s))
  expect_true("estimand" %in% names(s))
})

test_that("print.summary.patp prints the full curves", {
  fit <- make_simple_fit()
  s <- summary(fit)
  out <- capture.output(print(s))
  expect_true(any(grepl("Summary of patp fit", out)))
  expect_true(any(grepl("Full curve", out)))
})

test_that("print.patp works with no clustering (n_clusters = NA)", {
  set.seed(2)
  d <- data.frame(
    pid = 1:20,
    t0  = 0,
    t1  = rexp(20, 0.5),
    s0  = 1L,
    s1  = 2L
  )
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))
  fit <- patp(msm(t0, t1, s0, s1) ~ 1,
              data = d, tmat = tmat, id = "pid",
              h = 1, j = 2, s = 0, B = 30, seed = 2)
  out <- capture.output(print(fit))
  expect_true(any(grepl("no clustering", out)))
})
