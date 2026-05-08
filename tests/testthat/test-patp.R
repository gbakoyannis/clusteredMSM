# Helper: build a clustered illness-death dataset in interval format
make_interval_data <- function(n_clusters = 20, per_cluster = 8, seed = 1) {
  set.seed(seed)
  rows <- list()
  k <- 1L
  for (c in seq_len(n_clusters)) {
    grp <- c %% 2
    for (i in seq_len(per_cluster)) {
      ill_rate <- if (grp == 1) 0.8 else 0.4
      t_ill   <- rexp(1, rate = ill_rate)
      t_dth_h <- rexp(1, rate = 0.3)
      t_dth_i <- rexp(1, rate = 0.6)
      t_cens  <- runif(1, 1, 5)

      ill_first <- t_ill < t_dth_h && t_ill < t_cens
      dth_first <- t_dth_h < t_ill && t_dth_h < t_cens

      if (ill_first) {
        # H -> I, then I -> D or censored ill
        rows[[length(rows) + 1L]] <- data.frame(
          pid = k, site = c, treatment = grp,
          t0 = 0, t1 = t_ill, s0 = 1L, s1 = 2L
        )
        dth_t <- t_ill + t_dth_i
        if (dth_t < t_cens) {
          rows[[length(rows) + 1L]] <- data.frame(
            pid = k, site = c, treatment = grp,
            t0 = t_ill, t1 = dth_t, s0 = 2L, s1 = 3L
          )
        } else {
          rows[[length(rows) + 1L]] <- data.frame(
            pid = k, site = c, treatment = grp,
            t0 = t_ill, t1 = t_cens, s0 = 2L, s1 = 2L
          )
        }
      } else if (dth_first) {
        rows[[length(rows) + 1L]] <- data.frame(
          pid = k, site = c, treatment = grp,
          t0 = 0, t1 = t_dth_h, s0 = 1L, s1 = 3L
        )
      } else {
        # Censored healthy
        rows[[length(rows) + 1L]] <- data.frame(
          pid = k, site = c, treatment = grp,
          t0 = 0, t1 = t_cens, s0 = 1L, s1 = 1L
        )
      }
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}


test_that("patp one-sample formula returns expected structure", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ 1,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 50, seed = 1)

  expect_s3_class(fit, "patp")
  expect_named(fit, c("curves", "test", "n_clusters", "groups",
                      "call", "formula", "h", "j", "s", "B",
                      "n_subjects"),
               ignore.order = TRUE)
  expect_null(fit$test)
  expect_true(all(c("time", "P", "se", "ll", "ul") %in% names(fit$curves)))
})

test_that("patp two-sample formula returns test slot", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 50, seed = 1)

  expect_s3_class(fit, "patp")
  expect_false(is.null(fit$test))
  expect_named(fit$test, c("statistic", "p.value", "type"))
  expect_true(fit$test$p.value >= 0 && fit$test$p.value <= 1)
  expect_true("group" %in% names(fit$curves))
  expect_setequal(unique(fit$curves$group), c(0, 1))
})

test_that("patp B = 0 returns point estimate only", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ 1,
              data = d, tmat = tmat, id = "pid",
              h = 1, j = 2, s = 0, B = 0)

  expect_named(fit$curves, c("time", "P"))
  expect_equal(fit$B, 0)
})

test_that("patp errors on two-sample formula with B = 0", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat, id = "pid",
         h = 1, j = 2, s = 0, B = 0),
    "two-sample formulas require B > 0"
  )
})

test_that("patp errors on weighted = TRUE without cluster", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(
    patp(msm(t0, t1, s0, s1) ~ 1,
         data = d, tmat = tmat, id = "pid",
         h = 1, j = 2, s = 0, weighted = TRUE, B = 50),
    "weighted = TRUE requires"
  )
})

test_that("patp without cluster uses individual-level bootstrap", {
  d <- make_interval_data(n_clusters = 5, per_cluster = 5)
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ 1,
              data = d, tmat = tmat, id = "pid",
              h = 1, j = 2, s = 0, B = 50, seed = 1)

  expect_true(is.na(fit$n_clusters))
  expect_true(all(c("se", "ll", "ul") %in% names(fit$curves)))
})

test_that("patp probabilities are in [0, 1]", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ 1,
              data = d, tmat = tmat, id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 0)

  expect_true(all(fit$curves$P >= -1e-10 & fit$curves$P <= 1 + 1e-10))
})

test_that("patp seed produces reproducible results", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit1 <- patp(msm(t0, t1, s0, s1) ~ 1,
               data = d, tmat = tmat,
               id = "pid", cluster = "site",
               h = 1, j = 2, s = 0, B = 30, seed = 42)
  fit2 <- patp(msm(t0, t1, s0, s1) ~ 1,
               data = d, tmat = tmat,
               id = "pid", cluster = "site",
               h = 1, j = 2, s = 0, B = 30, seed = 42)

  expect_equal(fit1$curves$se, fit2$curves$se)
})

test_that("patp two-sample detects substantial group differences", {
  d <- make_interval_data(n_clusters = 30, per_cluster = 12, seed = 5)
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 200, seed = 5)

  # Group 1 has illness hazard 2x group 0
  expect_lt(fit$test$p.value, 0.10)
})

test_that("patp runs with cband = TRUE and produces band columns", {
  d <- make_interval_data()
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- suppressWarnings(
    patp(msm(t0, t1, s0, s1) ~ 1,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 100, cband = TRUE, seed = 1)
  )

  expect_true(all(c("ll.band", "ul.band") %in% names(fit$curves)))
  # Some band values are NA outside trim; some are not
  expect_true(any(!is.na(fit$curves$ll.band)))
})
