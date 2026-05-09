# Helper: build a clustered illness-death dataset in interval format.
#
# `design` controls the cluster/group structure of `treatment`:
#   * "case_i"   (default): every cluster carries both groups
#                (groups mixed within cluster). Bakoyannis (2021).
#   * "case_ii":  each cluster carries exactly one group
#                (cluster-randomized). Bakoyannis & Bandyopadhyay (2022).
#   * "case_iii": some clusters mixed, some single-group (NOT supported
#                by patp(); used to test validation errors).
make_interval_data <- function(n_clusters = 20, per_cluster = 8, seed = 1,
                               design = c("case_i", "case_ii", "case_iii")) {
  design <- match.arg(design)
  set.seed(seed)
  rows <- list()
  k <- 1L
  for (c in seq_len(n_clusters)) {
    cluster_grp <- c %% 2     # used for case_ii and to seed case_iii decisions
    for (i in seq_len(per_cluster)) {
      grp <- switch(design,
        case_i   = i %% 2,
        case_ii  = cluster_grp,
        # case_iii: even-numbered clusters mixed (alternate by subject),
        # odd-numbered clusters single-group (cluster_grp).
        case_iii = if (c %% 2 == 0L) i %% 2 else cluster_grp
      )
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
  d <- make_interval_data(n_clusters = 30, per_cluster = 12, seed = 3)
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 200, seed = 3)

  # Group 1 has illness hazard 2x group 0
  expect_lt(fit$test$p.value, 0.10)
})

test_that("patp design = 'shared' works on case (i) data", {
  d <- make_interval_data(design = "case_i")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 50, seed = 1,
              design = "shared")

  expect_false(is.null(fit$test))
  expect_identical(fit$design, "shared")
})

test_that("patp design = 'cluster_random' works on case (ii) data", {
  d <- make_interval_data(design = "case_ii")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 50, seed = 1,
              design = "cluster_random")

  expect_false(is.null(fit$test))
  expect_identical(fit$design, "cluster_random")
})

test_that("patp design = 'indep_random' works on case (ii) data", {
  d <- make_interval_data(design = "case_ii")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 50, seed = 1,
              design = "indep_random")

  expect_false(is.null(fit$test))
  expect_identical(fit$design, "indep_random")
})

test_that("patp design = 'auto' picks 'shared' on case (i) data", {
  d <- make_interval_data(design = "case_i")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
              data = d, tmat = tmat,
              id = "pid", cluster = "site",
              h = 1, j = 2, s = 0, B = 50, seed = 1)

  expect_identical(fit$design, "shared")
})

test_that("patp design = 'auto' picks 'indep_random' on case (ii) data and warns", {
  d <- make_interval_data(design = "case_ii")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_warning(
    fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
                data = d, tmat = tmat,
                id = "pid", cluster = "site",
                h = 1, j = 2, s = 0, B = 50, seed = 1),
    "set design = 'cluster_random' explicitly"
  )
  expect_identical(fit$design, "indep_random")
})

test_that("patp design = 'auto' on case (i) does not warn", {
  d <- make_interval_data(design = "case_i")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_warning(
    fit <- patp(msm(t0, t1, s0, s1) ~ treatment,
                data = d, tmat = tmat,
                id = "pid", cluster = "site",
                h = 1, j = 2, s = 0, B = 50, seed = 1),
    NA
  )
  expect_identical(fit$design, "shared")
})

test_that("patp design = 'indep_random' explicit does not warn", {
  d <- make_interval_data(design = "case_ii")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_warning(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "indep_random"),
    NA
  )
})

test_that("patp design = 'shared' errors on case (ii) data", {
  d <- make_interval_data(design = "case_ii")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "shared"),
    "design = \"shared\" requires every cluster"
  )
})

test_that("patp design = 'cluster_random' errors on case (i) data", {
  d <- make_interval_data(design = "case_i")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "cluster_random"),
    "design = \"cluster_random\" requires each cluster"
  )
})

test_that("patp design = 'indep_random' errors on case (i) data", {
  d <- make_interval_data(design = "case_i")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "indep_random"),
    "design = \"indep_random\" requires each cluster"
  )
})

test_that("patp errors on case (iii) data for every design", {
  d <- make_interval_data(design = "case_iii")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1),
    "Mixed cluster structure"
  )
  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "shared"),
    "design = \"shared\" requires every cluster"
  )
  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "cluster_random"),
    "design = \"cluster_random\" requires each cluster"
  )
  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid", cluster = "site",
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "indep_random"),
    "design = \"indep_random\" requires each cluster"
  )
})

test_that("patp design = 'cluster_random' errors without cluster column", {
  d <- make_interval_data(design = "case_i")
  tmat <- trans_mat(list(c(2, 3), 3, integer(0)))

  expect_error(
    patp(msm(t0, t1, s0, s1) ~ treatment,
         data = d, tmat = tmat,
         id = "pid",                    # no cluster
         h = 1, j = 2, s = 0, B = 50, seed = 1,
         design = "cluster_random"),
    "design = \"cluster_random\" requires a cluster column"
  )
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
