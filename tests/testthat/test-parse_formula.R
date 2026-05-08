make_test_data <- function() {
  data.frame(
    pid    = c(1, 1, 2, 3),
    site   = c(1, 1, 1, 2),
    arm    = c("A", "A", "B", "B"),
    t0     = c(0,   0.5, 0,   0),
    t1     = c(0.5, 1.5, 1.0, 2.0),
    s0     = c(1,   2,   1,   1),
    s1     = c(2,   3,   3,   1),
    stringsAsFactors = FALSE
  )
}

test_that("parse_msm_formula handles one-sample formulas", {
  d <- make_test_data()
  out <- parse_msm_formula(msm(t0, t1, s0, s1) ~ 1, d, id = "pid")

  expect_named(out, c("data", "group", "has_cluster", "group_name"))
  expect_null(out$group)
  expect_false(out$has_cluster)
  expect_named(out$data, c("id", "Tstart", "Tstop", "Sstart", "Sstop"))
  expect_equal(out$data$Tstart, c(0, 0.5, 0, 0))
})

test_that("parse_msm_formula handles two-sample formulas", {
  d <- make_test_data()
  out <- parse_msm_formula(msm(t0, t1, s0, s1) ~ arm, d, id = "pid")

  expect_equal(out$group, c("A", "A", "B", "B"))
  expect_equal(out$group_name, "arm")
  expect_true("group" %in% names(out$data))
})

test_that("parse_msm_formula picks up cluster column when supplied", {
  d <- make_test_data()
  out <- parse_msm_formula(msm(t0, t1, s0, s1) ~ 1, d,
                           id = "pid", cluster = "site")

  expect_true(out$has_cluster)
  expect_true("cluster" %in% names(out$data))
  expect_equal(out$data$cluster, c(1, 1, 1, 2))
})

test_that("parse_msm_formula treats NA cluster as absent", {
  d <- make_test_data()
  out <- parse_msm_formula(msm(t0, t1, s0, s1) ~ 1, d,
                           id = "pid", cluster = NA)
  expect_false(out$has_cluster)
  expect_false("cluster" %in% names(out$data))
})

test_that("parse_msm_formula errors when LHS is not msm()", {
  d <- make_test_data()
  expect_error(parse_msm_formula(t0 ~ 1, d, id = "pid"),
               "left-hand side must be a call to msm")
  expect_error(parse_msm_formula(Surv(t0, t1) ~ 1, d, id = "pid"),
               "left-hand side must be a call to msm")
})

test_that("parse_msm_formula errors on wrong number of msm arguments", {
  d <- make_test_data()
  expect_error(parse_msm_formula(msm(t0, t1, s0) ~ 1, d, id = "pid"),
               "exactly 4 arguments")
})

test_that("parse_msm_formula errors when id missing or invalid", {
  d <- make_test_data()
  expect_error(parse_msm_formula(msm(t0, t1, s0, s1) ~ 1, d, id = NULL),
               "single character string")
  expect_error(parse_msm_formula(msm(t0, t1, s0, s1) ~ 1, d,
                                 id = "missing_col"),
               "not found")
})

test_that("parse_msm_formula errors on multiple RHS terms", {
  d <- make_test_data()
  d$z <- 1:4
  expect_error(parse_msm_formula(msm(t0, t1, s0, s1) ~ arm + z, d,
                                 id = "pid"),
               "single grouping variable")
})

test_that("parse_msm_formula errors when group has != 2 levels", {
  d <- make_test_data()
  d$arm <- "A"
  expect_error(parse_msm_formula(msm(t0, t1, s0, s1) ~ arm, d, id = "pid"),
               "exactly 2 distinct")
})

test_that("parse_msm_formula errors when grouping variable not in data", {
  d <- make_test_data()
  expect_error(parse_msm_formula(msm(t0, t1, s0, s1) ~ ghost, d, id = "pid"),
               "not found")
})
