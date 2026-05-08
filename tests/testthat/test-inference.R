test_that("ci_cloglog returns symmetric-ish bounds around point", {
  point <- c(0.2, 0.5, 0.8)
  se    <- c(0.05, 0.05, 0.05)
  ci <- ci_cloglog(point, se)

  expect_named(ci, c("ll", "ul"))
  expect_true(all(ci$ll < point))
  expect_true(all(ci$ul > point))
  expect_true(all(ci$ll > 0 & ci$ul < 1))
})

test_that("ci_cloglog handles boundary point estimates with NA", {
  point <- c(0, 0.5, 1)
  se    <- c(0.1, 0.1, 0.1)
  ci <- ci_cloglog(point, se)

  expect_true(is.na(ci$ll[1]))
  expect_true(is.na(ci$ul[1]))
  expect_false(is.na(ci$ll[2]))
  expect_true(is.na(ci$ll[3]))
})

test_that("ci_cloglog narrows with smaller SE", {
  point <- rep(0.5, 3)
  ci_wide   <- ci_cloglog(point, rep(0.1,  3))
  ci_narrow <- ci_cloglog(point, rep(0.01, 3))

  expect_true(all(ci_narrow$ul - ci_narrow$ll <
                  ci_wide$ul   - ci_wide$ll))
})

test_that("ci_cloglog level argument works", {
  point <- 0.5
  se    <- 0.05
  ci_95 <- ci_cloglog(point, se, level = 0.95)
  ci_99 <- ci_cloglog(point, se, level = 0.99)

  expect_true(ci_99$ul - ci_99$ll > ci_95$ul - ci_95$ll)
})

test_that("ci_cloglog validates input lengths", {
  expect_error(ci_cloglog(c(0.3, 0.5), c(0.1)), "same length")
})


test_that("confidence_band returns NA outside trim range and valid inside", {
  set.seed(1)
  times <- seq(0.1, 1.0, length.out = 20)
  point <- pmax(0.01, pmin(0.99, sort(runif(20, 0.1, 0.9))))
  boot  <- matrix(rnorm(20 * 200, 0, 0.05), 20, 200) + point

  band <- confidence_band(point, boot, times = times, n = 100)
  expect_named(band, c("ll.band", "ul.band"))

  # At least some entries non-NA (in trim range)
  expect_true(any(!is.na(band$ll.band)))
  # Where defined, lower < point < upper
  ok <- !is.na(band$ll.band)
  expect_true(all(band$ll.band[ok] <= point[ok] + 1e-8))
  expect_true(all(band$ul.band[ok] >= point[ok] - 1e-8))
})

test_that("confidence_band validates dimensions", {
  expect_error(confidence_band(c(0.3, 0.5),
                               matrix(0, 3, 100),
                               times = c(1, 2), n = 50),
               "nrow\\(boot\\)")
})


test_that("ks_pvalue is in [0, 1] and reflects observed statistic", {
  set.seed(2)
  diff_point <- rnorm(50, 0, 0.05)
  diff_boot  <- matrix(rnorm(50 * 500, 0, 0.05), 50, 500)

  ks <- ks_pvalue(diff_point, diff_boot, n = 100)
  expect_named(ks, c("statistic", "p.value"))
  expect_true(ks$p.value >= 0 && ks$p.value <= 1)
  expect_true(ks$statistic >= 0)
})

test_that("ks_pvalue rejects under large observed difference", {
  set.seed(3)
  # Big observed difference, small bootstrap noise. Centered-bootstrap
  # convention: diff_boot replicates the observed statistic, so its
  # mean tracks diff_point (0.5), not 0.
  diff_point <- rep(0.5, 30)
  diff_boot  <- matrix(rnorm(30 * 200, 0.5, 0.01), 30, 200)

  ks <- ks_pvalue(diff_point, diff_boot, n = 100)
  expect_lt(ks$p.value, 0.05)
})

test_that("ks_pvalue does not reject under small observed difference", {
  set.seed(4)
  # Small observed difference, larger bootstrap variability
  diff_point <- rnorm(30, 0, 0.01)
  diff_boot  <- matrix(rnorm(30 * 500, 0, 0.1), 30, 500)

  ks <- ks_pvalue(diff_point, diff_boot, n = 100)
  expect_gt(ks$p.value, 0.05)
})
