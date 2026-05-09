test_that("ci_cloglog returns bounds straddling the point estimate", {
  point <- c(0.2, 0.5, 0.8)
  se    <- c(0.05, 0.05, 0.05)
  ci <- ci_cloglog(point, se)

  expect_named(ci, c("ll", "ul"))
  expect_true(all(ci$ll < point))
  expect_true(all(ci$ul > point))
  expect_true(all(ci$ll > 0 & ci$ul < 1))
})

test_that("ci_cloglog matches hand-computed delta-method interval", {
  ## P_hat = 0.1, SE on probability scale = 0.02:
  ##   g(0.1)        ~= 0.8340
  ##   |g'(0.1)|     = 1/(0.1 * |log 0.1|) ~= 4.3429
  ##   se_g          ~= 4.3429 * 0.02 ~= 0.0869
  ##   95% CI(cloglog) = 0.8340 +/- 1.96 * 0.0869
  ##   Back-transform  -> (0.0652, 0.1434)
  ci <- ci_cloglog(point = 0.1, se = 0.02, level = 0.95)
  expect_equal(ci$ll, 0.0652, tolerance = 1e-3)
  expect_equal(ci$ul, 0.1434, tolerance = 1e-3)
  expect_true(ci$ll < 0.1 && 0.1 < ci$ul)
  expect_true(ci$ll > 0  && ci$ul < 1)
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
  ci_wide   <- ci_cloglog(point, rep(0.10, 3))
  ci_narrow <- ci_cloglog(point, rep(0.01, 3))

  expect_true(all(ci_narrow$ul - ci_narrow$ll <
                  ci_wide$ul   - ci_wide$ll))
})

test_that("ci_cloglog level argument widens the interval", {
  point <- 0.5
  se    <- 0.05
  ci_95 <- ci_cloglog(point, se, level = 0.95)
  ci_99 <- ci_cloglog(point, se, level = 0.99)

  expect_true(ci_99$ul - ci_99$ll > ci_95$ul - ci_95$ll)
})

test_that("ci_cloglog: doubling SE roughly doubles cloglog-scale width", {
  point <- 0.3
  ci1 <- ci_cloglog(point, se = 0.02)
  ci2 <- ci_cloglog(point, se = 0.04)

  g       <- function(p) log(-log(p))
  width_g <- function(ci) g(ci$ll) - g(ci$ul)  # g is decreasing
  expect_equal(width_g(ci2) / width_g(ci1), 2, tolerance = 1e-12)
})

test_that("ci_cloglog: bounds in (0,1) and bracket the point", {
  set.seed(7)
  point <- stats::runif(20, 0.05, 0.95)
  se    <- rep(0.05, 20)
  ci    <- ci_cloglog(point, se)

  expect_true(all(ci$ll > 0 & ci$ll < point))
  expect_true(all(ci$ul > point & ci$ul < 1))
})

test_that("ci_cloglog validates input lengths", {
  expect_error(ci_cloglog(c(0.3, 0.5), c(0.1)), "same length")
})


## Helper for confidence_band tests: T x B replicate matrix.
make_boot <- function(point, se, B = 500, seed = 1) {
  set.seed(seed)
  T_ <- length(point)
  if (length(se) == 1L) se <- rep(se, T_)
  matrix(stats::rnorm(T_ * B,
                      mean = rep(point, B),
                      sd   = rep(se,    B)),
         nrow = T_, ncol = B)
}

test_that("confidence_band returns NA outside trim range and valid inside", {
  set.seed(1)
  times <- seq(0.1, 1.0, length.out = 20)
  point <- pmax(0.01, pmin(0.99, sort(stats::runif(20, 0.1, 0.9))))
  boot  <- make_boot(point, 0.05, B = 200, seed = 13)

  band <- confidence_band(point, boot, times = times)
  expect_named(band, c("ll.band", "ul.band"))

  expect_true(any(!is.na(band$ll.band)))
  ok <- !is.na(band$ll.band)
  expect_true(all(band$ll.band[ok] <= point[ok] + 1e-8))
  expect_true(all(band$ul.band[ok] >= point[ok] - 1e-8))
})

test_that("confidence_band is at least as wide as pointwise CI", {
  set.seed(2)
  times <- seq(0.1, 1.0, length.out = 30)
  point <- pmax(0.05, pmin(0.95, sort(stats::runif(30, 0.1, 0.9))))
  boot  <- make_boot(point, 0.05, B = 500, seed = 5)
  se    <- apply(boot, 1L, stats::sd, na.rm = TRUE)

  ci   <- ci_cloglog(point, se)
  band <- confidence_band(point, boot, times = times,
                          trim = c(0, 1))

  ok <- !is.na(band$ll.band) & !is.na(ci$ll)
  expect_true(all(band$ll.band[ok] <= ci$ll[ok] + 1e-8))
  expect_true(all(band$ul.band[ok] >= ci$ul[ok] - 1e-8))
})

test_that("confidence_band validates dimensions", {
  expect_error(confidence_band(c(0.3, 0.5),
                               matrix(0, 3, 100),
                               times = c(1, 2)),
               "nrow\\(boot\\)")
})


test_that("ks_pvalue is in [0, 1] and reflects observed statistic", {
  set.seed(2)
  diff_point <- stats::rnorm(50, 0, 0.05)
  diff_boot  <- matrix(stats::rnorm(50 * 500, 0, 0.05), 50, 500)

  ks <- ks_pvalue(diff_point, diff_boot, scale = sqrt(100))
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
  diff_boot  <- matrix(stats::rnorm(30 * 200, 0.5, 0.01), 30, 200)

  ks <- ks_pvalue(diff_point, diff_boot, scale = sqrt(100))
  expect_lt(ks$p.value, 0.05)
})

test_that("ks_pvalue does not reject under small observed difference", {
  set.seed(4)
  # Small observed difference, larger bootstrap variability
  diff_point <- stats::rnorm(30, 0, 0.01)
  diff_boot  <- matrix(stats::rnorm(30 * 500, 0, 0.1), 30, 500)

  ks <- ks_pvalue(diff_point, diff_boot, scale = sqrt(100))
  expect_gt(ks$p.value, 0.05)
})
