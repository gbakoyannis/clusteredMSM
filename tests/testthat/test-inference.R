## Helper: build a T x B matrix of bootstrap replicates with given
## per-row mean (point) and sd (se), under a fixed seed.
make_boot <- function(point, se, B = 2000, seed = 1) {
  set.seed(seed)
  T_ <- length(point)
  if (length(se) == 1L) se <- rep(se, T_)
  matrix(stats::rnorm(T_ * B,
                      mean = rep(point, B),
                      sd   = rep(se,    B)),
         nrow = T_, ncol = B)
}

test_that("ci_cloglog returns bounds straddling the point estimate", {
  point <- c(0.2, 0.5, 0.8)
  boot  <- make_boot(point, 0.05)
  ci <- ci_cloglog(point, boot)

  expect_named(ci, c("ll", "ul"))
  expect_true(all(ci$ll < point))
  expect_true(all(ci$ul > point))
  expect_true(all(ci$ll > 0 & ci$ul < 1))
})

test_that("ci_cloglog matches hand-computed cloglog-scale interval", {
  ## Hand-computation for P_hat = 0.1, replicates ~ N(0.1, 0.02^2):
  ##   g(0.1) = log(-log(0.1)) ~= 0.8340
  ##   |g'(0.1)| = 1/(0.1 * |log(0.1)|) ~= 4.3429
  ##   SE_g ~= 4.3429 * 0.02 ~= 0.0869
  ##   95% CI on cloglog scale: 0.8340 +/- 1.96 * 0.0869
  ##   Back-transform: ~ (0.066, 0.143)
  P_hat <- 0.1
  boot  <- make_boot(P_hat, 0.02, B = 1000, seed = 42)

  ci <- ci_cloglog(P_hat, boot, level = 0.95)
  ## Tolerances loose enough to absorb finite-sample (B=1000) noise
  ## around the asymptotic delta-method values 0.066 / 0.143.
  expect_equal(ci$ll, 0.066, tolerance = 0.10)
  expect_equal(ci$ul, 0.143, tolerance = 0.10)
  expect_true(ci$ll < P_hat && P_hat < ci$ul)
  expect_true(ci$ll > 0 && ci$ul < 1)
})

test_that("ci_cloglog handles boundary point estimates with NA", {
  point <- c(0, 0.5, 1)
  boot  <- make_boot(point, 0.05)
  ci <- ci_cloglog(point, boot)

  expect_true(is.na(ci$ll[1]))
  expect_true(is.na(ci$ul[1]))
  expect_false(is.na(ci$ll[2]))
  expect_true(is.na(ci$ll[3]))
})

test_that("ci_cloglog narrows with smaller bootstrap variability", {
  point <- rep(0.5, 3)
  ci_wide   <- ci_cloglog(point, make_boot(point, 0.10, seed = 7))
  ci_narrow <- ci_cloglog(point, make_boot(point, 0.01, seed = 7))

  expect_true(all(ci_narrow$ul - ci_narrow$ll <
                  ci_wide$ul   - ci_wide$ll))
})

test_that("ci_cloglog level argument widens the interval", {
  point <- 0.5
  boot  <- make_boot(point, 0.05, seed = 11)
  ci_95 <- ci_cloglog(point, boot, level = 0.95)
  ci_99 <- ci_cloglog(point, boot, level = 0.99)

  expect_true(ci_99$ul - ci_99$ll > ci_95$ul - ci_95$ll)
})

test_that("ci_cloglog: doubling bootstrap SE roughly doubles cloglog-scale width", {
  point  <- 0.3
  ## Reuse the same standardized noise so the only difference is sd
  set.seed(123)
  z      <- stats::rnorm(2000)
  boot1  <- matrix(point + 0.02 * z, nrow = 1L)
  boot2  <- matrix(point + 0.04 * z, nrow = 1L)

  ci1 <- ci_cloglog(point, boot1)
  ci2 <- ci_cloglog(point, boot2)

  g       <- function(p) log(-log(p))
  width_g <- function(ci) g(ci$ll) - g(ci$ul)  # g is decreasing
  ratio   <- width_g(ci2) / width_g(ci1)
  expect_equal(ratio, 2, tolerance = 0.10)
})

test_that("ci_cloglog: bounds always lie in (0,1) and bracket the point", {
  set.seed(7)
  point <- stats::runif(20, 0.05, 0.95)
  boot  <- make_boot(point, 0.05, B = 500, seed = 17)
  ci    <- ci_cloglog(point, boot)

  expect_true(all(ci$ll > 0 & ci$ll < point))
  expect_true(all(ci$ul > point & ci$ul < 1))
})

test_that("ci_cloglog validates input dimensions", {
  expect_error(ci_cloglog(c(0.3, 0.5), matrix(0, 3, 100)),
               "nrow\\(boot\\)")
  expect_error(ci_cloglog(0.5, c(0.4, 0.5, 0.6)),
               "matrix")
})


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

test_that("confidence_band is wider than pointwise CI on cloglog scale", {
  set.seed(2)
  times <- seq(0.1, 1.0, length.out = 30)
  point <- pmax(0.05, pmin(0.95, sort(stats::runif(30, 0.1, 0.9))))
  boot  <- make_boot(point, 0.05, B = 500, seed = 5)

  ci   <- ci_cloglog(point, boot)
  band <- confidence_band(point, boot, times = times,
                          trim = c(0, 1))

  ok <- !is.na(band$ll.band) & !is.na(ci$ll)
  ## Bands are at least as wide as pointwise intervals at every t in range
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
  diff_boot  <- matrix(stats::rnorm(30 * 200, 0.5, 0.01), 30, 200)

  ks <- ks_pvalue(diff_point, diff_boot, n = 100)
  expect_lt(ks$p.value, 0.05)
})

test_that("ks_pvalue does not reject under small observed difference", {
  set.seed(4)
  # Small observed difference, larger bootstrap variability
  diff_point <- stats::rnorm(30, 0, 0.01)
  diff_boot  <- matrix(stats::rnorm(30 * 500, 0, 0.1), 30, 500)

  ks <- ks_pvalue(diff_point, diff_boot, n = 100)
  expect_gt(ks$p.value, 0.05)
})
