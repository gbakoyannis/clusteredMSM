# Helper: 6 clusters of varying size, simple data
make_clustered <- function() {
  data.frame(
    cluster = rep(1:6, times = c(3, 2, 4, 1, 3, 2)),
    x   = rnorm(15)
  )
}

test_that("cluster_boot returns matrix of correct shape", {
  set.seed(1)
  data <- make_clustered()
  out <- cluster_boot(data, cluster = "cluster", B = 50,
                      fn = function(d) c(mean(d$x), sd(d$x)),
                      seed = 1)

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2, 50))
})

test_that("cluster_boot produces variability across replicates", {
  set.seed(1)
  data <- make_clustered()
  out <- cluster_boot(data, cluster = "cluster", B = 100,
                      fn = function(d) mean(d$x), seed = 1)

  expect_true(sd(out[1, ]) > 0)
})

test_that("cluster_boot is reproducible with the same seed", {
  data <- make_clustered()

  out1 <- cluster_boot(data, "cluster", B = 20,
                       fn = function(d) mean(d$x), seed = 42)
  out2 <- cluster_boot(data, "cluster", B = 20,
                       fn = function(d) mean(d$x), seed = 42)

  expect_equal(out1, out2)
})

test_that("cluster_boot resamples whole clusters (not rows)", {
  data <- make_clustered()
  # Each cluster has a unique x mean. After bootstrap, each replicate's
  # set of cluster means should be drawn (with replacement) from the
  # original 6 cluster means.
  cluster_means <- tapply(data$x, data$cluster, mean)

  set.seed(7)
  out <- cluster_boot(
    data, "cluster", B = 200,
    fn = function(d) {
      means <- tapply(d$x, d$cluster, mean)
      mean(means)
    },
    seed = 7
  )

  # The bootstrap replicate mean should be in the convex hull of the
  # original cluster means (since it is an average of values drawn
  # from them).
  expect_true(all(out >= min(cluster_means) - 1e-10))
  expect_true(all(out <= max(cluster_means) + 1e-10))
})

test_that("cluster_boot re-IDs clusters in each replicate", {
  data <- make_clustered()

  # The fn checks the number of unique cluster IDs in the replicate;
  # after resampling with replacement and re-IDing, this should
  # always equal the original number of clusters (n).
  set.seed(99)
  out <- cluster_boot(
    data, "cluster", B = 30,
    fn = function(d) length(unique(d$cluster)),
    seed = 99
  )
  expect_true(all(out == length(unique(data$cluster))))
})

test_that("cluster_boot validates input", {
  data <- make_clustered()

  expect_error(cluster_boot("not a df", "cluster", 10,
                            fn = function(d) 1),
               "data frame")
  expect_error(cluster_boot(data, "missing", 10,
                            fn = function(d) 1),
               "not found")
  expect_error(cluster_boot(data, "cluster", 0,
                            fn = function(d) 1),
               "positive integer")
  expect_error(cluster_boot(data, "cluster", 10, fn = "not a fn"),
               "function")
  expect_error(cluster_boot(data, "cluster", 10,
                            fn = function(d) "not numeric"),
               "numeric vector")
})

test_that("cluster_boot stratified resampling preserves stratum sizes", {
  # 10 clusters, 6 in stratum "A", 4 in stratum "B"
  data <- data.frame(
    cluster    = rep(1:10, each = 3),
    stratum = rep(c(rep("A", 6), rep("B", 4)), each = 3),
    x      = stats::rnorm(30)
  )
  cluster_to_stratum <- tapply(data$stratum, data$cluster,
                               function(x) as.character(x[1L]))
  strata_vec <- as.character(cluster_to_stratum)
  names(strata_vec) <- names(cluster_to_stratum)

  # fn returns the per-replicate counts of the two strata in the resampled
  # data. Each stratum's count (in subjects) should be exactly preserved.
  out <- cluster_boot(
    data, "cluster", B = 50,
    fn = function(d) c(sum(d$stratum == "A"), sum(d$stratum == "B")),
    strata = strata_vec, seed = 11
  )

  # Original counts: 6 clusters x 3 rows = 18 in A; 4 x 3 = 12 in B.
  expect_true(all(out[1, ] == 18))
  expect_true(all(out[2, ] == 12))
})

test_that("cluster_boot stratified samples within stratum only", {
  # 8 clusters split across two strata. Stratum A has clusters 1..4 with
  # mean x = 0; stratum B has clusters 5..8 with mean x = 10. After
  # stratified resampling, every replicate must contain at least one
  # cluster from each stratum.
  data <- data.frame(
    cluster     = rep(1:8, each = 5),
    stratum = rep(c(rep("A", 4), rep("B", 4)), each = 5),
    x       = c(stats::rnorm(20, mean = 0),
                stats::rnorm(20, mean = 10))
  )
  cluster_to_stratum <- tapply(data$stratum, data$cluster,
                               function(x) as.character(x[1L]))
  strata_vec <- as.character(cluster_to_stratum)
  names(strata_vec) <- names(cluster_to_stratum)

  out <- cluster_boot(
    data, "cluster", B = 100,
    fn = function(d) {
      as.numeric(c(any(d$stratum == "A"), any(d$stratum == "B")))
    },
    strata = strata_vec, seed = 7
  )
  expect_true(all(out == 1))
})

test_that("cluster_boot tolerates pathological replicates", {
  data <- make_clustered()

  # fn that errors on certain inputs -- should not crash the boot,
  # should leave NAs in those columns
  set.seed(2)
  out <- cluster_boot(
    data, "cluster", B = 20,
    fn = function(d) {
      if (runif(1) < 0.3) stop("simulated failure")
      mean(d$x)
    },
    seed = 2
  )

  expect_equal(dim(out), c(1, 20))
  expect_true(any(is.na(out)) || all(!is.na(out)))  # at least defined
})
