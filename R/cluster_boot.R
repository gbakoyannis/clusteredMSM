#' Generic Cluster Bootstrap Engine
#'
#' Performs the nonparametric cluster bootstrap: resamples whole clusters
#' with replacement and applies a user-supplied statistic function to
#' each bootstrap replicate. Replaces \code{mstate::msboot()} with a
#' function that does not assume a particular form for the statistic.
#'
#' @param data A data frame containing the cluster ID column named by
#'   \code{cid}. All other columns are passed through unchanged to
#'   \code{fn}.
#' @param cid Character. Name of the cluster ID column.
#' @param B Integer. Number of bootstrap replications.
#' @param fn Function. Called as \code{fn(boot_data, ...)} once per
#'   bootstrap replicate. Must return a numeric vector of fixed length
#'   (the length must be the same across replicates -- typically achieved
#'   by evaluating the statistic on a fixed grid of times).
#' @param ... Additional arguments passed to \code{fn}.
#' @param strata Optional named vector mapping each unique cluster ID
#'   (as character) to a stratum label. When supplied, resampling is
#'   stratified: in each replicate, the original number of clusters in
#'   each stratum is drawn with replacement from the clusters in that
#'   stratum (so per-stratum cluster counts are preserved across
#'   replicates). When \code{NULL} (default), unstratified resampling
#'   is used.
#' @param seed Optional integer. If non-NULL, \code{set.seed(seed)} is
#'   called before bootstrapping for reproducibility.
#' @param verbose Logical. If \code{TRUE}, prints a progress message
#'   every 100 replicates. Default \code{FALSE}.
#'
#' @return A numeric matrix with \code{length(fn(data, ...))} rows and
#'   \code{B} columns. Column \code{b} is the result of \code{fn} applied
#'   to bootstrap replicate \code{b}.
#'
#' @details
#' The bootstrap is non-parametric and at the cluster level: in each
#' replicate, n clusters are sampled with replacement (where n is the
#' number of unique values of \code{data[[cid]]}), and the resulting
#' cluster-bound rows are stacked into a new data frame.
#'
#' Resampled clusters are re-IDed (1, 2, ..., n) so that downstream
#' code treating cluster IDs as distinct labels (e.g., \code{tapply}
#' aggregations) works correctly even when the same original cluster
#' appears multiple times in a replicate.
#'
#' Stratified resampling (\code{strata}) supports cluster-randomized
#' designs in which each cluster belongs to a single group (case ii of
#' Bakoyannis & Bandyopadhyay 2022): the per-stratum cluster counts
#' are fixed across replicates, so a replicate can never wipe out one
#' of the groups.
#'
#' The bootstrap is statistic-agnostic: the same engine is used for
#' point-estimator standard errors, two-sample tests, and confidence
#' bands. Whatever \code{fn} returns is what gets bootstrapped.
#'
#' @examples
#' data(example_msm)
#'
#' # Cluster-bootstrap a simple summary: the mean of Tstop across
#' # rows. The user-supplied fn can return any fixed-length numeric
#' # vector; the typical use case is a transition probability curve
#' # evaluated on a fixed time grid (which is exactly what patp() does
#' # under the hood).
#' boot <- cluster_boot(
#'   data = example_msm, cid = "cluster", B = 50,
#'   fn   = function(d) mean(d$Tstop),
#'   seed = 1
#' )
#' c(point = mean(boot, na.rm = TRUE),
#'   se    = stats::sd(boot, na.rm = TRUE))
#'
#' @export
cluster_boot <- function(data, cid, B, fn, ...,
                         strata = NULL, seed = NULL, verbose = FALSE) {

  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!cid %in% names(data)) stop("column '", cid, "' not found in data")
  if (!is.numeric(B) || length(B) != 1L || B < 1L) {
    stop("'B' must be a positive integer")
  }
  if (!is.function(fn)) stop("'fn' must be a function")

  if (!is.null(seed)) set.seed(seed)

  cluster_ids <- unique(data[[cid]])
  n <- length(cluster_ids)

  ## Build the per-replicate cluster sampler. Two modes:
  ##   - unstratified: sample n clusters with replacement from all clusters
  ##   - stratified:   for each stratum, sample its cluster count with
  ##                   replacement from the clusters in that stratum
  if (is.null(strata)) {
    sample_clusters <- function() sample(cluster_ids, n, replace = TRUE)
  } else {
    if (is.null(names(strata))) {
      stop("'strata' must be a named vector (names = cluster IDs)")
    }
    cid_chr <- as.character(cluster_ids)
    missing_in_strata <- setdiff(cid_chr, names(strata))
    if (length(missing_in_strata) > 0L) {
      stop("'strata' is missing entries for cluster IDs: ",
           paste(utils::head(missing_in_strata, 5L), collapse = ", "))
    }
    strata_aligned <- unname(strata[cid_chr])
    by_stratum <- split(cluster_ids, strata_aligned)
    sample_clusters <- function() {
      drawn <- lapply(by_stratum, function(ids)
        sample(ids, length(ids), replace = TRUE))
      unlist(drawn, use.names = FALSE)
    }
  }

  # Pre-split data by cluster for fast subsetting in the loop
  cluster_rows <- split(seq_len(nrow(data)), data[[cid]])

  # Probe call: try the original data first to determine out_length and
  # validate that fn returns numeric. Tolerate failure -- some fn may
  # only succeed on resamples; in that case defer length determination
  # to the first successful replicate.
  probe <- tryCatch(fn(data, ...), error = function(e) NULL)
  out <- NULL
  out_length <- NA_integer_
  if (!is.null(probe)) {
    if (!is.numeric(probe)) {
      stop("'fn' must return a numeric vector")
    }
    out_length <- length(probe)
    out <- matrix(NA_real_, nrow = out_length, ncol = B)
  }

  for (b in seq_len(B)) {
    drawn <- sample_clusters()
    n_drawn <- length(drawn)

    # Build the bootstrap data frame. Re-ID clusters 1..n_drawn so duplicates
    # are kept distinct for downstream cluster-aware operations.
    pieces <- vector("list", n_drawn)
    for (k in seq_len(n_drawn)) {
      rows <- cluster_rows[[as.character(drawn[k])]]
      d <- data[rows, , drop = FALSE]
      d[[cid]] <- k
      pieces[[k]] <- d
    }
    boot_data <- do.call(rbind, pieces)

    result <- tryCatch(fn(boot_data, ...), error = function(e) NULL)
    if (is.null(result)) next
    if (!is.numeric(result)) {
      if (is.null(out)) stop("'fn' must return a numeric vector")
      next
    }
    if (is.null(out)) {
      out_length <- length(result)
      out <- matrix(NA_real_, nrow = out_length, ncol = B)
    }
    if (length(result) != out_length) next
    out[, b] <- result

    if (verbose && b %% 100 == 0) {
      message("cluster_boot: completed ", b, " of ", B, " replicates")
    }
  }

  if (is.null(out)) {
    stop("'fn' failed on every bootstrap replicate and on the probe call")
  }

  out
}
