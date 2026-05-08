#' Population-Averaged Transition Probabilities for Multistate Data
#'
#' Estimates the working-independence Aalen-Johansen transition
#' probability \eqn{P(X(t) = j \mid X(s) = h)} for clustered or
#' independent multistate process data, with cluster-bootstrap
#' standard errors, pointwise confidence intervals, and (optionally)
#' simultaneous confidence bands. If a grouping variable is supplied
#' on the right-hand side of the formula, also conducts a two-sample
#' Kolmogorov-Smirnov-type test of curve equality.
#'
#' @param formula A formula of the form
#'   \code{msm(Tstart, Tstop, Sstart, Sstop) ~ 1} for a one-sample
#'   analysis, or \code{msm(Tstart, Tstop, Sstart, Sstop) ~ group} for
#'   a two-sample analysis (estimate + test).
#' @param data A data frame containing the variables in \code{formula}
#'   plus the \code{id} (and optional \code{cluster}) columns.
#' @param tmat A K x K transition matrix, typically built with
#'   \code{\link{trans_mat}()}.
#' @param id Character. Name of the subject ID column. Required.
#' @param cluster Character or \code{NA}. Name of the cluster ID column.
#'   When \code{NA} (default), the bootstrap resamples individuals;
#'   when supplied, it resamples whole clusters.
#' @param h Integer in 1..K. Starting state.
#' @param j Integer in 1..K. Ending state.
#' @param s Numeric scalar. Conditioning time. Default 0.
#' @param weighted Logical. Inverse-cluster-size weighting (requires
#'   \code{cluster}). Default \code{FALSE}.
#' @param LMAJ Logical. Use the landmark Aalen-Johansen estimator
#'   (recommended when \code{s > 0} and the Markov assumption is
#'   implausible). Default \code{FALSE}.
#' @param B Integer. Number of bootstrap replications. \code{B = 0}
#'   skips inference and is permitted only for one-sample formulas.
#'   Default 1000.
#' @param cband Logical. If \code{TRUE}, return 95\% simultaneous
#'   confidence band limits (recommended \code{B >= 1000}).
#' @param level Confidence level. Default 0.95.
#' @param seed Optional integer for reproducible bootstrap.
#'
#' @return An S3 object of class \code{patp} containing:
#' \itemize{
#'   \item \code{call}, \code{formula}: the original call and formula.
#'   \item \code{curves}: a data frame with columns \code{time},
#'     \code{P}, \code{group} (if two-sample), and (if \code{B > 0})
#'     \code{se}, \code{ll}, \code{ul} -- and \code{ll.band},
#'     \code{ul.band} if \code{cband = TRUE}.
#'   \item \code{test}: \code{NULL} (one-sample) or a list with the
#'     observed K-S statistic and bootstrap p-value (two-sample).
#'   \item \code{n_subjects}, \code{n_clusters}, \code{groups},
#'     \code{h}, \code{j}, \code{s}, \code{B}.
#' }
#'
#' @references
#' Bakoyannis, G. (2021). Nonparametric analysis of nonhomogeneous
#' multistate processes with clustered observations. \emph{Biometrics},
#' 77(2), 533-546. \doi{10.1111/biom.13327}
#'
#' @examples
#' \dontrun{
#' tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)),
#'                   names = c("Healthy", "Ill", "Dead"))
#'
#' # One-sample
#' fit <- patp(msm(time0, time1, state0, state1) ~ 1,
#'             data = mydata, tmat = tmat,
#'             id = "subj_id", cluster = "site",
#'             h = 1, j = 2, s = 0,
#'             B = 1000, cband = TRUE)
#' fit
#'
#' # Two-sample
#' tt <- patp(msm(time0, time1, state0, state1) ~ treatment,
#'            data = mydata, tmat = tmat,
#'            id = "subj_id", cluster = "site",
#'            h = 1, j = 2, B = 1000)
#' tt
#' }
#'
#' @seealso \code{\link{msm}}, \code{\link{trans_mat}},
#'   \code{\link{validate_intervals}}.
#' @export
patp <- function(formula, data, tmat,
                 id, cluster = NA,
                 h, j, s = 0,
                 weighted = FALSE, LMAJ = FALSE,
                 B = 1000, cband = FALSE,
                 level = 0.95, seed = NULL) {

  call_ <- match.call()

  ## ---- Parse and validate ----
  parsed <- parse_msm_formula(formula, data, id = id, cluster = cluster)
  validate_intervals(parsed$data, tmat)
  long <- intervals_to_long(parsed$data, tmat)

  if (weighted && !parsed$has_cluster) {
    stop("weighted = TRUE requires the 'cluster' argument")
  }

  ## ---- Dispatch based on RHS ----
  is_two_sample <- !is.null(parsed$group)
  if (is_two_sample && B == 0L) {
    stop("two-sample formulas require B > 0")
  }

  fit <- if (is_two_sample) {
    .patp_twosample(long, tmat,
                    has_cluster = parsed$has_cluster,
                    h = h, j = j, s = s,
                    weighted = weighted, LMAJ = LMAJ,
                    B = B, cband = cband, level = level, seed = seed,
                    group_name = parsed$group_name)
  } else {
    .patp_onesample(long, tmat,
                    has_cluster = parsed$has_cluster,
                    h = h, j = j, s = s,
                    weighted = weighted, LMAJ = LMAJ,
                    B = B, cband = cband, level = level, seed = seed)
  }

  fit$call    <- call_
  fit$formula <- formula
  fit$h <- h; fit$j <- j; fit$s <- s; fit$B <- B
  fit$n_subjects <- length(unique(long$id))

  structure(fit, class = c("patp", "list"))
}
