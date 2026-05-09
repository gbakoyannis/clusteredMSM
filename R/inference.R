#' Pointwise Confidence Intervals via Complementary Log-Log Transformation
#'
#' Computes pointwise confidence intervals for a transition probability
#' estimate by bootstrapping directly on the cloglog scale: each
#' bootstrap replicate is transformed by \eqn{g(p) = \log(-\log p)}, the
#' standard error is taken as the across-replicate \code{sd} of the
#' transformed values, a symmetric interval is built on the cloglog
#' scale, and the result is back-transformed to the probability scale.
#'
#' @param point Numeric vector of point estimates in (0, 1), length T.
#' @param boot  Numeric matrix of bootstrap replicates of \code{point}:
#'   \code{nrow(boot) == length(point)} (rows index time points,
#'   columns index replicates), as produced by \code{cluster_boot()}.
#' @param level Confidence level. Default 0.95.
#'
#' @return A list with elements \code{ll} and \code{ul}: numeric vectors
#'   of lower and upper confidence limits, the same length as
#'   \code{point}.
#'
#' @details
#' Bootstrapping on the cloglog scale avoids the delta-method
#' approximation that an earlier implementation used (and which had
#' incorrect coverage when the bootstrap SE was passed in already on
#' the probability scale). The cloglog transformation is monotone
#' \emph{decreasing} on (0, 1), so the upper end of the interval on
#' the cloglog scale corresponds to the \emph{lower} end on the
#' probability scale, and vice versa.
#'
#' For point estimates exactly equal to 0 or 1 the cloglog is
#' undefined, and \code{NA} is returned for both limits at those
#' positions. Bootstrap replicates outside (0, 1) are masked to
#' \code{NA} when computing the cloglog-scale standard error.
#'
#' @keywords internal
#' @export
ci_cloglog <- function(point, boot, level = 0.95) {

  if (!is.matrix(boot)) {
    stop("'boot' must be a numeric matrix (rows = time points, cols = replicates)")
  }
  if (nrow(boot) != length(point)) {
    stop("nrow(boot) must equal length(point)")
  }

  z <- stats::qnorm(1 - (1 - level) / 2)

  ll <- ul <- rep(NA_real_, length(point))
  ok <- !is.na(point) & point > 0 & point < 1

  if (!any(ok)) return(list(ll = ll, ul = ul))

  g_boot <- boot[ok, , drop = FALSE]
  invalid <- !is.finite(g_boot) | g_boot <= 0 | g_boot >= 1
  g_boot[invalid] <- NA_real_
  g_boot <- log(-log(g_boot))

  g_point <- log(-log(point[ok]))
  se_g    <- apply(g_boot, 1L, stats::sd, na.rm = TRUE)

  # g(p) = log(-log p) is monotone decreasing on (0, 1):
  # the upper end on the cloglog scale gives the lower end on the
  # probability scale after back-transformation by exp(-exp(.)).
  ll[ok] <- exp(-exp(g_point + z * se_g))
  ul[ok] <- exp(-exp(g_point - z * se_g))

  list(ll = ll, ul = ul)
}


#' Simultaneous Confidence Band via Cluster-Bootstrap Quantile
#'
#' Constructs a simultaneous (uniform-over-time) confidence band for a
#' transition probability curve by working on the cloglog scale: each
#' bootstrap replicate is transformed, the studentized supremum
#' \eqn{\sup_t |g(\hat P_b(t)) - g(\hat P(t))| / \mathrm{SE}_g(t)} is
#' computed for every replicate, the \code{level} quantile of those
#' suprema is taken as the critical value, and the resulting band is
#' back-transformed to the probability scale.
#'
#' @param point Numeric vector of point estimates over a time grid.
#' @param boot Numeric matrix of bootstrap replicates: rows correspond
#'   to time points (matching \code{point}), columns to bootstrap
#'   replicates.
#' @param times Numeric vector of times matching the rows of \code{boot}
#'   and the entries of \code{point}.
#' @param level Confidence level. Default 0.95.
#' @param trim Numeric vector of length 2. Lower and upper quantiles of
#'   the jump-time distribution at which to clip the band (avoids the
#'   "fans out" behavior near the extremes of follow-up). Default
#'   \code{c(0.05, 0.95)}.
#'
#' @return A list with elements \code{ll.band} and \code{ul.band}:
#'   numeric vectors of band limits, with \code{NA} outside the trimmed
#'   range or where \code{point} is at the boundary 0/1.
#'
#' @details
#' The bootstrap SD on the cloglog scale already estimates the
#' standard error of \eqn{g(\hat P(t))} directly, so no further
#' \eqn{\sqrt{n}} rescaling is applied. The studentized-supremum
#' construction yields a band that is automatically variance-adapted
#' over time.
#'
#' @keywords internal
#' @export
confidence_band <- function(point, boot, times,
                            level = 0.95, trim = c(0.05, 0.95)) {

  if (!is.matrix(boot)) {
    stop("'boot' must be a numeric matrix (rows = time points, cols = replicates)")
  }
  if (nrow(boot) != length(point)) {
    stop("nrow(boot) must equal length(point)")
  }
  if (length(times) != length(point)) {
    stop("'times' must have the same length as 'point'")
  }

  na_band <- function() list(ll.band = rep(NA_real_, length(point)),
                             ul.band = rep(NA_real_, length(point)))

  jump_times <- times[c(TRUE, diff(point) != 0)]
  if (length(jump_times) == 0L) return(na_band())
  qs <- stats::quantile(jump_times, probs = trim)

  in_range <- !is.na(point) & point > 0 & point < 1 &
              times >= qs[1] & times <= qs[2]
  if (!any(in_range)) return(na_band())

  g_point <- log(-log(point[in_range]))
  g_boot  <- boot[in_range, , drop = FALSE]
  invalid <- !is.finite(g_boot) | g_boot <= 0 | g_boot >= 1
  g_boot[invalid] <- NA_real_
  g_boot <- log(-log(g_boot))

  se_g <- apply(g_boot, 1L, stats::sd, na.rm = TRUE)

  # Studentized residuals on cloglog scale, replicate-by-replicate sup
  resid <- abs(sweep(g_boot, 1L, g_point, FUN = "-")) / se_g
  sups  <- apply(resid, 2L, max, na.rm = TRUE)
  sups  <- sups[is.finite(sups)]
  if (length(sups) == 0L) return(na_band())
  c_a <- stats::quantile(sups, probs = level, na.rm = TRUE)

  ll.band <- ul.band <- rep(NA_real_, length(point))
  # g is decreasing => upper-cloglog = lower-probability
  ll.band[in_range] <- exp(-exp(g_point + c_a * se_g))
  ul.band[in_range] <- exp(-exp(g_point - c_a * se_g))

  list(ll.band = ll.band, ul.band = ul.band)
}


#' Two-Sample Kolmogorov-Smirnov-Type Test via Cluster Bootstrap
#'
#' Computes a p-value for the null hypothesis that a transition
#' probability curve is equal across two groups, based on the supremum
#' of the absolute difference between group-specific Aalen-Johansen
#' estimators and a cluster-bootstrap reference distribution.
#'
#' @param diff_point Numeric vector. Observed group difference of the
#'   point estimator (\eqn{\hat P_1 - \hat P_0}) on a fixed time grid.
#' @param diff_boot Numeric matrix. Cluster bootstrap replicates of the
#'   group difference: rows = time points (matching \code{diff_point}),
#'   columns = replicates.
#' @param n Integer. Number of independent clusters in the pooled sample.
#'
#' @return A list with elements \code{statistic} (the observed
#'   sup-statistic) and \code{p.value}.
#'
#' @details
#' The test statistic is \eqn{T_n = \sqrt{n} \sup_t | \hat P_1(t) -
#' \hat P_0(t) |}, with the null distribution approximated by the
#' bootstrap analogue applied to centered bootstrap differences.
#'
#' @keywords internal
#' @export
ks_pvalue <- function(diff_point, diff_boot, n) {

  if (nrow(diff_boot) != length(diff_point)) {
    stop("nrow(diff_boot) must equal length(diff_point)")
  }

  T_obs <- sqrt(n) * max(abs(diff_point), na.rm = TRUE)

  # Bootstrap analogue: center each bootstrap diff around the observed
  # diff (Bakoyannis 2020 Section 3.2)
  centered <- sweep(diff_boot, 1, diff_point, FUN = "-")
  T_boot   <- sqrt(n) * apply(abs(centered), 2, max, na.rm = TRUE)

  list(statistic = T_obs,
       p.value   = mean(T_boot >= T_obs, na.rm = TRUE))
}
