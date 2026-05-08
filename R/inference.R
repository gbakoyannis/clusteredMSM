#' Pointwise Confidence Intervals via Complementary Log-Log Transformation
#'
#' Computes 95% (or arbitrary level) pointwise confidence intervals for a
#' transition probability estimate by applying the cloglog transformation
#' \eqn{g(p) = \log(-\log p)} on the standard error scale and inverting
#' back to the probability scale. Used internally by \code{patp()}.
#'
#' @param point Numeric vector of point estimates in (0, 1).
#' @param se Numeric vector of standard errors, same length as \code{point}.
#' @param level Confidence level. Default 0.95.
#'
#' @return A list with elements \code{ll} and \code{ul}: numeric vectors
#'   of lower and upper confidence limits.
#'
#' @details
#' The cloglog transformation is preferred over the identity for
#' transition probabilities because it (a) keeps confidence limits in
#' (0, 1) without ad hoc truncation and (b) gives more accurate coverage
#' for probabilities near the boundaries.
#'
#' For point estimates exactly equal to 0 or 1, the cloglog is undefined
#' and \code{NA} is returned for both limits at those positions.
#'
#' @keywords internal
#' @export
ci_cloglog <- function(point, se, level = 0.95) {

  if (length(point) != length(se)) {
    stop("'point' and 'se' must have the same length")
  }

  z <- stats::qnorm(1 - (1 - level) / 2)

  ll <- ul <- rep(NA_real_, length(point))
  ok <- point > 0 & point < 1 & !is.na(point) & !is.na(se)

  log_log <- log(-log(point[ok]))
  scale   <- se[ok] / (point[ok] * log(point[ok]))

  ll[ok] <- exp(-exp(log_log - z * scale))
  ul[ok] <- exp(-exp(log_log + z * scale))

  list(ll = ll, ul = ul)
}


#' Simultaneous Confidence Band via Cluster-Bootstrap Quantile
#'
#' Constructs a simultaneous (uniform-over-time) confidence band for a
#' transition probability curve, using the equal-precision-style
#' construction of Bakoyannis (2021) on the cluster bootstrap output.
#'
#' @param point Numeric vector of point estimates over a time grid.
#' @param boot Numeric matrix of bootstrap replicates: rows correspond
#'   to time points (matching \code{point}), columns to bootstrap
#'   replicates.
#' @param times Numeric vector of times matching the rows of \code{boot}
#'   and the entries of \code{point}.
#' @param n Integer. The number of independent clusters in the original
#'   sample. Required for the band scaling.
#' @param level Confidence level. Default 0.95.
#' @param trim Numeric vector of length 2. Lower and upper quantiles of
#'   the jump-time distribution at which to clip the band (avoids the
#'   "fans out" behavior near the extremes of follow-up). Default
#'   \code{c(0.05, 0.95)}.
#'
#' @return A list with elements \code{ll.band} and \code{ul.band}:
#'   numeric vectors of band limits, with \code{NA} outside the trimmed
#'   range.
#'
#' @details
#' The band is computed on the cloglog scale using a variance-stabilizing
#' weight \eqn{q(t) = 1 / (1 + \sigma^2(t))}, then back-transformed.
#' This matches the construction in the original \code{patp()} code.
#'
#' @keywords internal
#' @export
confidence_band <- function(point, boot, times, n,
                            level = 0.95, trim = c(0.05, 0.95)) {

  if (nrow(boot) != length(point)) {
    stop("nrow(boot) must equal length(point)")
  }
  if (length(times) != length(point)) {
    stop("'times' must have the same length as 'point'")
  }

  sigma <- apply(boot, 1, stats::sd, na.rm = TRUE)
  q_t   <- 1 / (1 + sigma^2)

  jump_times <- times[c(TRUE, diff(point) != 0)]
  if (length(jump_times) == 0L) {
    return(list(ll.band = rep(NA_real_, length(point)),
                ul.band = rep(NA_real_, length(point))))
  }
  qs <- stats::quantile(jump_times, probs = trim)

  in_range <- times >= qs[1] & times <= qs[2] & point > 0 & point < 1

  if (sum(in_range) == 0L) {
    return(list(ll.band = rep(NA_real_, length(point)),
                ul.band = rep(NA_real_, length(point))))
  }

  B_t <- q_t[in_range] * boot[in_range, , drop = FALSE] /
           (log(point[in_range]) * point[in_range])
  B_t <- abs(B_t)
  c_a <- stats::quantile(apply(B_t, 2, max, na.rm = TRUE),
                         probs = level, na.rm = TRUE)

  ll.band <- ul.band <- rep(NA_real_, length(point))
  log_log <- log(-log(point[in_range]))
  shift   <- c_a / (sqrt(n) * q_t[in_range])

  ll.band[in_range] <- exp(-exp(log_log + shift))
  ul.band[in_range] <- exp(-exp(log_log - shift))

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
