#' Pointwise Confidence Intervals via Complementary Log-Log Transformation
#'
#' Computes pointwise confidence intervals for a transition probability
#' estimate using the cloglog transformation \eqn{g(p) = \log(-\log p)}
#' with a delta-method standard error scaling. The bootstrap is run
#' once, on the original probability scale, in line with Bakoyannis
#' (2021) Theorem 2; the cloglog SE is then derived analytically.
#'
#' @param point Numeric vector of point estimates in (0, 1), length T.
#' @param se Numeric vector of bootstrap standard errors of \code{point}
#'   on the probability scale, same length as \code{point} (i.e.
#'   \code{apply(boot_matrix, 1, sd, na.rm = TRUE)}).
#' @param level Confidence level. Default 0.95.
#'
#' @return A list with elements \code{ll} and \code{ul}: numeric vectors
#'   of lower and upper confidence limits, the same length as
#'   \code{point}.
#'
#' @details
#' By the delta method,
#' \eqn{\mathrm{SE}(g(\hat P)) = \mathrm{SE}(\hat P) / |\hat P \log \hat P|}
#' (here \eqn{|g'(p)| = 1/|p \log p|} on (0,1)). The interval is built
#' symmetrically on the cloglog scale and back-transformed via
#' \eqn{p = \exp(-\exp(\cdot))}. Because \eqn{g} is monotone
#' decreasing on (0, 1), the upper end on the cloglog scale becomes
#' the \emph{lower} end on the probability scale, and vice versa.
#'
#' For point estimates exactly equal to 0 or 1 the cloglog is
#' undefined and \code{NA} is returned at those positions.
#'
#' @keywords internal
#' @export
ci_cloglog <- function(point, se, level = 0.95) {

  if (length(point) != length(se)) {
    stop("'point' and 'se' must have the same length")
  }

  z <- stats::qnorm(1 - (1 - level) / 2)

  ll <- ul <- rep(NA_real_, length(point))
  ok <- !is.na(point) & point > 0 & point < 1 &
        !is.na(se)   & is.finite(se)

  if (!any(ok)) return(list(ll = ll, ul = ul))

  P    <- point[ok]
  g_P  <- log(-log(P))
  se_g <- se[ok] / abs(P * log(P))   # delta-method scaling, |g'(p)| = 1/|p log p|

  # g is decreasing on (0, 1):
  # upper end on cloglog scale -> lower end on probability scale.
  ll[ok] <- exp(-exp(g_P + z * se_g))
  ul[ok] <- exp(-exp(g_P - z * se_g))

  list(ll = ll, ul = ul)
}


#' Simultaneous Confidence Band via Cluster-Bootstrap Quantile
#'
#' Constructs a simultaneous (uniform-over-time) confidence band for a
#' transition probability curve. The band is built on the cloglog
#' scale: for each bootstrap replicate, the studentized supremum
#' \eqn{\sup_t |g(\hat P_b(t)) - g(\hat P(t))| / \mathrm{SE}_g(t)} is
#' computed, where \eqn{\mathrm{SE}_g(t)} is the delta-method SE on
#' the cloglog scale. The \code{level} quantile of those suprema is
#' the critical value, and the band is back-transformed to the
#' probability scale.
#'
#' @param point Numeric vector of point estimates over a time grid.
#' @param boot Numeric matrix of bootstrap replicates of \code{point}
#'   on the original probability scale: rows correspond to time
#'   points (matching \code{point}), columns to bootstrap replicates.
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
#' Only one bootstrap on the original (probability) scale is required:
#' \eqn{\mathrm{SE}_g(t)} is obtained by the delta method from
#' \eqn{\mathrm{SE}(\hat P(t)) = \mathrm{sd}(\hat P^*(t))}, and the
#' cloglog values \eqn{g(\hat P_b(t))} are computed by transforming
#' the existing replicates pointwise -- no separate cloglog-scale
#' bootstrap is run.
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

  P    <- point[in_range]
  B_in <- boot[in_range, , drop = FALSE]

  # Probability-scale SE of P_hat, then delta-method cloglog SE.
  se   <- apply(B_in, 1L, stats::sd, na.rm = TRUE)
  se_g <- se / abs(P * log(P))

  # Cloglog of P_hat and of replicates; mask invalid replicate values.
  g_P     <- log(-log(P))
  invalid <- !is.finite(B_in) | B_in <= 0 | B_in >= 1
  G       <- B_in
  G[invalid] <- NA_real_
  G       <- log(-log(G))

  # Studentized residuals on cloglog scale; replicate-by-replicate sup
  resid <- abs(sweep(G, 1L, g_P, FUN = "-")) / se_g
  sups  <- apply(resid, 2L, max, na.rm = TRUE)
  sups  <- sups[is.finite(sups)]
  if (length(sups) == 0L) return(na_band())
  c_a <- stats::quantile(sups, probs = level, na.rm = TRUE)

  ll.band <- ul.band <- rep(NA_real_, length(point))
  # g is decreasing => upper-cloglog = lower-probability
  ll.band[in_range] <- exp(-exp(g_P + c_a * se_g))
  ul.band[in_range] <- exp(-exp(g_P - c_a * se_g))

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
#' @param scale Numeric scalar. Asymptotic scaling factor applied to
#'   both the observed sup-statistic and the bootstrap analogue.
#'   Use \eqn{\sqrt{n}} (\eqn{n} = total clusters) for the
#'   dependent-groups design (Bakoyannis 2021, Theorem 3); use
#'   \eqn{\sqrt{n_1 n_2 / (n_1 + n_2)}} for the cluster-randomized
#'   design (Bakoyannis & Bandyopadhyay 2022, Theorem 2).
#'
#' @return A list with elements \code{statistic} (the observed
#'   sup-statistic, on the chosen scale) and \code{p.value}.
#'
#' @details
#' The test statistic is \eqn{T = c_n \sup_t | \hat P_1(t) -
#' \hat P_0(t) |} with \eqn{c_n =} \code{scale}, and the null
#' distribution is approximated by the bootstrap analogue applied to
#' centered bootstrap differences. The empirical p-value is invariant
#' to \code{scale} (it appears on both sides of the comparison); the
#' role of \code{scale} is to put the reported statistic on the
#' correct asymptotic scale per the relevant theorem.
#'
#' @keywords internal
#' @export
ks_pvalue <- function(diff_point, diff_boot, scale) {

  if (nrow(diff_boot) != length(diff_point)) {
    stop("nrow(diff_boot) must equal length(diff_point)")
  }

  T_obs <- scale * max(abs(diff_point), na.rm = TRUE)

  # Bootstrap analogue: center each bootstrap diff around the observed
  # diff (Bakoyannis & Bandyopadhyay 2022, Section 3).
  centered <- sweep(diff_boot, 1, diff_point, FUN = "-")
  T_boot   <- scale * apply(abs(centered), 2, max, na.rm = TRUE)

  list(statistic = T_obs,
       p.value   = mean(T_boot >= T_obs, na.rm = TRUE))
}
