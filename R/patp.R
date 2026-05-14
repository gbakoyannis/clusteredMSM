#' Population-Averaged Transition Probabilities for Multistate Process Data
#'
#' Computes the working-independence Aalen-Johansen estimator of the
#' transition probability \eqn{P(X(t) = j \mid X(s) = h)} for clustered or
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
#' @param weighted Logical. Inverse-cluster-size weighting to account for 
#'  potentially informative cluster size (requires \code{cluster}). 
#'  Default \code{FALSE}.
#' @param LMAJ Logical. Use the landmark Aalen-Johansen estimator
#'   (recommended when \code{s > 0} and the Markov assumption is
#'   questionable). Default \code{FALSE}.
#' @param B Integer. Number of bootstrap replications. \code{B = 0}
#'   skips inference and is permitted only for one-sample formulas.
#'   Default 1000.
#' @param cband Logical. If \code{TRUE}, return simultaneous
#'   confidence band limits at the level set by \code{level}
#'   (recommended \code{B >= 1000}).
#' @param design Character, two-sample only. One of
#'   \code{"auto"} (default), \code{"shared"},
#'   \code{"cluster_random"}, \code{"indep_random"}. Selects the
#'   two-sample regime; see the \emph{Two-sample designs} section.
#'   \code{"auto"} infers the regime from the cluster/group
#'   structure of the data, defaulting to \code{"indep_random"} when
#'   each cluster carries a single group (a safer default than
#'   assuming cluster randomization). When that fallback fires, a warning
#'   is emitted asking the user to set \code{"cluster_random"}
#'   explicitly if the per-group cluster counts were fixed by cluster
#'   randomization.
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
#' @section Two-sample designs:
#' Three regimes are supported, each with its own bootstrap scheme
#' and asymptotic scaling factor; the test statistic in every case is
#' \eqn{T = c_n \sup_t |\hat P_1(t) - \hat P_0(t)|}.
#' \describe{
#'   \item{\code{"shared"} (case i): dependent groups (multicenter
#'     trial).}{Every cluster contributes observations from both
#'     groups -- e.g., each hospital randomizes its patients to
#'     treatment or control. Unstratified cluster bootstrap;
#'     \eqn{c_n = \sqrt{n}} (Bakoyannis 2021, Theorem 3), with
#'     \eqn{n} the total number of clusters.}
#'   \item{\code{"cluster_random"} (case ii.a): cluster-randomized
#'     trial.}{Each cluster is assigned to exactly one group, with
#'     per-group cluster counts \eqn{n_1, n_2} fixed by the
#'     randomization. Stratified cluster bootstrap (resample within
#'     each group); \eqn{c_n = \sqrt{n_1 n_2 / (n_1 + n_2)}}
#'     (Bakoyannis & Bandyopadhyay 2022, Theorem 2). This regime
#'     must be opted into; \code{"auto"} will not select it.}
#'   \item{\code{"indep_random"} (case ii.b): independent
#'     observational comparison.}{Each cluster carries one group,
#'     but \eqn{n_1, n_2} are random (the population, not the
#'     analyst, decided who landed in which group). Unstratified
#'     cluster bootstrap; \eqn{c_n = \sqrt{n_1 n_2 / (n_1 + n_2)}}.
#'     The asymptotic regime is the same two-independent-samples
#'     limit as \code{"cluster_random"}; only the bootstrap
#'     differs. This is the \code{"auto"} default when each cluster
#'     carries a single group.}
#' }
#' Mixed structures (some clusters carry both groups, some only one)
#' are not supported in this version.
#'
#' \code{patp()} validates the supplied data against the requested
#' design and stops with an informative error if they disagree.
#'
#' The K-S statistic in v0.1 uses unit weight at every time point.
#' Weighted variants -- following the harmonic-mean weight
#' \eqn{W(t) = \prod_p Y_p(t) / \sum_p Y_p(t)} of Bakoyannis (2021)
#' Section 2.5 and the related construction in Bakoyannis & Bandyopadhyay
#' (2022) -- are planned for v0.2.
#'
#' @section Standard errors and confidence intervals:
#' The cluster bootstrap is run \emph{once}, on the original
#' probability scale, in line with Bakoyannis (2021) Theorem 2. The
#' \code{se} column in \code{curves} is the bootstrap standard
#' deviation of \eqn{\hat P(t)} on that scale, i.e.
#' \code{apply(boot_matrix, 1, sd, na.rm = TRUE)}.
#'
#' Pointwise CIs (\code{ll}, \code{ul}) are then built on the cloglog
#' scale \eqn{g(p) = \log(-\log p)} using the delta-method
#' standardization \eqn{\mathrm{SE}_g(t) = \mathrm{SE}(\hat P(t)) /
#' |\hat P(t) \log \hat P(t)|}, and back-transformed via
#' \eqn{p = \exp(-\exp(\cdot))}. Simultaneous bands
#' (\code{ll.band}, \code{ul.band}) use the same \eqn{\mathrm{SE}_g};
#' the \code{level}-quantile of the supremum gives the critical value.
#'
#' Because the cloglog transformation is nonlinear, \code{se} and the
#' resulting CI widths are \emph{not} equal in general -- the
#' CI is asymmetric on the probability scale. Report \code{se} for
#' descriptive purposes; use \code{ll}/\code{ul} for inference.
#'
#' The simultaneous band is trimmed to the central 90\% of observed
#' jump times by default to avoid the band fanning out near the
#' extremes of follow-up; outside that window \code{ll.band} and
#' \code{ul.band} are \code{NA}. See \code{\link{confidence_band}()}
#' for the \code{trim} argument.
#'
#' Note: in this version the simultaneous band uses an
#' equal-precision-type construction on the cloglog scale, which is
#' asymptotically valid but differs from the Hall-Wellner-type
#' construction of Bakoyannis (2021) Section 2.3 (planned for a
#' future release). See \code{\link{confidence_band}()} for the full
#' disclosure.
#'
#' @references
#' Bakoyannis, G. (2021). Nonparametric analysis of nonhomogeneous
#' multistate processes with clustered observations.
#' \emph{Biometrics}, 77(2), 533-546. \doi{10.1111/biom.13327}
#'
#' Bakoyannis, G., & Bandyopadhyay, D. (2022). Nonparametric tests
#' for multistate processes with clustered data. \emph{Annals of the
#' Institute of Statistical Mathematics}, 74(5), 837-867.
#' \doi{10.1007/s10463-021-00819-x}
#' 
#' Putter H, Spitoni C (2018). Non-parametric estimation of transition
#' probabilities in non-Markov multi-state models: the landmark
#' Aalen-Johansen estimator. \emph{Statistical Methods in Medical
#' Research} 27(7):2081-2092. \doi{10.1177/0962280216674497}
#'
#' @examples
#' data(example_msm)
#' tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)),
#'                   names = c("Healthy", "Ill", "Dead"))
#'
#' # One-sample: P(Ill at t | Healthy at 0). B is small here for speed;
#' # use B >= 1000 for reported results.
#' fit <- patp(msm(Tstart, Tstop, Sstart, Sstop) ~ 1,
#'             data = example_msm, tmat = tmat,
#'             id = "id", cluster = "cluster",
#'             h = 1, j = 2, s = 0,
#'             B = 50, seed = 1)
#' fit
#'
#' # Two-sample (estimate + test in one call). Each cluster in
#' # example_msm carries both treatment levels, so design = "shared".
#' tt <- patp(msm(Tstart, Tstop, Sstart, Sstop) ~ treatment,
#'            data = example_msm, tmat = tmat,
#'            id = "id", cluster = "cluster",
#'            h = 1, j = 2, B = 50,
#'            design = "shared", seed = 1)
#' tt
#'
#' @seealso \code{\link{msm}}, \code{\link{trans_mat}},
#'   \code{\link{validate_intervals}}.
#' @export
patp <- function(formula, data, tmat,
                 id, cluster = NA,
                 h, j, s = 0,
                 weighted = FALSE, LMAJ = FALSE,
                 B = 1000, cband = FALSE,
                 design = c("auto", "shared",
                            "cluster_random", "indep_random"),
                 level = 0.95, seed = NULL) {

  design <- match.arg(design)

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

  if (!is_two_sample && design != "auto") {
    stop("'design' applies only to two-sample formulas")
  }

  fit <- if (is_two_sample) {
    .patp_twosample(long, tmat,
                    has_cluster = parsed$has_cluster,
                    h = h, j = j, s = s,
                    weighted = weighted, LMAJ = LMAJ,
                    B = B, cband = cband, level = level, seed = seed,
                    group_name = parsed$group_name,
                    design = design)
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
