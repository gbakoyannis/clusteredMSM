#' Fit Cox Model Stratified by Transition and Return Cumulative Hazards
#'
#' Fits a Cox proportional hazards model stratified by transition type and
#' returns the Breslow-estimated cumulative transition hazards in long
#' format.
#'
#' @param data A data frame in long multistate format with columns
#'   \code{Tstart}, \code{Tstop}, \code{status}, and \code{trans}.
#'   Typically the output of \code{\link{intervals_to_long}()}.
#' @param tmat A K x K transition matrix from \code{trans_mat()}. Used
#'   only to validate that all transition IDs in \code{data} are defined
#'   in \code{tmat}; not used in the fit itself.
#' @param weights Optional numeric vector of length \code{nrow(data)}
#'   with observation weights. Pass \code{1 / data$clust.size} for the
#'   inverse-cluster-size-weighted estimator. \code{NULL} (default)
#'   means unweighted.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{time}: jump time of the transition.
#'   \item \code{Haz}: Breslow estimate of the cumulative transition
#'     hazard at \code{time}.
#'   \item \code{trans}: integer transition ID (matches \code{tmat}).
#' }
#' Sorted by \code{trans} then \code{time}.
#'
#' @details
#' The model fitted is
#' \code{coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
#' method = "breslow")}.
#' This is the Andersen-Gill counting-process formulation, which handles
#' multiple intervals per subject -- including non-monotone processes
#' with recovery -- provided the input data was correctly constructed
#' (one row per sojourn x candidate transition).
#'
#' @examples
#' data(example_msm)
#' tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)),
#'                   names = c("Healthy", "Ill", "Dead"))
#'
#' # fit_chaz takes long-format data; intervals_to_long() converts
#' # the interval format users supply. Most users reach fit_chaz
#' # indirectly via patp() -- this example shows the manual pipeline.
#' long <- intervals_to_long(example_msm, tmat)
#' haz  <- fit_chaz(long, tmat)
#' head(haz)
#'
#' @importFrom survival coxph Surv strata basehaz
#' @export
fit_chaz <- function(data, tmat, weights = NULL) {

  required <- c("Tstart", "Tstop", "status", "trans")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stop("'data' is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(data)) {
      stop("'weights' must be a numeric vector of length nrow(data)")
    }
    if (any(weights < 0, na.rm = TRUE)) {
      stop("'weights' must be non-negative")
    }
  }

  # Validate transition IDs against tmat
  if (!missing(tmat) && !is.null(tmat)) {
    valid_ids <- sort(unique(as.vector(tmat[!is.na(tmat)])))
    bad <- setdiff(unique(data$trans), valid_ids)
    if (length(bad)) {
      stop("data$trans contains IDs not in tmat: ",
           paste(bad, collapse = ", "))
    }
  }

  fit <- if (is.null(weights)) {
    survival::coxph(
      survival::Surv(Tstart, Tstop, status) ~ survival::strata(trans),
      data   = data,
      method = "breslow"
    )
  } else {
    survival::coxph(
      survival::Surv(Tstart, Tstop, status) ~ survival::strata(trans),
      data    = data,
      weights = weights,
      method  = "breslow"
    )
  }

  bh <- survival::basehaz(fit, centered = FALSE)

  # survival::basehaz formats strata levels as "trans=1", "trans=2", ...
  # in modern releases; older releases used plain "1", "2", .... Detect
  # which is in use and parse accordingly (avoids a coercion warning
  # from running the wrong cast first).
  raw <- as.character(bh$strata)
  trans <- if (any(grepl("=", raw, fixed = TRUE))) {
    as.integer(sub(".*=", "", raw))
  } else {
    as.integer(raw)
  }

  out <- data.frame(
    time  = bh$time,
    Haz   = bh$hazard,
    trans = trans
  )

  out <- out[order(out$trans, out$time), ]
  rownames(out) <- NULL
  out
}
