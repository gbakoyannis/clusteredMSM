#' Fit Cox Model Stratified by Transition and Return Cumulative Hazards
#'
#' Fits a Cox proportional hazards model stratified by transition type and
#' returns the Breslow-estimated cumulative transition hazards in long
#' format. Replaces the combination of \code{mstate::msfit()} and the
#' \code{coxph() + basehaz()} sequence used in the original
#' \code{clustered-multistate} code.
#'
#' @param data A data frame in long multistate format with columns
#'   \code{Tstart}, \code{Tstop}, \code{status}, and \code{trans}.
#'   Typically the output of \code{ms_prep()} or \code{mstate::msprep()}.
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
#' Sorted by \code{trans} then \code{time}. Format is identical to the
#' \code{$Haz} element of an \code{mstate::msfit} object, so it can be
#' passed directly to \code{prodint_AJ()}.
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
#' The Breslow tie-breaking method is used to match the convention of
#' \code{mstate::msfit()} and ensure backward compatibility with the
#' original implementation.
#'
#' @examples
#' \dontrun{
#' tmat <- trans_mat(list(c(2, 3), 3, integer(0)),
#'                   names = c("Healthy", "Ill", "Dead"))
#' haz <- fit_chaz(msd, tmat)
#' head(haz)
#' }
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
