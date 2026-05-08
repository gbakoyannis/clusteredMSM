#' Construct a Multistate Interval Object for Use in patp() Formulas
#'
#' Packages four columns -- interval start time, interval end time,
#' starting state, and ending state -- into a single object suitable
#' for use as the response in a \code{\link{patp}} formula.
#' Conceptually analogous to \code{\link[survival]{Surv}()} in survival
#' analysis: the user calls it inside a formula
#' (\code{msm(Tstart, Tstop, Sstart, Sstop) ~ ...}).
#'
#' @param Tstart Numeric vector. Start time of each interval.
#' @param Tstop Numeric vector. End time of each interval.
#' @param Sstart Integer-valued numeric vector. State occupied at
#'   \code{Tstart}. States must be whole numbers; the matrix backing
#'   stores them as double (matrices have a single storage mode), but
#'   values are validated to be integer-valued.
#' @param Sstop Integer-valued numeric vector. State occupied at
#'   \code{Tstop}. Same storage caveat as \code{Sstart}. If
#'   \code{Sstart == Sstop}, the row represents a censored interval;
#'   this is permitted only on the final row of a subject's record.
#'
#' @return A four-column matrix of class \code{"msm"} with columns
#'   \code{Tstart}, \code{Tstop}, \code{Sstart}, \code{Sstop}.
#'
#' @details
#' Used inside a model formula passed to \code{\link{patp}}, e.g.
#' \code{patp(msm(Tstart, Tstop, Sstart, Sstop) ~ 1, ...)}. The four
#' arguments can be named anything in the user's data -- the formula
#' machinery looks them up in the supplied \code{data}.
#'
#' @examples
#' Tstart <- c(0,   1.5, 0,   0)
#' Tstop  <- c(1.5, 3.0, 2.0, 1.2)
#' Sstart <- c(1,   2,   1,   1)
#' Sstop  <- c(2,   3,   3,   1)   # last row censored (S unchanged)
#' obj <- msm(Tstart, Tstop, Sstart, Sstop)
#' head(obj)
#'
#' @seealso \code{\link{patp}}, \code{\link{validate_intervals}}.
#' @export
msm <- function(Tstart, Tstop, Sstart, Sstop) {

  if (missing(Tstart) || missing(Tstop) ||
      missing(Sstart) || missing(Sstop)) {
    stop("msm() requires four arguments: Tstart, Tstop, Sstart, Sstop")
  }

  n <- length(Tstart)
  if (length(Tstop)  != n || length(Sstart) != n || length(Sstop) != n) {
    stop("All arguments to msm() must have the same length")
  }
  if (!is.numeric(Tstart) || !is.numeric(Tstop)) {
    stop("Tstart and Tstop must be numeric")
  }
  if (any(Tstart >= Tstop, na.rm = TRUE)) {
    stop("Tstart must be strictly less than Tstop on every row")
  }
  if (!is.numeric(Sstart) || !is.numeric(Sstop)) {
    stop("Sstart and Sstop must be numeric (integer-valued)")
  }
  if (any(Sstart != as.integer(Sstart), na.rm = TRUE) ||
      any(Sstop  != as.integer(Sstop),  na.rm = TRUE)) {
    stop("Sstart and Sstop must be integer-valued (whole numbers)")
  }

  out <- cbind(Tstart = Tstart, Tstop = Tstop,
               Sstart = Sstart, Sstop = Sstop)
  class(out) <- c("msm", class(out))
  out
}


#' @export
print.msm <- function(x, ...) {
  cat("msm object: ", nrow(x), " interval",
      if (nrow(x) != 1L) "s" else "", "\n", sep = "")
  print(unclass(utils::head(x, 6L)), ...)
  if (nrow(x) > 6L) cat("... (", nrow(x) - 6L, " more)\n", sep = "")
  invisible(x)
}
