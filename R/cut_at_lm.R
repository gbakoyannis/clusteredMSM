#' Truncate a Long-Format Multistate Dataset at a Landmark Time
#'
#' Restricts a long-format multistate dataset to follow-up after a
#' landmark time \code{s}. Intervals ending before \code{s} are removed;
#' intervals straddling \code{s} are truncated to start at \code{s}.
#' Drop-in replacement for \code{mstate::cutLMms()}.
#'
#' @param data A data frame in long format with columns \code{Tstart},
#'   \code{Tstop}, and \code{status}. Typically the output of
#'   \code{ms_prep()} or \code{mstate::msprep()}.
#' @param s Numeric scalar. The landmark time.
#'
#' @return A data frame of the same shape as \code{data}, restricted and
#'   truncated as described above.
#'
#' @details
#' This function is used in landmark Aalen-Johansen estimation, which
#' relaxes the Markov assumption by re-fitting hazards on the cohort
#' still under observation at time \code{s}.
#'
#' Important: this function only truncates the time axis. To compute the
#' landmark estimator, you must also restrict to subjects who are in the
#' relevant starting state at time \code{s} -- use \code{state_at()} for
#' that step.
#'
#' Intervals where \code{Tstop == s} exactly are dropped (the subject has
#' already transitioned by \code{s}); intervals where \code{Tstart >= s}
#' are kept unchanged; intervals straddling \code{s}
#' (\code{Tstart < s < Tstop}) have their \code{Tstart} reset to \code{s}.
#' If a straddling interval was going to end in an event
#' (\code{status == 1}) at \code{Tstop > s}, the event is preserved.
#'
#' @examples
#' \dontrun{
#' msd_lm <- cut_at_lm(msd, s = 1.5)
#' }
#'
#' @export
cut_at_lm <- function(data, s) {

  required <- c("Tstart", "Tstop", "status")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stop("'data' is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(s) || length(s) != 1L) {
    stop("'s' must be a numeric scalar")
  }

  # Drop intervals that end at or before s
  data <- data[data$Tstop > s, , drop = FALSE]

  # Truncate intervals that straddle s
  straddling <- data$Tstart < s
  data$Tstart[straddling] <- s

  rownames(data) <- NULL
  data
}
