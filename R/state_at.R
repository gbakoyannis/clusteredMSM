#' Identify Each Subject's State at a Given Time
#'
#' Given a long-format multistate dataset, returns the state each subject
#' occupies at time \code{s}. Used by the landmark Aalen-Johansen estimator
#' to identify the risk set at a conditioning time.
#'
#' @param data A data frame in long format with columns \code{Tstart},
#'   \code{Tstop}, \code{from}, \code{to}, \code{status}, and the subject
#'   ID column named by \code{id}.
#' @param s Numeric scalar. The time at which to evaluate states.
#' @param id Character. Name of the subject ID column. Default \code{"id"}.
#'
#' @return A data frame with columns \code{<id>} and \code{state}, one row
#'   per subject who is under observation at time \code{s}. Subjects whose
#'   follow-up has not yet started (\code{min(Tstart) > s}) or who have
#'   already entered an absorbing state (\code{max(Tstop) <= s} with
#'   \code{status == 1}) are excluded.
#'
#' @details
#' A subject's state at time \code{s} is the \code{from} state of the
#' interval that contains \code{s} (i.e., \code{Tstart <= s < Tstop}).
#' Because the long format duplicates intervals across competing
#' transitions, the result is deduplicated to one row per subject.
#'
#' For exact ties (\code{s == Tstop}), the subject is treated as having
#' transitioned: they appear in their new state if a subsequent interval
#' exists, otherwise they are excluded as absorbed/censored.
#'
#' @references
#' Putter H, Spitoni C (2018). Non-parametric estimation of transition
#' probabilities in non-Markov multi-state models: the landmark
#' Aalen-Johansen estimator. \emph{Statistical Methods in Medical
#' Research} 27(7):2081-2092. \doi{10.1177/0962280216674497}
#' 
#' @examples
#' \dontrun{
#' # After running intervals_to_long on illness-death data:
#' state_at(msd, s = 1.5, id = "id")
#' }
#'
#' @export
state_at <- function(data, s, id = "id") {

  required <- c("Tstart", "Tstop", "from", "to", "status", id)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stop("'data' is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(s) || length(s) != 1L) {
    stop("'s' must be a numeric scalar")
  }

  data <- data[order(data[[id]], data$Tstart, data$Tstop), ]

  # Subjects whose follow-up hasn't reached time s
  starts  <- tapply(data$Tstart, data[[id]], min)
  stops   <- tapply(data$Tstop,  data[[id]], max)
  in_obs_window <- names(starts)[starts <= s & stops > s]

  # Soft "absorbed" rule: keep subjects with any interval extending past s.
  has_later <- unique(data[[id]][data$Tstop > s])
  keep_ids <- intersect(in_obs_window, as.character(has_later))

  sub <- data[as.character(data[[id]]) %in% as.character(keep_ids) &
              data$Tstart <= s & data$Tstop > s,
              c(id, "from")]

  out <- unique(sub)
  names(out) <- c(id, "state")
  rownames(out) <- NULL
  out
}
