#' Strictly Validate User-Supplied Interval Data
#'
#' Performs strict validation on a data frame of multistate interval
#' data (the output of \code{\link{parse_msm_formula}()}). Stops with
#' an informative error on the first violation.
#'
#' @param data A data frame with columns \code{id}, \code{Tstart},
#'   \code{Tstop}, \code{Sstart}, \code{Sstop}.
#' @param tmat A K x K transition matrix from \code{\link{trans_mat}()}.
#' @param tol Numeric tolerance for temporal-contiguity checks.
#'   Default 1e-9.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass. Stops with
#'   an error otherwise.
#'
#' @details
#' Validation rules, applied in order:
#' \enumerate{
#'   \item Required columns are present and of correct type.
#'   \item Every row has \code{Tstart < Tstop}.
#'   \item All states are integers in \code{1..K}.
#'   \item Within each subject, rows are sorted by time.
#'   \item Within each subject, rows are temporally contiguous:
#'     \code{Tstop[k] == Tstart[k+1]}.
#'   \item Within each subject, rows are spatially contiguous:
#'     \code{Sstop[k] == Sstart[k+1]}.
#'   \item No row has \code{Sstart} equal to an absorbing state.
#'   \item No row follows one whose \code{Sstop} is absorbing.
#'   \item Each non-censored transition (\code{Sstart != Sstop}) is
#'     allowed by \code{tmat}.
#'   \item \code{Sstart == Sstop} (censoring) is permitted only on the
#'     final row of a subject's record.
#' }
#'
#' @keywords internal
#' @export
validate_intervals <- function(data, tmat, tol = 1e-9) {

  required <- c("id", "Tstart", "Tstop", "Sstart", "Sstop")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stop("'data' is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(data$Tstart) || !is.numeric(data$Tstop)) {
    stop("'Tstart' and 'Tstop' must be numeric")
  }
  if (!is.numeric(data$Sstart) || !is.numeric(data$Sstop)) {
    stop("'Sstart' and 'Sstop' must be integer-valued numerics")
  }
  if (any(data$Sstart != as.integer(data$Sstart)) ||
      any(data$Sstop  != as.integer(data$Sstop))) {
    stop("'Sstart' and 'Sstop' must be integer-valued")
  }

  K <- nrow(tmat)
  if (any(data$Sstart < 1L | data$Sstart > K) ||
      any(data$Sstop  < 1L | data$Sstop  > K)) {
    stop("State values must be integers in 1..", K)
  }

  if (any(data$Tstart >= data$Tstop)) {
    bad <- which(data$Tstart >= data$Tstop)[1L]
    stop("Tstart must be < Tstop on every row; first offender is row ",
         bad)
  }

  absorbing <- which(rowSums(!is.na(tmat)) == 0L)

  ## ---- Per-subject contiguity and absorbing-state rules ----
  data <- data[order(data$id, data$Tstart), , drop = FALSE]
  ids <- data$id
  starts <- which(!duplicated(ids))
  stops  <- c(starts[-1L] - 1L, length(ids))

  for (s in seq_along(starts)) {
    rows <- starts[s]:stops[s]
    sub  <- data[rows, , drop = FALSE]
    n    <- nrow(sub)
    sid  <- sub$id[1L]

    # Sstart not absorbing on any row
    if (any(sub$Sstart %in% absorbing)) {
      stop("subject '", sid, "': Sstart is an absorbing state on row ",
           which(sub$Sstart %in% absorbing)[1L], " (this should be ",
           "the previous row's Sstop only)")
    }

    if (n >= 2L) {
      # Temporal contiguity
      gaps <- abs(sub$Tstop[-n] - sub$Tstart[-1L])
      if (any(gaps > tol)) {
        bad <- which(gaps > tol)[1L]
        stop("subject '", sid, "': intervals are not temporally ",
             "contiguous (Tstop[", bad, "] = ", sub$Tstop[bad],
             ", Tstart[", bad + 1L, "] = ", sub$Tstart[bad + 1L], ")")
      }
      # Spatial contiguity
      mismatches <- which(sub$Sstop[-n] != sub$Sstart[-1L])
      if (length(mismatches)) {
        bad <- mismatches[1L]
        stop("subject '", sid, "': states are not spatially contiguous ",
             "(Sstop[", bad, "] = ", sub$Sstop[bad],
             ", Sstart[", bad + 1L, "] = ", sub$Sstart[bad + 1L], ")")
      }
      # No row after one whose Sstop is absorbing
      if (any(sub$Sstop[-n] %in% absorbing)) {
        stop("subject '", sid, "': has a row after entering absorbing ",
             "state ", sub$Sstop[-n][sub$Sstop[-n] %in% absorbing][1L])
      }
    }

    # Censoring (Sstart == Sstop) only on the final row
    self <- which(sub$Sstart == sub$Sstop)
    if (length(self)) {
      if (any(self != n)) {
        stop("subject '", sid, "': Sstart == Sstop is permitted only ",
             "on the last row (censoring); offender is row ",
             self[self != n][1L])
      }
    }

    # Non-censored transitions must be allowed by tmat
    transitioning <- sub$Sstart != sub$Sstop
    if (any(transitioning)) {
      tr_rows <- which(transitioning)
      for (r in tr_rows) {
        if (is.na(tmat[sub$Sstart[r], sub$Sstop[r]])) {
          stop("subject '", sid, "': transition ", sub$Sstart[r],
               " -> ", sub$Sstop[r], " on row ", r,
               " is not allowed by tmat")
        }
      }
    }
  }

  invisible(TRUE)
}


#' Convert Validated Interval Data to the Internal Long Format
#'
#' Internal helper. After \code{\link{validate_intervals}()} has
#' confirmed a dataset is well-formed, this function expands each
#' interval row into one row per allowed transition out of the
#' interval's state, with \code{status = 1} only on the row whose
#' destination matches the actual transition.
#'
#' @param data A validated interval data frame (output of
#'   \code{\link{parse_msm_formula}()}).
#' @param tmat A K x K transition matrix.
#'
#' @return A data frame in long multistate format with columns
#'   \code{id}, \code{Tstart}, \code{Tstop}, \code{from}, \code{to},
#'   \code{trans}, \code{status}, plus \code{cluster} and \code{group}
#'   if present.
#'
#' @keywords internal
intervals_to_long <- function(data, tmat) {

  K <- nrow(tmat)
  pieces <- vector("list", nrow(data))

  for (i in seq_len(nrow(data))) {
    row_i <- data[i, , drop = FALSE]
    h <- row_i$Sstart
    censored <- row_i$Sstart == row_i$Sstop

    out_idx <- which(!is.na(tmat[h, ]))
    if (length(out_idx) == 0L) next  # absorbing -- shouldn't happen post-validation

    block <- data.frame(
      id     = rep(row_i$id,     length(out_idx)),
      Tstart = rep(row_i$Tstart, length(out_idx)),
      Tstop  = rep(row_i$Tstop,  length(out_idx)),
      from   = rep(h,            length(out_idx)),
      to     = out_idx,
      trans  = tmat[h, out_idx],
      status = as.integer(!censored & out_idx == row_i$Sstop)
    )

    if (!is.null(row_i$cluster)) block$cluster <- row_i$cluster
    if (!is.null(row_i$group))   block$group   <- row_i$group

    pieces[[i]] <- block
  }

  out <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1L))])
  rownames(out) <- NULL
  structure(out, class = c("cmsdata", "data.frame"))
}
