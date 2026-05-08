#' Build a Transition Matrix for a Multistate Process
#'
#' Constructs a square transition matrix in the format expected by
#' \code{\link{patp}} and \code{\link{prodint_AJ}}. Drop-in replacement for
#' \code{mstate::transMat()} that supports cyclic transitions
#' (e.g., illness-death with recovery).
#'
#' @param x A list of length K. Element \code{k} is an integer vector of
#'   states reachable directly from state \code{k}. Use
#'   \code{integer(0)} or \code{NULL} for absorbing states.
#' @param names Optional character vector of length K with state names.
#'   If \code{NULL}, states are named \code{"1"}, \code{"2"}, ..., \code{"K"}.
#'
#' @return A K x K integer matrix. Entry \code{[h, j]} is the integer
#'   transition ID for the direct transition h -> j, or \code{NA} if no
#'   such transition exists. Transition IDs are assigned in row-major
#'   order (all transitions out of state 1 first, then state 2, etc.).
#'   The matrix has \code{dimnames} \code{list(from = names, to = names)}.
#'
#' @examples
#' # Illness-death without recovery (progressive)
#' trans_mat(list(c(2, 3), 3, integer(0)),
#'           names = c("Healthy", "Ill", "Dead"))
#'
#' # Illness-death WITH recovery (non-monotone)
#' trans_mat(list(c(2, 3), c(1, 3), integer(0)),
#'           names = c("Healthy", "Ill", "Dead"))
#'
#' # Four-state competing risks with recovery from two intermediate states
#' trans_mat(list(c(2, 3, 4), c(1, 4), c(1, 4), integer(0)))
#'
#' @export
trans_mat <- function(x, names = NULL) {

  if (!is.list(x)) stop("'x' must be a list")
  K <- length(x)
  if (K < 2L) stop("'x' must have at least 2 elements (need >= 2 states)")

  if (is.null(names)) {
    names <- as.character(seq_len(K))
  } else {
    if (!is.character(names) || length(names) != K) {
      stop("'names' must be a character vector of length ", K)
    }
    if (anyDuplicated(names)) stop("'names' must be unique")
  }

  M <- matrix(NA_integer_, nrow = K, ncol = K)
  dimnames(M) <- list(from = names, to = names)

  trans_id <- 1L
  for (from in seq_len(K)) {
    targets <- x[[from]]
    if (is.null(targets) || length(targets) == 0L) next

    if (!is.numeric(targets) || any(targets != as.integer(targets))) {
      stop("element ", from, " of 'x' must contain integer state indices")
    }
    targets <- as.integer(targets)

    if (any(targets < 1L | targets > K)) {
      stop("element ", from, " of 'x' contains state index outside 1:", K)
    }
    if (any(targets == from)) {
      stop("self-transition ", from, " -> ", from,
           " is not allowed (states must change)")
    }
    if (anyDuplicated(targets)) {
      stop("element ", from, " of 'x' has duplicated target states")
    }

    for (to in targets) {
      M[from, to] <- trans_id
      trans_id <- trans_id + 1L
    }
  }

  if (trans_id == 1L) stop("no transitions defined")

  M
}
