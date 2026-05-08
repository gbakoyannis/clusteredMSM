#' Parse a patp Formula into Validated Interval Data
#'
#' Internal helper. Takes a formula of the form
#' \code{msm(Tstart, Tstop, Sstart, Sstop) ~ 1} or
#' \code{msm(Tstart, Tstop, Sstart, Sstop) ~ group_var}, plus a data
#' frame and ID column names, and returns a list with the assembled
#' interval data and an optional grouping vector.
#'
#' @param formula A formula. LHS must be a call to \code{msm()}; RHS
#'   must be either \code{1} (one-sample) or a single variable name
#'   (two-sample).
#' @param data A data frame to evaluate the formula in.
#' @param id Character. Name of the subject ID column in \code{data}.
#'   Required.
#' @param cluster Character or \code{NA}. Name of the cluster ID column
#'   in \code{data}. \code{NA} means no clustering (individual-level
#'   bootstrap).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{data}: a data frame with columns \code{id}, \code{Tstart},
#'     \code{Tstop}, \code{Sstart}, \code{Sstop}, plus \code{cluster} if
#'     supplied and \code{group} if the formula has an RHS.
#'   \item \code{group}: \code{NULL} (one-sample) or a vector aligned with
#'     the rows of \code{data} (two-sample).
#'   \item \code{has_cluster}: logical, whether a cluster column was
#'     supplied.
#' }
#'
#' @keywords internal
parse_msm_formula <- function(formula, data, id, cluster = NA) {

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (missing(id) || is.null(id) || !is.character(id) || length(id) != 1L) {
    stop("'id' must be a single character string naming a column in 'data'")
  }
  if (!id %in% names(data)) {
    stop("id column '", id, "' not found in 'data'")
  }

  cluster_supplied <- !is.null(cluster) && !is.na(cluster)
  if (cluster_supplied) {
    if (!is.character(cluster) || length(cluster) != 1L) {
      stop("'cluster' must be a single character string or NA")
    }
    if (!cluster %in% names(data)) {
      stop("cluster column '", cluster, "' not found in 'data'")
    }
  }

  ## ---- Validate the LHS is a call to msm() ----
  lhs <- formula[[2L]]
  if (!is.call(lhs) || as.character(lhs[[1L]]) != "msm") {
    stop("formula's left-hand side must be a call to msm(), e.g. ",
         "msm(Tstart, Tstop, Sstart, Sstop) ~ 1")
  }
  if (length(lhs) != 5L) {  # msm + 4 args
    stop("msm() must be called with exactly 4 arguments: ",
         "Tstart, Tstop, Sstart, Sstop")
  }

  ## ---- Evaluate LHS in data ----
  msm_obj <- eval(lhs, envir = data, enclos = environment(formula))
  if (!inherits(msm_obj, "msm")) {
    stop("LHS evaluation did not produce an msm object")
  }
  if (nrow(msm_obj) != nrow(data)) {
    stop("msm object has ", nrow(msm_obj), " rows but 'data' has ",
         nrow(data), " rows")
  }

  ## ---- Examine RHS for grouping variable ----
  term_labels <- attr(stats::terms(formula), "term.labels")
  if (length(term_labels) > 1L) {
    stop("patp supports a single grouping variable on the right-hand ",
         "side; got ", length(term_labels))
  }

  group <- NULL
  if (length(term_labels) == 1L) {
    if (!term_labels %in% names(data)) {
      stop("grouping variable '", term_labels, "' not found in 'data'")
    }
    group <- data[[term_labels]]
    if (length(unique(stats::na.omit(group))) != 2L) {
      stop("grouping variable '", term_labels, "' must take exactly 2 ",
           "distinct non-NA values; got ",
           length(unique(stats::na.omit(group))))
    }
  }

  ## ---- Assemble output data frame ----
  out <- data.frame(
    id     = data[[id]],
    Tstart = msm_obj[, "Tstart"],
    Tstop  = msm_obj[, "Tstop"],
    Sstart = msm_obj[, "Sstart"],
    Sstop  = msm_obj[, "Sstop"]
  )

  if (cluster_supplied) out$cluster <- data[[cluster]]
  if (!is.null(group))  out$group   <- group

  list(
    data        = out,
    group       = group,
    has_cluster = cluster_supplied,
    group_name  = if (length(term_labels) == 1L) term_labels else NULL
  )
}
