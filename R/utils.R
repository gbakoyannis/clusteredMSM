#' Validate Independent-Clusters Assumption
#'
#' Checks that no individual ID appears in more than one cluster.
#' Stops with an informative error if the assumption is violated.
#'
#' @param data A data frame.
#' @param cid Character. Cluster ID column name.
#' @param id Character. Individual ID column name.
#'
#' @keywords internal
check_clusters <- function(data, cid, id) {
  by_id <- stats::aggregate(data[[cid]], by = list(data[[id]]),
                            FUN = function(x) length(unique(x)))
  bad <- by_id$x > 1
  if (any(bad)) {
    stop("same individual(s) appear in more than 1 cluster ",
         "(violation of the independent clusters assumption)")
  }
  invisible(TRUE)
}


#' Add Cluster-Size Column for Inverse-Cluster-Size Weighting
#'
#' Computes the number of distinct individuals per cluster and merges
#' it back into the long-format data as a \code{clust.size} column.
#'
#' @param data A long-format data frame containing \code{cid} and \code{id}.
#' @param cid Character. Cluster ID column name.
#' @param id Character. Individual ID column name.
#'
#' @return The input data with an added \code{clust.size} column.
#'
#' @keywords internal
add_cluster_sizes <- function(data, cid, id) {
  per_subject <- stats::aggregate(rep(1, nrow(data)),
                                  by = list(data[[cid]], data[[id]]),
                                  FUN = mean)
  per_cluster <- stats::aggregate(per_subject$x,
                                  by = list(per_subject$Group.1),
                                  FUN = sum)
  names(per_cluster) <- c(cid, "clust.size")

  data <- merge(data, per_cluster, by = cid)
  data <- data[order(data[[cid]], data[[id]]), ]
  rownames(data) <- NULL
  data
}


#' Interpolate a Step Function onto a Time Grid
#'
#' Right-continuous step interpolation of a (time, value) sequence onto
#' a target grid. Used to align bootstrap replicates onto a common time
#' grid so they can be assembled into a matrix.
#'
#' @param times Numeric. Jump times where the function changes.
#' @param values Numeric. Function values; \code{values[k]} is the
#'   value on \code{[times[k], times[k+1])}.
#' @param grid Numeric. Output grid.
#' @param init Numeric scalar. Value to return for grid points before
#'   \code{times[1]}. Default 0.
#'
#' @return Numeric vector of values on \code{grid}.
#'
#' @keywords internal
step_interp <- function(times, values, grid, init = 0) {
  if (length(times) == 0L) return(rep(init, length(grid)))
  fn <- stats::stepfun(times, c(init, values))
  fn(grid)
}
