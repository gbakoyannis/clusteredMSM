#' Print Method for patp Objects
#'
#' Concise summary of a \code{\link{patp}} fit. For two-sample fits,
#' includes the K-S statistic and p-value.
#'
#' @param x A \code{patp} object.
#' @param digits Integer. Number of significant digits for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return The object \code{x}, invisibly.
#'
#' @export
print.patp <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nclusteredMSM: Population-Averaged Transition Probability\n")
  cat(strrep("-", 56), "\n", sep = "")

  cat("Call:\n  ")
  print(x$call)

  cat("\nEstimand: P(X(t) =", x$j, "| X(", x$s, ") =", x$h, ")\n")
  cat("Subjects:", x$n_subjects)
  if (!is.na(x$n_clusters)) {
    cat("   Clusters:", x$n_clusters, "\n")
  } else {
    cat("   (no clustering)\n")
  }

  if (x$B > 0L) {
    cat("Bootstrap replications: B =", x$B, "\n")
  } else {
    cat("Bootstrap: NONE (B = 0; standard errors and CIs not computed)\n")
  }

  if (!is.null(x$test)) {
    cat("\nTwo-sample test (", x$test$type, "):\n", sep = "")
    cat("  Groups:    ", paste(x$groups, collapse = " vs "), "\n", sep = "")
    cat("  Statistic: ", format(x$test$statistic, digits = digits), "\n",
        sep = "")
    cat("  p-value:   ", format.pval(x$test$p.value, digits = digits), "\n",
        sep = "")
  }

  cat("\nCurve summary (", nrow(x$curves), " rows total):\n", sep = "")
  print(utils::head(x$curves, 6L), digits = digits, row.names = FALSE)
  if (nrow(x$curves) > 6L) {
    cat("... (", nrow(x$curves) - 6L, " more rows; access via $curves)\n",
        sep = "")
  }

  invisible(x)
}


#' Summary Method for patp Objects
#'
#' Returns the full curve(s) and (if applicable) the test result.
#'
#' @param object A \code{patp} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with named components for further inspection.
#'
#' @export
summary.patp <- function(object, ...) {
  out <- list(
    call         = object$call,
    formula      = object$formula,
    estimand     = sprintf("P(X(t) = %d | X(%g) = %d)",
                           object$j, object$s, object$h),
    n_subjects   = object$n_subjects,
    n_clusters   = object$n_clusters,
    B            = object$B,
    curves       = object$curves,
    test         = object$test,
    groups       = object$groups
  )
  class(out) <- "summary.patp"
  out
}


#' @export
print.summary.patp <- function(x, ...) {
  cat("\nSummary of patp fit\n")
  cat(strrep("-", 56), "\n", sep = "")
  cat("Call:\n  "); print(x$call)
  cat("\nEstimand:  ", x$estimand, "\n", sep = "")
  cat("Subjects:  ", x$n_subjects, "\n", sep = "")
  if (!is.na(x$n_clusters)) cat("Clusters:  ", x$n_clusters, "\n", sep = "")
  cat("Bootstrap: B = ", x$B, "\n", sep = "")
  if (!is.null(x$test)) {
    cat("\nTwo-sample test:\n")
    cat("  Type:      ", x$test$type, "\n")
    cat("  Statistic: ", format(x$test$statistic), "\n")
    cat("  p-value:   ", format.pval(x$test$p.value), "\n")
    cat("  Groups:    ", paste(x$groups, collapse = " vs "), "\n")
  }
  cat("\nFull curve(s):\n")
  print(x$curves, row.names = FALSE)
  invisible(x)
}
