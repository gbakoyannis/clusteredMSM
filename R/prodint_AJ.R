#' Aalen-Johansen Product Integral for General Multistate Processes
#'
#' Computes the row of the transition probability matrix
#' \deqn{P(s, t) = \prod_{u \in (s, t]} (I + dA(u))}
#' starting from state \code{h}, given Nelson-Aalen / Breslow cumulative
#' transition hazards. Drop-in replacement for \code{mstate::probtrans()}
#' that supports general (possibly non-monotone) transition structures
#' such as illness-death with recovery.
#'
#' @param haz_df Data frame with columns \code{time}, \code{Haz},
#'   \code{trans}, typically the output of \code{fit_chaz()} or the
#'   \code{$Haz} element of an \code{mstate::msfit} object. \code{Haz}
#'   is the cumulative (not incremental) hazard.
#' @param tmat A K x K transition matrix from \code{trans_mat()}.
#'   Cyclic transitions (both \code{[h, j]} and \code{[j, h]} non-NA)
#'   are permitted.
#' @param predt Numeric scalar. Starting time s. Probabilities are
#'   computed for t > predt. Default 0.
#' @param h Integer in 1..K. Starting state.
#'
#' @return A data frame with columns \code{time}, \code{pstate1}, ...,
#'   \code{pstateK}. Row k gives \eqn{P(predt, time_k)[h, ]}, the
#'   probabilities of being in each state at \code{time_k} given the
#'   process was in state \code{h} at \code{predt}.
#'
#' @details
#' At each unique jump time \eqn{u}, the increment matrix
#' \eqn{dA(u)} has off-diagonal entry \eqn{(h, j)} equal to the
#' Nelson-Aalen increment of the cumulative hazard for transition
#' \eqn{h \to j}, and diagonal entries equal to the negative row sums.
#' The transition probability matrix is then updated by
#' \eqn{P \leftarrow P (I + dA(u))}.
#'
#' Unlike \code{mstate::probtrans()}, which uses a forward recursion
#' that exploits acyclicity, this implementation works for any
#' transition structure including processes with cycles. For progressive
#' (acyclic) models, output matches \code{probtrans()} to numerical
#' precision.
#'
#' Computational complexity is O(J K^2) where J is the number of unique
#' jump times and K is the number of states.
#'
#' @examples
#' \dontrun{
#' # Illness-death with recovery
#' tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)),
#'                   names = c("Healthy", "Ill", "Dead"))
#' haz  <- fit_chaz(msd, tmat)
#' P    <- prodint_AJ(haz, tmat, predt = 0, h = 1)
#' head(P)
#' }
#'
#' @export
prodint_AJ <- function(haz_df, tmat, predt = 0, h = 1) {

  ## ---- input checks ----
  stopifnot(is.data.frame(haz_df),
            all(c("time", "Haz", "trans") %in% names(haz_df)))
  stopifnot(is.matrix(tmat), nrow(tmat) == ncol(tmat))
  K <- nrow(tmat)
  stopifnot(length(h) == 1L, h %in% seq_len(K))
  stopifnot(length(predt) == 1L, is.numeric(predt))

  ntrans <- suppressWarnings(max(tmat, na.rm = TRUE))
  if (!is.finite(ntrans)) stop("tmat has no defined transitions")
  if (!all(sort(unique(haz_df$trans)) %in% seq_len(ntrans))) {
    stop("haz_df$trans contains ids not present in tmat")
  }

  ## map each transition id to its (from, to) pair
  idx <- which(!is.na(tmat), arr.ind = TRUE)
  ord <- order(tmat[!is.na(tmat)])
  idx <- idx[ord, , drop = FALSE]
  rownames(idx) <- NULL

  ## ---- restrict to t > predt and find unique jump times ----
  haz_df <- haz_df[haz_df$time > predt, , drop = FALSE]
  haz_df <- haz_df[order(haz_df$trans, haz_df$time), , drop = FALSE]

  jt <- sort(unique(haz_df$time))
  out <- matrix(0, length(jt) + 1L, K)
  P   <- diag(K)
  out[1L, ] <- P[h, ]

  ## previous cumulative hazard per transition (for differencing)
  prev_haz <- numeric(ntrans)

  ## ---- main loop: product integrate over jump times ----
  for (k in seq_along(jt)) {
    rows <- haz_df[haz_df$time == jt[k], , drop = FALSE]
    dA <- matrix(0, K, K)

    for (i in seq_len(nrow(rows))) {
      tr  <- rows$trans[i]
      inc <- rows$Haz[i] - prev_haz[tr]
      if (inc < 0) {
        if (inc < -1e-10) {
          warning("negative hazard increment at t=", jt[k],
                  " trans=", tr)
        }
        inc <- 0
      }
      prev_haz[tr] <- rows$Haz[i]
      dA[idx[tr, 1L], idx[tr, 2L]] <- inc
    }
    diag(dA) <- -rowSums(dA)

    P <- P %*% (diag(K) + dA)
    out[k + 1L, ] <- P[h, ]
  }

  ## ---- assemble output ----
  res <- data.frame(time = c(predt, jt), out)
  names(res)[-1L] <- paste0("pstate", seq_len(K))
  res
}
