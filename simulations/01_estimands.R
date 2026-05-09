#' True estimand values for the Bakoyannis (2021) simulation scenario.
#'
#' To compute empirical coverage and bias, we need the true values of
#' the population-averaged transition and state occupation probabilities
#' under the simulator's data-generating mechanism. Because the model is
#' non-Markov (frailty + ICS), there is no neat closed form -- but we
#' can compute the truth to arbitrary precision via a *very large*
#' Monte Carlo simulation that does not invoke the package's estimator.
#'
#' Rationale: the true population-averaged probability is
#'   P_{0,2}(t) = E_v[P(X(t) = 2 | v)]
#' where the inner expectation is over the data-generating mechanism
#' conditional on the frailty. We approximate by simulating a very
#' large dataset (no censoring, no truncation) and computing the
#' empirical proportion of subjects in state 2 at time t.
#'
#' For the ACM population, larger clusters are over-represented (each
#' subject contributes one observation, weighted naturally by m_i).
#' For the TCM population, each cluster contributes equally.
#'
#' Run this once and cache the result; it does not depend on (n, F_M).

compute_truth <- function(t_grid,
                          n_clusters_truth = 5000,
                          cluster_size_range = c(5, 15),
                          seed = 20260508) {

  ## Simulate a very large clean dataset (no censoring, no truncation)
  ## Using the same simulator with right_censor = FALSE
  d <- simulate_clusters(
    n_clusters         = n_clusters_truth,
    cluster_size_range = cluster_size_range,
    two_sample         = FALSE,
    right_censor       = FALSE,
    left_truncate      = FALSE,
    tau                = max(t_grid) + 1,
    seed               = seed
  )

  ## For each subject and each time in t_grid, compute their state at t.
  ## A subject's state at time t is the state of the interval containing t,
  ## or the absorbing state if t is past their last interval.
  ids       <- unique(d$id)
  clusters  <- d$cluster[match(ids, d$id)]    # cluster ID per subject

  ## State at t for each subject across grid
  state_mat <- matrix(NA_integer_, nrow = length(ids), ncol = length(t_grid))

  for (k in seq_along(ids)) {
    sub <- d[d$id == ids[k], , drop = FALSE]
    sub <- sub[order(sub$Tstart), , drop = FALSE]

    for (g in seq_along(t_grid)) {
      t <- t_grid[g]
      ## Find interval containing t
      idx <- which(sub$Tstart <= t & sub$Tstop > t)
      if (length(idx) > 0L) {
        state_mat[k, g] <- sub$Sstart[idx[1]]
      } else if (max(sub$Tstop) <= t) {
        ## t is past last interval; subject is in the final Sstop state
        state_mat[k, g] <- sub$Sstop[nrow(sub)]
      }
      ## Else (t before first Tstart) leave NA
    }
  }

  ## ACM-population P(X(t) = 2): proportion of subjects in state 2 at t
  p_acm <- apply(state_mat, 2L, function(col) mean(col == 2L, na.rm = TRUE))

  ## TCM-population: each cluster contributes equally. Average within
  ## cluster first, then across clusters.
  cluster_means <- matrix(NA_real_, nrow = length(unique(clusters)),
                          ncol = length(t_grid))
  unique_clusters <- unique(clusters)
  for (g in seq_along(t_grid)) {
    in_state_2 <- as.integer(state_mat[, g] == 2L)
    for (c_idx in seq_along(unique_clusters)) {
      cluster_subjs <- which(clusters == unique_clusters[c_idx])
      cluster_means[c_idx, g] <- mean(in_state_2[cluster_subjs], na.rm = TRUE)
    }
  }
  p_tcm <- colMeans(cluster_means, na.rm = TRUE)

  data.frame(
    time   = t_grid,
    P_acm  = p_acm,
    P_tcm  = p_tcm
  )
}


#' Helper: percentile of follow-up time distribution.
#' The paper evaluates at tau_{0.4} and tau_{0.6}. We need to compute
#' these percentiles from the simulator's output.
follow_up_percentiles <- function(probs = c(0.4, 0.6),
                                  n_clusters = 1000,
                                  cluster_size_range = c(5, 15),
                                  right_censor = TRUE,
                                  left_truncate = FALSE,
                                  seed = 20260508) {
  d <- simulate_clusters(
    n_clusters = n_clusters,
    cluster_size_range = cluster_size_range,
    right_censor = right_censor,
    left_truncate = left_truncate,
    seed = seed
  )
  ## Follow-up time per subject: max(Tstop) - min(Tstart)
  fu_times <- aggregate(cbind(Tstart, Tstop) ~ id, data = d,
                        FUN = function(x) c(min = min(x), max = max(x)))
  ## Use Tstop max as the follow-up endpoint
  end_per_id <- tapply(d$Tstop, d$id, max)
  stats::quantile(end_per_id, probs = probs)
}
