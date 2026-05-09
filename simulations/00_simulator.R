#' Simulator for the Bakoyannis (2021) clustered illness-death scenario.
#'
#' Generates clustered non-Markov illness-death data with informative
#' cluster size, matching the data-generating mechanism of Section 3 of
#' Bakoyannis (2021), Biometrics 77(2):533-546, doi:10.1111/biom.13327.
#'
#' Process structure:
#'   States: 1 = Healthy, 2 = Ill, 3 = Dead (absorbing)
#'   All subjects start at state 1 at t = 0
#'   Allowed transitions: 1->2, 1->3, 2->3
#'
#' Conditional cumulative intensities (given frailty v_i and cluster size m_i):
#'   A_{0,12}(t; v_i) = [0.25 + 0.25 * I{m_i <= E(M_1)}] * v_i * t
#'   A_{0,13}(t; v_i) = 0.25 * v_i * t
#'   A_{0,23}(t; v_i) = 0.5  * v_i * t
#'
#' Frailty: v_i ~ Gamma(shape = 1, scale = 1) -- shared within cluster.
#' This induces non-Markov marginal dynamics and ICS via the
#' I{m_i <= E(M_1)} indicator.
#'
#' For two-sample data under H1, the 1->2 intensity becomes:
#'   A_{0,p12}(t; v_i) = [0.25 + 0.5 * I(p=2) + 0.25 * I{m_i <= E(M_1)}] * v_i * t
#' for treatment arm p in {1, 2}.
#'
#' Arm allocation:
#'   cluster_randomized = FALSE (default): 1:1 split within each cluster
#'     (case (i) in Bakoyannis 2021 -- every cluster carries both arms).
#'   cluster_randomized = TRUE: whole clusters assigned to arms in
#'     alternating order (cluster 1 -> arm 1, cluster 2 -> arm 2, ...);
#'     all subjects in a cluster share that cluster's arm. This produces
#'     case (ii.a) / (ii.b) data (each cluster carries exactly one arm).

simulate_clusters <- function(n_clusters,
                              cluster_size_range,         # c(low, high) for U{low, ..., high}
                              two_sample = FALSE,
                              cluster_randomized = FALSE,
                              under_alternative = FALSE,
                              right_censor = TRUE,
                              left_truncate = FALSE,
                              tau = 3,                    # max follow-up (cens upper bound)
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (length(cluster_size_range) != 2L) {
    stop("cluster_size_range must be c(low, high)")
  }
  m_low  <- cluster_size_range[1]
  m_high <- cluster_size_range[2]
  E_M    <- (m_low + m_high) / 2     # E(M_1) for the indicator I{m_i <= E(M_1)}

  ## --- Cluster-level draws ---
  cluster_sizes <- sample(m_low:m_high, n_clusters, replace = TRUE)
  frailties     <- stats::rgamma(n_clusters, shape = 1, scale = 1)
  small_cluster <- cluster_sizes <= E_M       # logical, length n_clusters

  ## --- Per-subject simulation ---
  subj_rows <- list()
  next_id   <- 1L

  for (i in seq_len(n_clusters)) {
    v_i      <- frailties[i]
    m_i      <- cluster_sizes[i]
    is_small <- small_cluster[i]

    if (two_sample) {
      if (cluster_randomized) {
        ## Whole-cluster assignment, alternating arms (case ii.a / ii.b).
        cluster_arm <- if (i %% 2L == 1L) 1L else 2L
        arms <- rep(cluster_arm, m_i)
      } else {
        ## 1:1 allocation within cluster (round up if m_i odd) -- case (i).
        n_arm2  <- floor(m_i / 2)
        n_arm1  <- m_i - n_arm2
        if (n_arm1 < 1L || n_arm2 < 1L) {
          # Force at least one in each arm
          n_arm1 <- max(1L, n_arm1)
          n_arm2 <- max(1L, n_arm2)
        }
        arms <- c(rep(1L, n_arm1), rep(2L, n_arm2))
      }
    } else {
      arms <- rep(NA_integer_, m_i)
    }

    for (m in seq_len(m_i)) {
      arm <- arms[m]

      ## --- Hazard slopes (rates) for this subject ---
      lambda_13 <- 0.25 * v_i
      lambda_23 <- 0.50 * v_i

      ## 1->2 rate depends on cluster size (ICS) and arm (under H1)
      base_12 <- 0.25 + 0.25 * as.integer(is_small)
      if (two_sample && under_alternative && arm == 2L) {
        base_12 <- base_12 + 0.50
      }
      lambda_12 <- base_12 * v_i

      ## --- Independent latent times for transitions out of state 1 ---
      T_12 <- stats::rexp(1, rate = lambda_12)
      T_13 <- stats::rexp(1, rate = lambda_13)

      ## --- Right censoring time ---
      C <- if (right_censor) stats::runif(1, 0, tau) else Inf

      ## --- Left truncation time ---
      ## For state-occupation analyses under truncation, the paper sets
      ## L = 0 with probability 2/3 (so ~33% truncated). For transition
      ## probability analyses (s = 0.5), all subjects have L ~ Beta(1, 2).
      ## Caller can override via left_truncate argument.
      L <- if (left_truncate) {
        if (stats::runif(1) < 1/3) stats::rbeta(1, 1, 2) else 0
      } else {
        0
      }

      ## --- Determine the trajectory ---
      rows_subj <- build_subject_intervals(
        T_12 = T_12, T_13 = T_13, lambda_23 = lambda_23,
        L = L, C = C, id = next_id, cluster = i, arm = arm
      )

      if (!is.null(rows_subj)) {
        subj_rows[[length(subj_rows) + 1L]] <- rows_subj
      }
      next_id <- next_id + 1L
    }
  }

  out <- do.call(rbind, subj_rows)
  rownames(out) <- NULL
  out
}


#' Build the interval-format rows for one subject given latent transition
#' times, censoring, and left truncation.
#'
#' Returns NULL if the subject is fully left-truncated (entered observation
#' after the latent absorbing event time would have occurred).
#'
#' @keywords internal
build_subject_intervals <- function(T_12, T_13, lambda_23, L, C,
                                    id, cluster, arm) {
  ## Determine which transition out of state 1 happens first.
  if (T_12 < T_13) {
    ## Transition 1 -> 2 occurs at T_12 (if not censored/truncated)
    T_first   <- T_12
    next_state <- 2L
    ## Then in state 2, transition 2 -> 3 at T_first + Exp(lambda_23)
    T_23_gap  <- stats::rexp(1, rate = lambda_23)
    T_absorb  <- T_first + T_23_gap
  } else {
    T_first   <- T_13
    next_state <- 3L
    T_absorb  <- T_first
  }

  ## --- Apply left truncation ---
  ## If the subject's entire trajectory ends before L, they're not observed.
  if (L >= T_absorb && next_state == 3L) {
    ## Subject reached absorbing state before entering observation.
    ## Treat as not observed.
    return(NULL)
  }
  if (L >= C) {
    ## Truncated past their right-censoring time.
    return(NULL)
  }

  ## Effective end of follow-up
  end_obs <- min(C, T_absorb)
  censored <- C < T_absorb

  ## --- Build interval rows ---
  rows <- list()

  if (next_state == 2L && T_first < end_obs) {
    ## Subject went 1 -> 2 within follow-up

    ## First sojourn: state 1 from L to T_first
    if (T_first > L) {
      rows[[length(rows) + 1L]] <- data.frame(
        id      = id,
        cluster = cluster,
        arm     = arm,
        Tstart  = L,
        Tstop   = T_first,
        Sstart  = 1L,
        Sstop   = 2L
      )
    }

    ## Second sojourn: state 2 from T_first to end_obs
    rows[[length(rows) + 1L]] <- data.frame(
      id      = id,
      cluster = cluster,
      arm     = arm,
      Tstart  = T_first,
      Tstop   = end_obs,
      Sstart  = 2L,
      Sstop   = if (censored) 2L else 3L
    )

  } else {
    ## Either (a) direct 1 -> 3 within follow-up,
    ##        (b) censored in state 1, or
    ##        (c) went 1 -> 2 but censored before T_first.
    ## In all cases, output a single interval starting from state 1.
    Sstop_value <- if (!censored && next_state == 3L && T_first <= end_obs) {
      3L
    } else {
      1L  # censored in state 1
    }

    ## Edge case: T_first <= L means the 1->2 transition was supposed to
    ## happen before observation began. Treat as truncated (NULL).
    if (T_first <= L && !censored) {
      return(NULL)
    }

    rows[[length(rows) + 1L]] <- data.frame(
      id      = id,
      cluster = cluster,
      arm     = arm,
      Tstart  = L,
      Tstop   = end_obs,
      Sstart  = 1L,
      Sstop   = Sstop_value
    )
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}
