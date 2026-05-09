# Generate `example_msm`: a synthetic clustered illness-death-with-recovery
# dataset used in README, examples, and tests.
#
# States: 1 = Healthy, 2 = Ill, 3 = Dead (absorbing).
# Transitions: 1<->2 (sickness/recovery), 1->3, 2->3.
# 8 clusters x 5 subjects = 40 subjects, balanced across two treatment arms.
# Treatment 1 has higher sickness hazard than treatment 0.
# Cluster-level multiplicative frailty induces within-cluster correlation.
#
# Run:  Rscript data-raw/example_msm.R
# Produces:
#   data/example_msm.rda           (lazy-loaded via data())
#   inst/extdata/example_data.csv  (loaded via system.file())

set.seed(2026)

n_clusters  <- 8L
per_cluster <- 5L

cluster_frailty <- exp(stats::rnorm(n_clusters, sd = 0.3))

sim_subject <- function(id, cluster, treatment) {
  fr <- cluster_frailty[cluster]
  rate_HI <- (if (treatment == 1L) 0.6 else 0.3) * fr  # Healthy -> Ill
  rate_HD <- 0.05 * fr                                 # Healthy -> Dead
  rate_IH <- 0.4  * fr                                 # Ill -> Healthy (recovery)
  rate_ID <- 0.3  * fr                                 # Ill -> Dead

  cens_time <- stats::runif(1, 2, 6)

  rows  <- list()
  state <- 1L
  t     <- 0

  repeat {
    if (state == 1L) {
      t_a <- stats::rexp(1, rate_HI); t_b <- stats::rexp(1, rate_HD)
      next_state <- if (t_a < t_b) 2L else 3L
    } else {
      t_a <- stats::rexp(1, rate_IH); t_b <- stats::rexp(1, rate_ID)
      next_state <- if (t_a < t_b) 1L else 3L
    }
    t_event <- min(t_a, t_b)

    if (t + t_event >= cens_time) {
      rows[[length(rows) + 1L]] <- data.frame(
        id = id, cluster = cluster, treatment = treatment,
        Tstart = t, Tstop = cens_time,
        Sstart = state, Sstop = state
      )
      break
    }

    rows[[length(rows) + 1L]] <- data.frame(
      id = id, cluster = cluster, treatment = treatment,
      Tstart = t, Tstop = t + t_event,
      Sstart = state, Sstop = next_state
    )
    t     <- t + t_event
    state <- next_state
    if (state == 3L) break
  }

  do.call(rbind, rows)
}

rows <- list()
sid  <- 1L
for (cl in seq_len(n_clusters)) {
  for (i in seq_len(per_cluster)) {
    rows[[length(rows) + 1L]] <- sim_subject(
      id = sid, cluster = cl, treatment = (sid %% 2L)
    )
    sid <- sid + 1L
  }
}
example_msm <- do.call(rbind, rows)
rownames(example_msm) <- NULL

# Round times for readability without sacrificing strict contiguity.
example_msm$Tstart <- round(example_msm$Tstart, 4)
example_msm$Tstop  <- round(example_msm$Tstop,  4)

# Sanity-check: validate against the package contract.
devtools::load_all(quiet = TRUE)
tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)),
                  names = c("Healthy", "Ill", "Dead"))
validate_intervals(example_msm, tmat)

save(example_msm,
     file = "data/example_msm.rda",
     compress = "bzip2", version = 2)
utils::write.csv(example_msm,
                 file = "inst/extdata/example_data.csv",
                 row.names = FALSE)

message(sprintf(
  "Wrote data/example_msm.rda and inst/extdata/example_data.csv (%d rows, %d subjects, %d clusters)",
  nrow(example_msm),
  length(unique(example_msm$id)),
  length(unique(example_msm$cluster))
))
