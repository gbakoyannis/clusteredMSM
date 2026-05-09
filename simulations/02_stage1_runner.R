#' Stage 1 sanity check: pre-CRAN simulation validation.
#'
#' Five calibration checks, modeled on Bakoyannis (2021) Section 3 but
#' scoped down for a single-evening run:
#'
#'   (A)    Pointwise CI coverage for P_{0,2}(tau_{0.4}) -- target ~0.95
#'   (B)    Simultaneous band coverage for P_{0,2}(.) -- target ~0.95
#'   (C.i)  K-S Type I error, design = "shared" -- target ~0.05
#'          (data: cluster_randomized = FALSE; every cluster carries both arms)
#'   (C.iia) K-S Type I error, design = "cluster_random" -- target ~0.05
#'          (data: cluster_randomized = TRUE; stratified-by-arm bootstrap)
#'   (C.iib) K-S Type I error, design = "indep_random" -- target ~0.05
#'          (same data as C.iia; unstratified bootstrap)
#'
#' Scope: one cell only.
#'   n_clusters = 40
#'   cluster size ~ U{5, 15}
#'   right censoring only (no left truncation)
#'   500 Monte Carlo replications
#'   500 bootstrap reps per analysis
#'
#' Expected runtime: 60-90 minutes on a laptop (three two-sample
#' bootstraps instead of one).
#'
#' Source the simulator and estimands first:
#'   source("00_simulator.R")
#'   source("01_estimands.R")
#' Then run this script.

suppressPackageStartupMessages({
  library(clusteredMSM)
})

## ---- Configuration ------------------------------------------------------

CONFIG <- list(
  n_clusters         = 40,
  cluster_size_range = c(5, 15),
  n_reps             = 500,
  B_boot             = 500,
  seed               = 20260508,
  results_dir        = "results"
)

if (!dir.exists(CONFIG$results_dir)) {
  dir.create(CONFIG$results_dir, recursive = TRUE)
}

tmat <- trans_mat(list(c(2, 3), 3, integer(0)),
                  names = c("Healthy", "Ill", "Dead"))


## ---- Compute true values for the estimands ------------------------------

cat("Step 0: Computing true values via large-sample Monte Carlo...\n")
cat("  (This is a one-time cost; takes ~1-2 minutes.)\n")

## Determine evaluation points
t_grid_truth <- seq(0.05, 2.5, by = 0.05)
truth_table <- compute_truth(
  t_grid             = t_grid_truth,
  n_clusters_truth   = 3000,
  cluster_size_range = CONFIG$cluster_size_range,
  seed               = CONFIG$seed
)

## Determine tau_{0.4} (the 40th percentile of follow-up time)
fu_quants <- follow_up_percentiles(
  probs              = c(0.4, 0.6),
  n_clusters         = 1000,
  cluster_size_range = CONFIG$cluster_size_range,
  right_censor       = TRUE,
  left_truncate      = FALSE,
  seed               = CONFIG$seed
)
tau_04 <- as.numeric(fu_quants["40%"])
cat(sprintf("  tau_0.4 = %.3f\n", tau_04))

## Find the closest grid point to tau_04
target_idx <- which.min(abs(truth_table$time - tau_04))
target_t   <- truth_table$time[target_idx]
true_P_acm <- truth_table$P_acm[target_idx]
cat(sprintf("  P_{0,2}(%.3f) (ACM) = %.4f (true value)\n", target_t, true_P_acm))


## ---- Helper: extract pointwise estimate, SE, and CI at a target time ----

approx_ci_at <- function(curves, t_target) {
  ## Step-function interpolation of the curve, SE, and CI to t_target
  ord <- order(curves$time)
  curves <- curves[ord, , drop = FALSE]
  idx <- max(which(curves$time <= t_target), 1L)
  list(
    P  = curves$P[idx],
    se = curves$se[idx],
    ll = curves$ll[idx],
    ul = curves$ul[idx]
  )
}


## ---- (A) Pointwise CI coverage at tau_0.4 -------------------------------

cat("\nStep 1: Pointwise CI coverage check (target: 0.95)\n")
cat(sprintf("  Reps: %d, B: %d, n_clusters: %d\n",
            CONFIG$n_reps, CONFIG$B_boot, CONFIG$n_clusters))
cat("  Running...\n")

set.seed(CONFIG$seed + 1)
results_A <- vector("list", CONFIG$n_reps)
pb <- txtProgressBar(min = 0, max = CONFIG$n_reps, style = 3)

for (r in seq_len(CONFIG$n_reps)) {
  d <- simulate_clusters(
    n_clusters         = CONFIG$n_clusters,
    cluster_size_range = CONFIG$cluster_size_range,
    right_censor       = TRUE,
    left_truncate      = FALSE
  )

  fit <- tryCatch(
    patp(msm(Tstart, Tstop, Sstart, Sstop) ~ 1,
         data = d, tmat = tmat,
         id = "id", cluster = "cluster",
         h = 1, j = 2, s = 0,
         B = CONFIG$B_boot, cband = FALSE,
         seed = CONFIG$seed + r),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    results_A[[r]] <- list(point = NA, se = NA, ll = NA, ul = NA, covered = NA)
  } else {
    ci <- approx_ci_at(fit$curves, target_t)
    covered <- !is.na(ci$ll) && !is.na(ci$ul) &&
               ci$ll <= true_P_acm && ci$ul >= true_P_acm
    results_A[[r]] <- list(point = ci$P, se = ci$se,
                           ll = ci$ll, ul = ci$ul,
                           covered = covered)
  }
  setTxtProgressBar(pb, r)
}
close(pb)

results_A_df <- do.call(rbind, lapply(results_A, as.data.frame))
coverage_A <- mean(results_A_df$covered, na.rm = TRUE)
bias_A     <- mean(results_A_df$point, na.rm = TRUE) - true_P_acm
mcsd_A     <- stats::sd(results_A_df$point, na.rm = TRUE)
ase_A      <- mean(results_A_df$se, na.rm = TRUE)

cat(sprintf("\n  Empirical coverage: %.3f (target 0.95, 95%% CI = %.3f - %.3f)\n",
            coverage_A,
            coverage_A - 1.96 * sqrt(coverage_A * (1 - coverage_A) / CONFIG$n_reps),
            coverage_A + 1.96 * sqrt(coverage_A * (1 - coverage_A) / CONFIG$n_reps)))
cat(sprintf("  Bias              : %+.4f\n", bias_A))
cat(sprintf("  MCSD              : %.4f\n", mcsd_A))
cat(sprintf("  ASE (avg SE)      : %.4f\n", ase_A))


## ---- (B) Simultaneous band coverage -------------------------------------

cat("\nStep 2: Simultaneous band coverage check (target: 0.95)\n")
cat("  Running...\n")

set.seed(CONFIG$seed + 2)
results_B <- logical(CONFIG$n_reps)
pb <- txtProgressBar(min = 0, max = CONFIG$n_reps, style = 3)

for (r in seq_len(CONFIG$n_reps)) {
  d <- simulate_clusters(
    n_clusters         = CONFIG$n_clusters,
    cluster_size_range = CONFIG$cluster_size_range,
    right_censor       = TRUE,
    left_truncate      = FALSE
  )

  fit <- tryCatch(
    suppressWarnings(
      patp(msm(Tstart, Tstop, Sstart, Sstop) ~ 1,
           data = d, tmat = tmat,
           id = "id", cluster = "cluster",
           h = 1, j = 2, s = 0,
           B = CONFIG$B_boot, cband = TRUE,
           seed = CONFIG$seed + r)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || !"ll.band" %in% names(fit$curves)) {
    results_B[r] <- NA
    setTxtProgressBar(pb, r)
    next
  }

  ## Check whether the band covers the truth at every grid point where
  ## both the band and the truth are defined.
  curves <- fit$curves
  in_range <- !is.na(curves$ll.band) & !is.na(curves$ul.band)
  if (sum(in_range) == 0L) {
    results_B[r] <- NA
    setTxtProgressBar(pb, r)
    next
  }

  ## Interpolate truth onto curves$time
  truth_at_grid <- approx(truth_table$time, truth_table$P_acm,
                          xout = curves$time, method = "linear",
                          rule = 2)$y

  ## Band covers iff at every in-range time, ll.band <= truth <= ul.band
  diffs_low  <- truth_at_grid[in_range] - curves$ll.band[in_range]
  diffs_high <- curves$ul.band[in_range] - truth_at_grid[in_range]
  covered_all <- all(diffs_low >= -1e-10, na.rm = TRUE) &&
                 all(diffs_high >= -1e-10, na.rm = TRUE)
  results_B[r] <- covered_all
  setTxtProgressBar(pb, r)
}
close(pb)

coverage_B <- mean(results_B, na.rm = TRUE)
cat(sprintf("\n  Empirical band coverage: %.3f (target 0.95, 95%% CI = %.3f - %.3f)\n",
            coverage_B,
            coverage_B - 1.96 * sqrt(coverage_B * (1 - coverage_B) / sum(!is.na(results_B))),
            coverage_B + 1.96 * sqrt(coverage_B * (1 - coverage_B) / sum(!is.na(results_B)))))


## ---- (C.i) Two-sample K-S Type I error: design = "shared" --------------

cat("\nStep 3a: Two-sample K-S Type I error, design = 'shared' (target: 0.05)\n")
cat("  Data: cluster_randomized = FALSE (every cluster carries both arms)\n")
cat("  Running...\n")

set.seed(CONFIG$seed + 3)
results_Ci <- logical(CONFIG$n_reps)
pb <- txtProgressBar(min = 0, max = CONFIG$n_reps, style = 3)

for (r in seq_len(CONFIG$n_reps)) {
  ## Generate two-sample data UNDER THE NULL (under_alternative = FALSE)
  d <- simulate_clusters(
    n_clusters         = CONFIG$n_clusters,
    cluster_size_range = CONFIG$cluster_size_range,
    two_sample         = TRUE,
    cluster_randomized = FALSE,
    under_alternative  = FALSE,
    right_censor       = TRUE,
    left_truncate      = FALSE
  )

  fit <- tryCatch(
    suppressWarnings(
      patp(msm(Tstart, Tstop, Sstart, Sstop) ~ arm,
           data = d, tmat = tmat,
           id = "id", cluster = "cluster",
           h = 1, j = 2, s = 0,
           design = "shared",
           B = CONFIG$B_boot,
           seed = CONFIG$seed + r)
    ),
    error = function(e) NULL
  )

  results_Ci[r] <- if (is.null(fit) || is.null(fit$test)) NA else
    fit$test$p.value < 0.05
  setTxtProgressBar(pb, r)
}
close(pb)

type1_error_Ci <- mean(results_Ci, na.rm = TRUE)
cat(sprintf("\n  Empirical Type I error: %.3f (target 0.05, 95%% CI = %.3f - %.3f)\n",
            type1_error_Ci,
            type1_error_Ci - 1.96 * sqrt(type1_error_Ci * (1 - type1_error_Ci) / sum(!is.na(results_Ci))),
            type1_error_Ci + 1.96 * sqrt(type1_error_Ci * (1 - type1_error_Ci) / sum(!is.na(results_Ci)))))


## ---- (C.iia + C.iib) Cluster-randomized data, two analyses --------------
## Generate one cluster-randomized dataset per replicate and run both
## design = "cluster_random" (stratified bootstrap) and design =
## "indep_random" (unstratified bootstrap) on it. This isolates the effect
## of the bootstrap scheme from data-generation noise.

cat("\nStep 3b: Two-sample K-S Type I error on cluster-randomized data\n")
cat("        C.iia: design = 'cluster_random' (stratified bootstrap)\n")
cat("        C.iib: design = 'indep_random'   (unstratified bootstrap)\n")
cat("  Data: cluster_randomized = TRUE (each cluster carries one arm)\n")
cat("  Running...\n")

set.seed(CONFIG$seed + 4)
results_Ciia <- logical(CONFIG$n_reps)
results_Ciib <- logical(CONFIG$n_reps)
pb <- txtProgressBar(min = 0, max = CONFIG$n_reps, style = 3)

for (r in seq_len(CONFIG$n_reps)) {
  d <- simulate_clusters(
    n_clusters         = CONFIG$n_clusters,
    cluster_size_range = CONFIG$cluster_size_range,
    two_sample         = TRUE,
    cluster_randomized = TRUE,
    under_alternative  = FALSE,
    right_censor       = TRUE,
    left_truncate      = FALSE
  )

  fit_iia <- tryCatch(
    suppressWarnings(
      patp(msm(Tstart, Tstop, Sstart, Sstop) ~ arm,
           data = d, tmat = tmat,
           id = "id", cluster = "cluster",
           h = 1, j = 2, s = 0,
           design = "cluster_random",
           B = CONFIG$B_boot,
           seed = CONFIG$seed + r)
    ),
    error = function(e) NULL
  )

  fit_iib <- tryCatch(
    suppressWarnings(
      patp(msm(Tstart, Tstop, Sstart, Sstop) ~ arm,
           data = d, tmat = tmat,
           id = "id", cluster = "cluster",
           h = 1, j = 2, s = 0,
           design = "indep_random",
           B = CONFIG$B_boot,
           seed = CONFIG$seed + r + 100000L)
    ),
    error = function(e) NULL
  )

  results_Ciia[r] <- if (is.null(fit_iia) || is.null(fit_iia$test)) NA else
    fit_iia$test$p.value < 0.05
  results_Ciib[r] <- if (is.null(fit_iib) || is.null(fit_iib$test)) NA else
    fit_iib$test$p.value < 0.05
  setTxtProgressBar(pb, r)
}
close(pb)

type1_error_Ciia <- mean(results_Ciia, na.rm = TRUE)
type1_error_Ciib <- mean(results_Ciib, na.rm = TRUE)

cat(sprintf("\n  C.iia (cluster_random): Type I error: %.3f (95%% CI = %.3f - %.3f)\n",
            type1_error_Ciia,
            type1_error_Ciia - 1.96 * sqrt(type1_error_Ciia * (1 - type1_error_Ciia) / sum(!is.na(results_Ciia))),
            type1_error_Ciia + 1.96 * sqrt(type1_error_Ciia * (1 - type1_error_Ciia) / sum(!is.na(results_Ciia)))))
cat(sprintf("  C.iib (indep_random):   Type I error: %.3f (95%% CI = %.3f - %.3f)\n",
            type1_error_Ciib,
            type1_error_Ciib - 1.96 * sqrt(type1_error_Ciib * (1 - type1_error_Ciib) / sum(!is.na(results_Ciib))),
            type1_error_Ciib + 1.96 * sqrt(type1_error_Ciib * (1 - type1_error_Ciib) / sum(!is.na(results_Ciib)))))


## ---- Summary ------------------------------------------------------------

summary_table <- data.frame(
  Check                     = c("Pointwise CI coverage",
                                "Simultaneous band coverage",
                                "K-S Type I error (shared)",
                                "K-S Type I error (cluster_random)",
                                "K-S Type I error (indep_random)"),
  Target                    = c(0.95, 0.95, 0.05, 0.05, 0.05),
  Empirical                 = c(coverage_A, coverage_B,
                                type1_error_Ci,
                                type1_error_Ciia,
                                type1_error_Ciib),
  N_reps                    = c(sum(!is.na(results_A_df$covered)),
                                sum(!is.na(results_B)),
                                sum(!is.na(results_Ci)),
                                sum(!is.na(results_Ciia)),
                                sum(!is.na(results_Ciib))),
  Pass                      = c(
    abs(coverage_A - 0.95) < 0.03,
    abs(coverage_B - 0.95) < 0.04,
    abs(type1_error_Ci - 0.05) < 0.025,
    abs(type1_error_Ciia - 0.05) < 0.025,
    abs(type1_error_Ciib - 0.05) < 0.025
  )
)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Stage 1 sanity check summary\n")
cat(strrep("=", 70), "\n", sep = "")
print(summary_table, row.names = FALSE)
cat("\n")

write.csv(summary_table,
          file.path(CONFIG$results_dir, "stage1_summary.csv"),
          row.names = FALSE)
saveRDS(list(config = CONFIG,
             truth = truth_table,
             tau_04 = tau_04,
             true_P_acm = true_P_acm,
             coverage_pointwise         = coverage_A,
             coverage_band              = coverage_B,
             type1_error_shared         = type1_error_Ci,
             type1_error_cluster_random = type1_error_Ciia,
             type1_error_indep_random   = type1_error_Ciib,
             pointwise_results          = results_A_df,
             band_results               = results_B,
             test_results_shared         = results_Ci,
             test_results_cluster_random = results_Ciia,
             test_results_indep_random   = results_Ciib),
        file.path(CONFIG$results_dir, "stage1_full.rds"))

cat(sprintf("Saved summary to %s/stage1_summary.csv\n", CONFIG$results_dir))
cat(sprintf("Saved full results to %s/stage1_full.rds\n", CONFIG$results_dir))
