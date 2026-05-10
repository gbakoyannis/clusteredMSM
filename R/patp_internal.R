## Internal orchestration functions for patp().
##
## Not exported. These are split out so the public patp() function
## stays a thin formula-parsing dispatcher.

## ---- One-sample ----------------------------------------------------------

.patp_onesample <- function(long, tmat, has_cluster,
                            h, j, s,
                            weighted, LMAJ,
                            B, cband, level, seed) {

  resample_col <- if (has_cluster) "cluster" else "id"

  if (weighted) long <- add_cluster_sizes(long, "cluster", "id")

  ## Point estimate
  point <- .patp_point(long, tmat, h = h, j = j, s = s,
                       weighted = weighted, lmaj = LMAJ)

  curves <- data.frame(time = point$time, P = point$pstate)

  if (B > 0L) {
    grid <- point$time
    init <- if (h == j) 1 else 0

    fn_boot <- function(boot_data, ...) {
      pb <- tryCatch(
        .patp_point(boot_data, tmat, h = h, j = j, s = s,
                    weighted = weighted, lmaj = LMAJ),
        error = function(e) NULL
      )
      if (is.null(pb)) return(rep(NA_real_, length(grid)))
      step_interp(pb$time, pb$pstate, grid, init = init)
    }

    boot_mat <- cluster_boot(long, cluster = resample_col, B = B,
                             fn = fn_boot, seed = seed)

    se <- apply(boot_mat, 1L, stats::sd, na.rm = TRUE)
    ci <- ci_cloglog(curves$P, se, level = level)

    curves$se <- se
    curves$ll <- ci$ll
    curves$ul <- ci$ul

    if (cband) {
      band <- confidence_band(curves$P, boot_mat,
                              times = grid, level = level)
      curves$ll.band <- band$ll.band
      curves$ul.band <- band$ul.band
    }
  }

  list(curves     = curves,
       test       = NULL,
       n_clusters = if (has_cluster) length(unique(long$cluster)) else NA_integer_,
       groups     = NULL)
}


## ---- Two-sample ---------------------------------------------------------

.patp_twosample <- function(long, tmat, has_cluster,
                            h, j, s,
                            weighted, LMAJ,
                            B, cband, level, seed,
                            group_name,
                            design = c("auto", "shared",
                                       "cluster_random", "indep_random")) {

  design <- match.arg(design)
  resample_col <- if (has_cluster) "cluster" else "id"

  if (weighted) long <- add_cluster_sizes(long, "cluster", "id")

  group_levels <- sort(unique(long$group))
  long_0 <- long[long$group == group_levels[1L], , drop = FALSE]
  long_1 <- long[long$group == group_levels[2L], , drop = FALSE]

  ## ---- Resolve / validate design against the cluster/group structure ----
  ##
  ## Three regimes:
  ##   shared         every cluster carries both groups (case i,
  ##                  Bakoyannis 2021).
  ##   cluster_random each cluster carries one group, n_1 and n_2 fixed
  ##                  by design (case ii.a, Bakoyannis & Bandyopadhyay
  ##                  2022). Stratified bootstrap.
  ##   indep_random   each cluster carries one group, n_1 and n_2 random
  ##                  (independent observational comparison; case ii.b).
  ##                  Unstratified bootstrap.
  ##
  ## The two-sample asymptotic regime applies for both ii.a and ii.b, so
  ## the scaling factor is sqrt(n_1 n_2 / (n_1 + n_2)) for either; only
  ## the bootstrap (stratified vs. unstratified) differs.
  ##
  ## "auto" never picks "cluster_random" -- that is a stronger claim
  ## about the data-generating process and must be opted into by the
  ## user.
  if (has_cluster) {
    cluster_groups <- tapply(long$group, long$cluster,
                             function(x) length(unique(x)))
    n_one  <- sum(cluster_groups == 1L)
    n_both <- sum(cluster_groups == 2L)
  } else {
    n_one  <- NA_integer_
    n_both <- NA_integer_
  }

  if (design == "auto") {
    if (!has_cluster) {
      design <- "indep_random"
    } else if (n_one == 0L) {
      design <- "shared"
    } else if (n_both == 0L) {
      design <- "indep_random"            # safer default than cluster_random
      warning("Each cluster contains observations from only one group. ",
              "Using design = 'indep_random' (independent observational ",
              "comparison, unstratified bootstrap). If this is a ",
              "cluster-randomized trial with fixed n_1 and n_2, set ",
              "design = 'cluster_random' explicitly for stratified ",
              "resampling.",
              call. = FALSE)
    } else {
      stop("Mixed cluster structure detected (some clusters carry both ",
           "groups, some only one). This regime is not supported in ",
           "the current version.")
    }
  } else if (design == "shared") {
    if (has_cluster && n_one > 0L) {
      stop("design = \"shared\" requires every cluster to contain ",
           "observations from both groups, but some clusters contain ",
           "only one. If this is a cluster-randomized trial, use ",
           "design = \"cluster_random\"; if it is an independent ",
           "observational comparison, use design = \"indep_random\".")
    }
  } else if (design == "cluster_random") {
    if (!has_cluster) {
      stop("design = \"cluster_random\" requires a cluster column. ",
           "Cluster-randomized inference is defined at the cluster ",
           "level; subject-level stratification would be a different ",
           "procedure.")
    }
    if (n_both > 0L) {
      stop("design = \"cluster_random\" requires each cluster to ",
           "contain observations from only one group, but some ",
           "clusters contain both. Use design = \"shared\" if this ",
           "is a multicenter trial.")
    }
  } else if (design == "indep_random") {
    if (has_cluster && n_both > 0L) {
      stop("design = \"indep_random\" requires each cluster to ",
           "contain observations from only one group, but some ",
           "clusters contain both. Use design = \"shared\" if this ",
           "is a multicenter trial.")
    }
  }

  p0 <- .patp_point(long_0, tmat, h = h, j = j, s = s,
                    weighted = weighted, lmaj = LMAJ)
  p1 <- .patp_point(long_1, tmat, h = h, j = j, s = s,
                    weighted = weighted, lmaj = LMAJ)

  grid <- sort(unique(c(p0$time, p1$time)))
  init <- if (h == j) 1 else 0
  p0_grid <- step_interp(p0$time, p0$pstate, grid, init = init)
  p1_grid <- step_interp(p1$time, p1$pstate, grid, init = init)
  diff_point <- p1_grid - p0_grid

  fn_diff <- function(boot_data, ...) {
    d0 <- boot_data[boot_data$group == group_levels[1L], , drop = FALSE]
    d1 <- boot_data[boot_data$group == group_levels[2L], , drop = FALSE]

    pb0 <- tryCatch(
      .patp_point(d0, tmat, h = h, j = j, s = s,
                  weighted = weighted, lmaj = LMAJ),
      error = function(e) NULL
    )
    pb1 <- tryCatch(
      .patp_point(d1, tmat, h = h, j = j, s = s,
                  weighted = weighted, lmaj = LMAJ),
      error = function(e) NULL
    )
    if (is.null(pb0) || is.null(pb1)) {
      return(rep(NA_real_, length(grid)))
    }

    g1 <- step_interp(pb1$time, pb1$pstate, grid, init = init)
    g0 <- step_interp(pb0$time, pb0$pstate, grid, init = init)
    g1 - g0
  }

  fn_curves <- function(boot_data, ...) {
    d0 <- boot_data[boot_data$group == group_levels[1L], , drop = FALSE]
    d1 <- boot_data[boot_data$group == group_levels[2L], , drop = FALSE]
    pb0 <- tryCatch(
      .patp_point(d0, tmat, h = h, j = j, s = s,
                  weighted = weighted, lmaj = LMAJ),
      error = function(e) NULL
    )
    pb1 <- tryCatch(
      .patp_point(d1, tmat, h = h, j = j, s = s,
                  weighted = weighted, lmaj = LMAJ),
      error = function(e) NULL
    )
    g0 <- if (is.null(pb0)) rep(NA_real_, length(grid)) else
            step_interp(pb0$time, pb0$pstate, grid, init = init)
    g1 <- if (is.null(pb1)) rep(NA_real_, length(grid)) else
            step_interp(pb1$time, pb1$pstate, grid, init = init)
    c(g0, g1)
  }

  ## ---- Bootstrap: stratified iff design == "cluster_random" ----
  if (design == "cluster_random") {
    cluster_to_group <- tapply(long$group, long$cluster,
                               function(x) as.character(x[1L]))
    strata_vec <- as.character(cluster_to_group)
    names(strata_vec) <- names(cluster_to_group)
  } else {
    strata_vec <- NULL
  }

  diff_boot   <- cluster_boot(long, cluster = resample_col, B = B,
                              fn = fn_diff, strata = strata_vec,
                              seed = seed)
  curves_boot <- cluster_boot(long, cluster = resample_col, B = B,
                              fn = fn_curves, strata = strata_vec,
                              seed = if (is.null(seed)) NULL else seed + 1L)

  n <- length(unique(long[[resample_col]]))

  ## ---- Asymptotic scaling factor ----
  ## sqrt(n) for case (i) [shared]; sqrt(n_1 n_2 / (n_1 + n_2)) for
  ## both two-independent-sample regimes (ii.a and ii.b). The factor
  ## is determined by the asymptotic regime, not by whether the
  ## bootstrap is stratified.
  if (design == "shared") {
    scale_factor <- sqrt(n)
  } else {
    n1 <- length(unique(long_0[[resample_col]]))
    n2 <- length(unique(long_1[[resample_col]]))
    scale_factor <- sqrt(n1 * n2 / (n1 + n2))
  }

  ## Per-group SEs and CIs
  ng <- length(grid)
  boot_0 <- curves_boot[seq_len(ng), , drop = FALSE]
  boot_1 <- curves_boot[(ng + 1L):(2L * ng), , drop = FALSE]
  se_0   <- apply(boot_0, 1L, stats::sd, na.rm = TRUE)
  se_1   <- apply(boot_1, 1L, stats::sd, na.rm = TRUE)
  ci_0   <- ci_cloglog(p0_grid, se_0, level = level)
  ci_1   <- ci_cloglog(p1_grid, se_1, level = level)

  curves <- rbind(
    data.frame(time = grid, P = p0_grid,
               se = se_0, ll = ci_0$ll, ul = ci_0$ul,
               group = group_levels[1L], stringsAsFactors = FALSE),
    data.frame(time = grid, P = p1_grid,
               se = se_1, ll = ci_1$ll, ul = ci_1$ul,
               group = group_levels[2L], stringsAsFactors = FALSE)
  )

  if (cband) {
    band_0 <- confidence_band(p0_grid, boot_0, times = grid, level = level)
    band_1 <- confidence_band(p1_grid, boot_1, times = grid, level = level)
    curves$ll.band <- c(band_0$ll.band, band_1$ll.band)
    curves$ul.band <- c(band_0$ul.band, band_1$ul.band)
  }

  test <- ks_pvalue(diff_point, diff_boot, scale = scale_factor)

  list(curves       = curves,
       test         = list(statistic = test$statistic,
                           p.value   = test$p.value,
                           type      = "Kolmogorov-Smirnov"),
       n_clusters   = if (has_cluster) length(unique(long$cluster)) else NA_integer_,
       groups       = group_levels,
       group_name   = group_name,
       design       = design)
}


## ---- The bare point estimator (formerly exported as patp_point) ----

.patp_point <- function(data, tmat, h, j, s = 0,
                        weighted = FALSE, lmaj = FALSE) {
  if (lmaj) return(.patp_lmaj(data, tmat, h, j, s, weighted))

  w <- if (weighted) {
    if (!"clust.size" %in% names(data)) {
      stop("weighted = TRUE requires a 'clust.size' column")
    }
    1 / data$clust.size
  } else NULL

  haz <- fit_chaz(data, tmat, weights = w)
  P   <- prodint_AJ(haz, tmat, predt = s, h = h)
  out <- P[, c("time", paste0("pstate", j))]
  names(out) <- c("time", "pstate")
  out
}


.patp_lmaj <- function(data, tmat, h, j, s, weighted = FALSE) {
  if (!"id" %in% names(data)) {
    stop("'data' must contain an 'id' column for landmark estimation")
  }
  in_h <- state_at(data, s, id = "id")
  in_h_ids <- in_h[["id"]][in_h$state %in% h]
  if (length(in_h_ids) == 0L) {
    stop("no subjects in state(s) ", paste(h, collapse = ","),
         " at time s = ", s)
  }
  data_lm <- cut_at_lm(data, s)
  data_lm <- data_lm[data_lm[["id"]] %in% in_h_ids, , drop = FALSE]

  w <- if (weighted) {
    if (!"clust.size" %in% names(data_lm)) {
      stop("weighted = TRUE requires a 'clust.size' column")
    }
    1 / data_lm$clust.size
  } else NULL

  haz <- fit_chaz(data_lm, tmat, weights = w)
  P   <- prodint_AJ(haz, tmat, predt = s, h = h)
  out <- P[, c("time", paste0("pstate", j))]
  names(out) <- c("time", "pstate")
  out
}
