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

    boot_mat <- cluster_boot(long, cid = resample_col, B = B,
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
                            group_name) {

  resample_col <- if (has_cluster) "cluster" else "id"

  if (weighted) long <- add_cluster_sizes(long, "cluster", "id")

  group_levels <- sort(unique(long$group))
  long_0 <- long[long$group == group_levels[1L], , drop = FALSE]
  long_1 <- long[long$group == group_levels[2L], , drop = FALSE]

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

  diff_boot   <- cluster_boot(long, cid = resample_col, B = B,
                              fn = fn_diff, seed = seed)
  curves_boot <- cluster_boot(long, cid = resample_col, B = B,
                              fn = fn_curves,
                              seed = if (is.null(seed)) NULL else seed + 1L)

  n <- length(unique(long[[resample_col]]))

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

  test <- ks_pvalue(diff_point, diff_boot, n = n)

  list(curves       = curves,
       test         = list(statistic = test$statistic,
                           p.value   = test$p.value,
                           type      = "Kolmogorov-Smirnov"),
       n_clusters   = if (has_cluster) length(unique(long$cluster)) else NA_integer_,
       groups       = group_levels,
       group_name   = group_name)
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
