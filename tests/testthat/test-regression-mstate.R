# Regression test: clusteredMSM's prodint_AJ() must agree with
# mstate::probtrans() on a progressive (acyclic) illness-death model
# to within numerical precision. mstate is in Suggests only, so this
# whole file is skipped when mstate is unavailable.

skip_if_not_installed("mstate")

# ---- Build a small progressive illness-death dataset --------------
# States: 1 = Healthy, 2 = Ill, 3 = Dead. 3 is absorbing.
# Allowed transitions: 1 -> 2, 1 -> 3, 2 -> 3 (no recovery).
# We simulate latent illness, healthy-death, ill-death and censoring
# times once, then materialise BOTH the mstate wide format and the
# clusteredMSM interval format from the *same* realised paths.
make_progressive_data <- function(n = 50, seed = 20260508) {
  set.seed(seed)

  t_ill   <- stats::rexp(n, rate = 0.6)   # latent H -> I
  t_dth_h <- stats::rexp(n, rate = 0.3)   # latent H -> D
  t_dth_i <- stats::rexp(n, rate = 0.5)   # latent I -> D (sojourn in I)
  t_cens  <- stats::runif(n, 1, 5)

  wide <- data.frame(
    id          = seq_len(n),
    illt        = NA_real_,
    illstatus   = 0L,
    deatht      = NA_real_,
    deathstatus = 0L
  )

  intervals <- vector("list", n)

  for (i in seq_len(n)) {
    if (t_ill[i] < t_dth_h[i] && t_ill[i] < t_cens[i]) {
      # H -> I, then either I -> D or censored ill
      ill_t <- t_ill[i]
      wide$illt[i]      <- ill_t
      wide$illstatus[i] <- 1L
      dth_t <- ill_t + t_dth_i[i]
      if (dth_t < t_cens[i]) {
        wide$deatht[i]      <- dth_t
        wide$deathstatus[i] <- 1L
        intervals[[i]] <- data.frame(
          id     = i,
          Tstart = c(0, ill_t),
          Tstop  = c(ill_t, dth_t),
          Sstart = c(1L, 2L),
          Sstop  = c(2L, 3L)
        )
      } else {
        wide$deatht[i]      <- t_cens[i]
        wide$deathstatus[i] <- 0L
        intervals[[i]] <- data.frame(
          id     = i,
          Tstart = c(0, ill_t),
          Tstop  = c(ill_t, t_cens[i]),
          Sstart = c(1L, 2L),
          Sstop  = c(2L, 2L)
        )
      }
    } else if (t_dth_h[i] < t_ill[i] && t_dth_h[i] < t_cens[i]) {
      # H -> D directly
      dth_t <- t_dth_h[i]
      wide$illt[i]        <- dth_t
      wide$illstatus[i]   <- 0L
      wide$deatht[i]      <- dth_t
      wide$deathstatus[i] <- 1L
      intervals[[i]] <- data.frame(
        id = i, Tstart = 0, Tstop = dth_t, Sstart = 1L, Sstop = 3L
      )
    } else {
      # Censored healthy
      wide$illt[i]        <- t_cens[i]
      wide$illstatus[i]   <- 0L
      wide$deatht[i]      <- t_cens[i]
      wide$deathstatus[i] <- 0L
      intervals[[i]] <- data.frame(
        id = i, Tstart = 0, Tstop = t_cens[i], Sstart = 1L, Sstop = 1L
      )
    }
  }

  list(wide = wide, intervals = do.call(rbind, intervals))
}


test_that("prodint_AJ agrees with mstate::probtrans on a progressive model", {
  dat <- make_progressive_data(n = 50)

  # ---- mstate pipeline -------------------------------------------
  tmat_m <- mstate::transMat(list(c(2, 3), 3, integer(0)),
                             names = c("Healthy", "Ill", "Dead"))

  ms_long <- mstate::msprep(
    time   = c(NA, "illt",      "deatht"),
    status = c(NA, "illstatus", "deathstatus"),
    data   = dat$wide,
    trans  = tmat_m,
    id     = "id"
  )

  cox_m <- survival::coxph(
    survival::Surv(Tstart, Tstop, status) ~ survival::strata(trans),
    data = ms_long, method = "breslow"
  )
  msf <- mstate::msfit(cox_m, trans = tmat_m, variance = FALSE)
  pt  <- mstate::probtrans(msf, predt = 0, variance = FALSE)
  # pt[[1]] is P(0, t) starting from state 1 (Healthy)
  ms_curve <- pt[[1]]   # columns: time, pstate1, pstate2, pstate3

  # ---- clusteredMSM pipeline --------------------------------------
  tmat_c <- trans_mat(list(c(2, 3), 3, integer(0)),
                      names = c("Healthy", "Ill", "Dead"))

  validate_intervals(dat$intervals, tmat_c)
  long_c <- intervals_to_long(dat$intervals, tmat_c)
  haz_c  <- fit_chaz(long_c, tmat_c)
  cm_curve <- prodint_AJ(haz_c, tmat_c, predt = 0, h = 1)

  # ---- Align both onto a common grid -----------------------------
  # Use a dense union of jump times so the comparison covers the
  # full support, not just one of the two outputs.
  grid <- sort(unique(c(ms_curve$time, cm_curve$time)))

  diffs_per_state <- vapply(1:3, function(j) {
    col <- paste0("pstate", j)
    # Initial value at t = 0 for column j is 1 if j == 1 (we start in
    # state 1) else 0 -- needed for grid points before the first jump.
    init <- if (j == 1L) 1 else 0
    ms_g <- step_interp(ms_curve$time, ms_curve[[col]], grid, init = init)
    cm_g <- step_interp(cm_curve$time, cm_curve[[col]], grid, init = init)
    max(abs(ms_g - cm_g))
  }, numeric(1))

  max_diff <- max(diffs_per_state)
  if (max_diff > 1e-10) {
    # Locate the worst time point across all three states for diagnostics.
    worst <- vapply(1:3, function(j) {
      col <- paste0("pstate", j)
      init <- if (j == 1L) 1 else 0
      ms_g <- step_interp(ms_curve$time, ms_curve[[col]], grid, init = init)
      cm_g <- step_interp(cm_curve$time, cm_curve[[col]], grid, init = init)
      d <- abs(ms_g - cm_g)
      c(state = j, t = grid[which.max(d)], diff = max(d))
    }, numeric(3))
    message("Per-state max abs differences:\n",
            paste(capture.output(print(t(worst))), collapse = "\n"))
  }

  expect_lt(max_diff, 1e-10)
})
