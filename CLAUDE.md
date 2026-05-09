# CLAUDE.md — clusteredMSM

This file is read at the start of every Claude Code session. It encodes
the key decisions and conventions for this package so they don't have
to be re-derived. When in doubt, follow what's written here; when
something here looks wrong, surface it before changing direction
silently.

## What this package is

`clusteredMSM` is the canonical R package implementing Bakoyannis (2021),
*Nonparametric analysis of nonhomogeneous multistate processes with
clustered observations*, Biometrics 77(2):533–546, doi:10.1111/biom.13327.

It supersedes the `clustered-multistate` GitHub repo, which was a
collection of scripts. This package:

- Estimates population-averaged transition probabilities
  P(X(t) = j | X(s) = h) for clustered or independent multistate data.
- Provides cluster-bootstrap pointwise CIs and simultaneous bands.
- Conducts two-sample Kolmogorov-Smirnov-type tests for curve equality
  (linear and L2 tests planned for v0.2).
- **Supports non-monotone (recovery) multistate models**, which the
  original `mstate`-based implementation could not.

## Public API: a single function with a formula interface

Users interact with the package through one function: `patp()`. The
formula's right-hand side determines whether it's a one-sample or
two-sample analysis.

```r
# One-sample: estimate the curve with CIs
patp(msm(Tstart, Tstop, Sstart, Sstop) ~ 1,
     data = mydata, tmat = tmat,
     id = "subj_id", cluster = "site",
     h = 1, j = 2, B = 1000, cband = TRUE)

# Two-sample: estimate both curves AND test their equality
patp(msm(Tstart, Tstop, Sstart, Sstop) ~ treatment,
     data = mydata, tmat = tmat,
     id = "subj_id", cluster = "site",
     h = 1, j = 2, B = 1000)
```

Argument conventions:

- **`id`**: required; subject identifier column.
- **`cluster`**: optional (default `NA`). When `NA`, bootstrap resamples
  individuals; when supplied, it resamples whole clusters.
- **`B`**: bootstrap replications, default 1000. `B = 0` returns point
  estimate only (one-sample only; two-sample requires `B > 0`).
- **`weighted = TRUE`**: requires `cluster`; uses inverse-cluster-size
  weighting.

Return value: an S3 object of class `patp` with `print()` and
`summary()` methods. The object has slots `$curves`, `$test` (NULL for
one-sample), `$call`, `$formula`, plus metadata.

## Input data format

Users supply data with **one row per mutually-exclusive time interval
per subject**. Required columns (any names; the formula maps them):

- `Tstart`, `Tstop`: numeric times bracketing the interval.
- `Sstart`, `Sstop`: integer states at start and end of the interval.
- A subject ID column (named via `id` argument).
- Optionally a cluster ID column (named via `cluster` argument).
- Optionally a binary grouping variable (named via formula RHS).

Within each subject, rows must be:
- **Temporally contiguous:** `Tstop[k] == Tstart[k+1]`.
- **Spatially contiguous:** `Sstop[k] == Sstart[k+1]`.

Censoring is encoded as `Sstart == Sstop` and is permitted **only on
the final row** of a subject's record. Absorbing states are the natural
endpoint and have no row after them.

Validation is **strict** — any violation triggers an error with an
informative message. See `validate_intervals()` for the full rule list.

## Architecture: layered, single-responsibility

```
trans_mat ────────────────────────┐
                                  │
msm + parse_msm_formula ──────────┤
                                  │
validate_intervals + intervals_to_long ─┐
                                        │
fit_chaz ────────┐                      │
                 ├──> .patp_point ──────┤
prodint_AJ ──────┤                      │
                 │                      │
state_at,        │                      │
cut_at_lm ────> .patp_lmaj              │
                                        ├──> patp() (public API)
cluster_boot ───────────────────────────┤
                                        │
inference (ci_cloglog,                  │
           confidence_band,             │
           ks_pvalue) ──────────────────┘

methods (print.patp, summary.patp): user-facing display
utils (check_clusters, add_cluster_sizes, step_interp): shared helpers
```

| Layer | File(s) | Responsibility |
|-------|---------|----------------|
| Constructor | `msm.R` | Class-tagged matrix for the formula LHS. |
| Formula parsing | `parse_formula.R` | Internal `parse_msm_formula()`. |
| Validation | `validate_intervals.R` | Strict rules + transform to long format. |
| Hazards | `fit_chaz.R` | Cox stratified by transition + Breslow hazards. |
| Product integral | `prodint_AJ.R` | Aalen-Johansen, recovery-capable. |
| Landmark helpers | `state_at.R`, `cut_at_lm.R` | Replace `mstate::xsect`/`cutLMms`. |
| Estimator (internal) | `patp_internal.R` | `.patp_onesample`, `.patp_twosample`, `.patp_point`, `.patp_lmaj`. |
| Bootstrap engine | `cluster_boot.R` | Generic, statistic-agnostic. |
| Inference | `inference.R` | `ci_cloglog`, `confidence_band`, `ks_pvalue`. |
| Public API | `patp.R` | The single user-facing entry point. |
| S3 methods | `methods.R` | `print.patp`, `summary.patp`. |
| Utilities | `utils.R` | Shared helpers. |

**Why this matters:** any one component can be swapped without touching
the others. We replaced `mstate::probtrans()` with `prodint_AJ()` to
get recovery support; nothing else changed. The single public function
hides internal complexity from the user.

### Adding new tests (planned: linear, L2)

The bootstrap engine is generic. Adding a new two-sample test means:

1. Write a pure `<name>_pvalue(diff_point, diff_boot, ...)` in
   `R/inference.R` (or its own file).
2. Add a `test_type` argument to `patp()` and dispatch in
   `.patp_twosample` based on it.
3. Add unit tests that don't require running the full bootstrap.

Do NOT duplicate the bootstrap construction. The bootstrap matrix is
the expensive part; computing different statistics from the same matrix
is essentially free.

## What replaced what (mstate removal)

The package depends only on `survival` (which ships with R itself).
Everything previously borrowed from `mstate` has a native replacement:

| `mstate` function | Replacement | Notes |
|---|---|---|
| `transMat()` | `trans_mat()` | Same API, supports cycles. |
| `msprep()` | (removed) | Users now supply interval-format data directly. |
| `msfit()` | `fit_chaz()` | Wraps `survival::coxph` + `basehaz`. |
| `probtrans()` | `prodint_AJ()` | **Key change** — handles cyclic transition matrices. |
| `xsect()` | `state_at()` | Subject's state at landmark time. |
| `cutLMms()` | `cut_at_lm()` | Truncate intervals at landmark time. |
| `msboot()` | `cluster_boot()` | Generic; takes any statistic function. |

`mstate` is in `Suggests`, not `Imports`. Used only for optional
regression tests that compare `prodint_AJ` against `probtrans` on
progressive examples.

## Naming conventions

- **Package name:** `clusteredMSM` (camel-case, no underscore, CRAN-safe).
- **Functions:** `snake_case` (`patp`, `fit_chaz`, `cluster_boot`).
- **Internal helpers:** leading dot (`.patp_onesample`, `.patp_point`).
  Not exported; not in NAMESPACE.
- **The public entry point is `patp`.** This name has brand equity from
  the paper and the original GitHub code — do not change it.
- **The S3 class returned by `patp` is `c("patp", "list")`.**
- **The S3 class for validated long-format data is `c("cmsdata", "data.frame")`.**

## Argument naming convention

- `id`: subject ID column name. Required, no default.
- `cluster`: cluster ID column name. Default `NA`.
- (Note: in early drafts this was called `cid`; renamed to `cluster`
  for clarity.)
- `tmat`: transition matrix from `trans_mat()`.
- `h`, `j`: starting and ending states for the estimand.
- `s`: conditioning time, default 0.
- `B`: bootstrap replications, default 1000.

## Testing strategy

- **Framework:** `testthat` 3rd edition.
- **One test file per source file:** `R/foo.R` ↔ `tests/testthat/test-foo.R`.
- **Unit tests** for pure functions (`prodint_AJ`, `ci_cloglog`,
  `cluster_boot`, `validate_intervals`, etc.).
- **Integration tests** for `patp` use a synthetic clustered
  illness-death dataset (helper `make_interval_data()` in `test-patp.R`).
- **Critical regression test:** `tests/testthat/test-regression-mstate.R`
  generates a 50-subject progressive illness-death dataset, runs both
  the `mstate` pipeline (`transMat` → `msprep` → `coxph` → `msfit` →
  `probtrans`) and the `clusteredMSM` pipeline (`trans_mat` →
  `validate_intervals` → `intervals_to_long` → `fit_chaz` →
  `prodint_AJ`) on the same realised paths, aligns both onto the union
  jump-time grid via `step_interp()`, and asserts element-wise agreement
  to within `1e-10`. Verified: the two pipelines agree to machine
  epsilon (max abs diff ~2e-16) on this acyclic example, so the 1e-10
  threshold has nine orders of magnitude of headroom.
- **Do not require `mstate` in regular tests** (only `Suggests`). Wrap
  `mstate`-dependent tests in `skip_if_not_installed("mstate")`.

## Coding rules

- **Always namespace external function calls:** `survival::coxph`,
  `stats::sd`. CRAN requires this.
- **No `library()` calls inside `R/*.R` files.** Use `@importFrom` in
  roxygen or `survival::` prefixes.
- **No `T` / `F` — always `TRUE` / `FALSE`.**
- **No writes outside `tempdir()`** in examples or tests.
- **No `cat()` for messages** in functions; use `message()` or `warning()`.
- **Examples should run in <5 seconds** or be wrapped in `\donttest{}`.
- **No `set.seed()` in package functions** (only in examples and tests).
  Bootstrap functions take an explicit `seed` argument instead.
- **Input validation up front.** Public functions check inputs and stop
  with informative errors; internal helpers may skip validation if the
  caller has already done it.

## Validation contract (validate_intervals)

Strict rules, applied in order — stop on first violation:

1. Required columns present (`id`, `Tstart`, `Tstop`, `Sstart`, `Sstop`).
2. Numeric / integer types where expected.
3. `Tstart < Tstop` on every row.
4. States are integers in `1..K`.
5. Within each subject, rows sorted and temporally contiguous
   (`Tstop[k] == Tstart[k+1]`, within `tol`).
6. Within each subject, rows spatially contiguous
   (`Sstop[k] == Sstart[k+1]`).
7. No row has `Sstart` equal to an absorbing state.
8. No row follows one whose `Sstop` is absorbing.
9. Each non-censored transition (`Sstart != Sstop`) is allowed by `tmat`.
10. `Sstart == Sstop` (censoring) only on the final row of a subject.

## Bootstrap conventions

- **Cluster-level resampling** if `cluster` was supplied; otherwise
  individual-level resampling.
- **Re-ID resampled clusters as 1..n** so duplicate draws stay distinct
  for downstream cluster-aware operations.
- **Pre-split data by cluster ID** for fast subsetting in the loop.
- **Tolerate pathological replicates.** Errors in `fn` produce NA
  columns rather than crashing the whole bootstrap.
- **Align jump times to a common grid** via `step_interp()` before
  assembling bootstrap matrices.
- **Bootstrap once, on the original probability scale.** The cluster
  bootstrap produces replicates of \eqn{P^*(t)}, matching Bakoyannis
  (2021) Theorem 2 (\eqn{\sqrt{n}\{P^*(t) - \hat P(t)\}}). The
  probability-scale SE is `apply(boot_mat, 1, sd, na.rm = TRUE)` --
  do **not** divide by `sqrt(n)` (the bootstrap SD is already an
  estimate of \eqn{\mathrm{SE}(\hat P(t))}).
- **Cloglog CIs use the delta method**, not a separate cloglog-scale
  bootstrap. `ci_cloglog(point, se)` takes the probability-scale SE
  and applies \eqn{\mathrm{SE}_g = \mathrm{SE}(\hat P) / |\hat P
  \log \hat P|}; the simultaneous band uses the same delta-method
  \eqn{\mathrm{SE}_g} to studentize cloglog residuals computed from
  the existing replicate matrix. One bootstrap, one scale.
- **Two-sample test scaling depends on the asymptotic regime, not on
  the bootstrap scheme.** `ks_pvalue(diff_point, diff_boot, scale)`
  multiplies both the observed sup-statistic and the
  centered-bootstrap analogue by `scale`; the empirical p-value is
  invariant to `scale` but the reported statistic is on the correct
  asymptotic scale:
    - `design = "shared"` (case i, all clusters carry both groups,
      multicenter): `scale = sqrt(n)`, n = total clusters.
      Bakoyannis (2021) Theorem 3.
    - `design = "cluster_random"` or `"indep_random"` (case ii.a /
      ii.b, two-independent-samples regime):
      `scale = sqrt(n_1 * n_2 / (n_1 + n_2))`. Bakoyannis &
      Bandyopadhyay (2022) Theorem 2.
- **The bootstrap differs across the two-sample regimes:**
    - `"shared"`: unstratified cluster bootstrap.
    - `"cluster_random"`: cluster bootstrap stratified by group
      (`cluster_boot(..., strata = cluster_to_group)`); fixes per-
      group cluster counts at \eqn{n_1, n_2}. Justified by the
      randomization-by-design assumption.
    - `"indep_random"`: unstratified cluster bootstrap (\eqn{n_1, n_2}
      may vary across replicates). Same scaling as
      `"cluster_random"`, but the bootstrap is *not* stratified
      because \eqn{n_1, n_2} are random in the data-generating
      process.
- **Stratified resampling** is opt-in via `cluster_boot(..., strata
  = ...)`. `strata` is a named vector mapping cluster ID (as
  character) to stratum label. Per-stratum cluster counts are
  preserved in every replicate.

### Two-sample designs

`patp()` accepts a `design` argument with values `"auto"` (default),
`"shared"`, `"cluster_random"`, `"indep_random"`. The supported
regimes:

- **Case (i): `"shared"` -- dependent groups (multicenter).** Every
  cluster carries both groups. Unstratified bootstrap; `sqrt(n)`
  scaling. Bakoyannis (2021).
- **Case (ii.a): `"cluster_random"` -- cluster-randomized trial.**
  Each cluster carries one group, \eqn{n_1, n_2} fixed by the
  randomization. Stratified-by-group bootstrap;
  `sqrt(n_1 n_2 / (n_1 + n_2))` scaling. Bakoyannis &
  Bandyopadhyay (2022). Must be opted into by the user --
  `"auto"` will not select it.
- **Case (ii.b): `"indep_random"` -- independent observational
  comparison.** Each cluster carries one group, but \eqn{n_1, n_2}
  are random. Unstratified bootstrap; same scaling as ii.a.
- **Case (iii): mixed cluster structure** (some clusters carry both
  groups, some only one) is not supported in v0.1; errors in every
  regime. Planned for v0.2 alongside the linear and \eqn{L^2} tests.

`design = "auto"` infers the regime from the data: all-clusters-
both -> `"shared"` (silent, unambiguous); all-clusters-one ->
`"indep_random"` *with a warning* nudging the user to opt into
`"cluster_random"` if the per-group cluster counts were fixed by
randomization; mixed -> error. Validation happens inside
`.patp_twosample` *before* any bootstrap work, so users get a fast
and informative error.

## Performance notes

- Hot path during bootstrap is `prodint_AJ()`. O(J·K²) per replicate.
- For typical problems (K ≤ 5, J ≤ a few thousand, B = 1000) runtime
  is fine in pure R.
- Two-sample analysis bootstraps both curves and the difference in
  one engine call (saves ~half the runtime vs. naive 2 separate
  bootstraps).
- If profiling identifies `prodint_AJ()` as the bottleneck, the first
  optimization is pre-splitting `haz_df` by jump time. After that,
  Rcpp.

## CRAN submission checklist

Pre-submission:
- [ ] `DESCRIPTION` has real maintainer email and ORCID
- [ ] `R/zzz.R` includes `utils::globalVariables(c("Tstart", "Tstop", ...))`
- [ ] `inst/extdata/example_data.csv` ships with the package
- [x] Regression test against `mstate::probtrans()` passes (skipped
      gracefully if `mstate` not installed)
- [ ] `devtools::document()` regenerates `man/` and `NAMESPACE`
- [ ] `devtools::check()` returns 0 errors, 0 warnings, 0 notes locally
- [ ] `devtools::check_win_devel()` clean
- [ ] `devtools::check_mac_release()` clean
- [ ] Examples run in <5 seconds (or use `\donttest{}`)
- [ ] Vignette builds cleanly via `devtools::build_vignettes()`

Submission:
- [ ] `cran-comments.md` updated
- [ ] `devtools::submit_cran()` invoked
- [ ] Confirmation email confirmed

Post-submission:
- [ ] Address reviewer comments within a week
- [ ] Bump `Version:` and update `NEWS.md` for each resubmission

## What citation("clusteredMSM") returns

`inst/CITATION` pins the canonical citation to the Biometrics 2021
paper. Do not change this.

## Things NOT to do

- Do not reintroduce `mstate` as `Imports`. It's `Suggests` only.
- Do not edit `prodint_AJ()` to "match `probtrans()` for performance"
  — the whole point is that it generalizes beyond `probtrans`.
- Do not split `patp` back into `patp` + `patp_test` — the unified
  formula API is the public surface.
- Do not change the public function name `patp`.
- Do not let `R/utils.R` become a junk drawer.
- Do not duplicate the bootstrap construction across helpers.

## Open TODOs

- [x] Verify `prodint_AJ()` matches `mstate::probtrans()` to ~1e-12 on
      the example data (done: agrees to ~2e-16 on a 50-subject
      progressive illness-death case; see `test-regression-mstate.R`)
- [x] Ship example data: `data/example_msm.rda` (lazy-loaded) +
      `inst/extdata/example_data.csv` (CSV mirror); both generated
      reproducibly from `data-raw/example_msm.R`
      (`set.seed(2026)`, 40 subjects, 8 clusters,
      illness-death-with-recovery). Documented in `R/data.R`.
      `LazyData: true` added to `DESCRIPTION`.
- [ ] Add `R/zzz.R` with `globalVariables(...)` for CRAN's
      "no visible binding" notes (not currently needed — `R CMD
      check` is clean with 0 notes; revisit if notes appear)
- [ ] Implement linear and L2 tests (planned v0.2)
- [ ] Add `plot.patp()` S3 method
- [ ] Decide whether to drop `mstate` from `Suggests` after regression
      test is in place

## Useful commands

```r
devtools::document()                # regenerate man/ and NAMESPACE
devtools::test()                    # run all tests
devtools::check()                   # full R CMD check
devtools::build()                   # build source tarball
devtools::check_win_devel()         # CRAN-relevant Windows R-devel check
devtools::install()                 # install locally
covr::package_coverage()            # test coverage report
```

## References

- Bakoyannis, G. (2021). Nonparametric analysis of nonhomogeneous
  multistate processes with clustered observations. *Biometrics*,
  77(2), 533-546. doi:10.1111/biom.13327
- Aalen, O. O., & Johansen, S. (1978). An empirical transition matrix
  for non-homogeneous Markov chains based on censored observations.
  *Scandinavian Journal of Statistics*, 5, 141-150.
- Putter, H., Fiocco, M., & Geskus, R. B. (2007). Tutorial in
  biostatistics: competing risks and multi-state models.
  *Statistics in Medicine*, 26(11), 2389-2430.
