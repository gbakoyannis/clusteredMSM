---
editor_options: 
  markdown: 
    wrap: 72
---

# clusteredMSM 0.1.0

## Initial CRAN release

-   Initial CRAN release of `clusteredMSM`.

## Highlights

-   Nonparametric estimation of population-averaged transition
    probabilities for clustered or independent multistate process data,
    with cluster-bootstrap pointwise confidence intervals and
    simultaneous confidence bands.
-   Two-sample Kolmogorov-Smirnov-type test, dispatched automatically
    when a grouping variable appears on the formula's right-hand side.
    Three regimes via a `design` argument: `"shared"` (multicenter,
    every cluster carries both groups), `"cluster_random"`
    (cluster-randomized trial, stratified cluster bootstrap), and
    `"indep_random"` (independent observational comparison, unstratified
    cluster bootstrap). `"auto"` infers the regime from the data.
-   Native support for non-monotone (recovery) multistate processes via
    a generic product-integral estimator (`prodint_AJ()`).
-   Self-contained: depends only on `survival` (which ships with R).
    Removes the dependency on `mstate` that limited earlier
    implementations to progressive models.
-   Estimation methodology from Bakoyannis (2021)
    <doi:10.1111/biom.13327>; two-sample inference for the
    cluster-randomized and independent-samples designs follows
    Bakoyannis and Bandyopadhyay (2022)
    <doi:10.1007/s10463-021-00819-x>.

## Public API

The package exposes a single user-facing function, `patp()`, with a
formula-based interface modelled after `survival::Surv()`:

``` r
# One-sample
patp(msm(Tstart, Tstop, Sstart, Sstop) ~ 1,
     data = mydata, tmat = tmat,
     id = "subj_id", cluster = "site",
     h = 1, j = 2, B = 1000, cband = TRUE)

# Two-sample (estimate + test)
patp(msm(Tstart, Tstop, Sstart, Sstop) ~ treatment,
     data = mydata, tmat = tmat,
     id = "subj_id", cluster = "site",
     h = 1, j = 2, B = 1000)
```

`patp()` returns an S3 object of class `patp` with `print()` and
`summary()` methods.

## Input format

Data are supplied in **interval format**: one row per mutually-exclusive
time interval per subject, with columns for interval start time, end
time, starting state, and ending state. Censoring is encoded as
`Sstart == Sstop` on the final row. Within each subject, intervals must
be temporally and spatially contiguous; this is enforced by strict
validation (`validate_intervals()`).

## Known limitations

-   The two-sample Kolmogorov-Smirnov-type test uses **unit weight** at
    every time point. Bakoyannis (2021) Section 2.5 recommends the
    harmonic-mean weight `W(t) = prod_p Y_p(t) / sum_p Y_p(t)`, which
    downweights regions where one group's at-risk set is small and tames
    tail instability of the difference estimator. A weighted variant
    following Bakoyannis (2021) Section 2.5 and Bakoyannis &
    Bandyopadhyay (2022) is planned for v0.2, alongside the linear and
    L2 tests.
-   Confidence bands use a standardized supremum construction on the
    cloglog scale; the q-weighted construction from Bakoyannis (2021)
    Section 2.3 will be added in v0.2.
