# cran-comments.md

## Test environments

* local: macOS Sonoma 14.x, R 4.4.x
* win-builder (devtools::check_win_devel()): R-devel
* R-hub (rhub::rhub_check()): Linux (Ubuntu 22.04, R-release and R-devel)
* GitHub Actions: ubuntu-latest, macos-latest, windows-latest

## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE is from `devtools::check_win_devel()` and has two parts, both
expected:

* "New submission" -- this is the first CRAN submission of the package.
* "Possibly misspelled words in DESCRIPTION": Aalen, Bakoyannis,
  Bandyopadhyay, Johansen, Multistate, multistate. The first four are
  author surnames (Aalen-Johansen estimator; Bakoyannis 2021;
  Bakoyannis & Bandyopadhyay 2022); "multistate" is the standard
  technical term in the survival-analysis literature. All spellings
  are intentional.

## This is a new release

This is the first submission of clusteredMSM.

The package implements methods from Bakoyannis (2021)
<doi:10.1111/biom.13327> for nonparametric analysis of clustered
multistate process data. It is intended as the canonical software
implementation accompanying that paper, and extends it to support
non-monotone (recovery) multistate processes.

## Inference calibration (simulation)

The cluster-bootstrap inference machinery was validated by Monte Carlo
simulation following the data-generating mechanism of Bakoyannis (2021)
Section 3: 500 replications, n = 40 clusters, cluster size U{5, 15},
progressive illness-death process (no recovery), right-censored,
B = 500 cluster
bootstrap replicates per replication. The pointwise estimand is the
all-cause-mortality transition probability evaluated at tau_0.4, the
40th percentile of the follow-up-time distribution under the
data-generating settings.

| Check                                  | Target | Empirical |
|----------------------------------------|--------|-----------|
| Pointwise CI coverage at tau_0.4       | 0.95   | 0.934     |
| Simultaneous band coverage             | 0.95   | 0.950     |
| K-S Type I error (shared clusters)     | 0.05   | 0.054     |
| K-S Type I error (cluster-randomized)  | 0.05   | 0.050     |
| K-S Type I error (independent obs.)    | 0.05   | 0.044     |

At tau_0.4 the empirical bias was -0.0050, and the average bootstrap
SE (0.0218) closely matched the Monte Carlo standard deviation
(0.0215), indicating that the point estimator is approximately
unbiased and that the bootstrap variance estimator accurately
captures the empirical sampling variability. All five coverage/error
rates are within Monte Carlo sampling variability of their nominal
levels.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
