# Simulations

This directory contains simulation studies for the `clusteredMSM`
package, modelled on the simulation design of Bakoyannis (2021),
Biometrics 77(2):533-546, doi:10.1111/biom.13327.

## Files

- `00_simulator.R` -- Data-generating function for the clustered
  non-Markov illness-death scenario with informative cluster size.
- `01_estimands.R` -- Computes true population-averaged probabilities
  via large-sample Monte Carlo (no closed form due to the frailty
  + ICS structure).
- `02_stage1_runner.R` -- Stage 1 sanity-check: three calibration
  checks before CRAN submission. Single-cell scope.

## Running the Stage 1 sanity check

From an R session:

```r
setwd("simulations")
source("00_simulator.R")
source("01_estimands.R")
source("02_stage1_runner.R")
```

Or from the command line:

```bash
cd simulations
Rscript 02_stage1_runner.R
```

(Make sure `clusteredMSM` is installed first.)

## What Stage 1 checks

Three calibration properties of the inference machinery, each with
500 Monte Carlo replications and 500 bootstrap reps:

1. **Pointwise CI coverage** of `P_{0,2}(tau_{0.4})` -- target 0.95
2. **Simultaneous band coverage** of `P_{0,2}(.)` -- target 0.95
3. **Two-sample K-S Type I error** under the null -- target 0.05

Single cell: `n_clusters = 40`, cluster size `~ U{5, 15}`, right
censoring only, no left truncation.

## What "passing" means

A check passes if the empirical rate is within 2 binomial standard
errors of the nominal rate (the 95% binomial confidence interval
based on 500 reps). Specifically:

- Coverage 0.93-0.97: pass
- Type I error 0.03-0.07: pass

Tighter tolerances would require more replications (paper used 1000).

## What to do if a check fails

If empirical coverage is well below 0.95 (say, < 0.92), or empirical
Type I error is well above 0.05 (say, > 0.08), the calibration of
the inference machinery is off. Possible causes:

- Bootstrap B is too small for the band check -> increase B and rerun
- The cloglog transformation has issues at boundary probabilities
- Sample size is too small for the asymptotic approximation
- A bug in `cluster_boot`, `prodint_AJ`, or `ci_cloglog`

Document any deviations in the package (vignette, NEWS) and consider
delaying CRAN submission until the cause is understood.

## Output

Results are saved to `results/`:

- `stage1_summary.csv` -- one-line summary per check
- `stage1_full.rds` -- complete results including per-replicate
  point estimates, CI bounds, and band/test outcomes

## Beyond Stage 1

The full Bakoyannis (2021) simulation study covers:

- 6 sample-size cells (n in {20, 40, 80} x F_M in {U{5,15}, U{10,30}})
- 4 estimands (state occupation at tau_0.4, tau_0.6 for both ACM and
  TCM populations) plus transition probabilities and bands
- 2 censoring scenarios (right-cens only, right-cens + left-trunc)
- Both H_0 and H_1 for the two-sample test

Full reproduction is ~8-12 hours of compute time. The simulator
function `simulate_clusters()` supports all combinations via its
arguments (`cluster_size_range`, `two_sample`, `under_alternative`,
`left_truncate`). Stage 2 expansion involves writing a runner that
sweeps over the parameter grid; the simulator infrastructure here
is reusable.

## Reference

Bakoyannis, G. (2021). Nonparametric analysis of nonhomogeneous
multistate processes with clustered observations. *Biometrics*,
77(2), 533-546. doi:10.1111/biom.13327
