# clusteredMSM

Nonparametric analysis of clustered multistate process data.

`clusteredMSM` provides population-averaged transition probability
estimates, pointwise confidence intervals, simultaneous confidence
bands, and two-sample Kolmogorov-Smirnov-type tests for multistate
process data with cluster-correlated observations. Methods are based
on [Bakoyannis (2021)](https://doi.org/10.1111/biom.13327) and use the
working-independence Aalen-Johansen estimator with a cluster-bootstrap
variance.

Unlike its predecessor (the `clustered-multistate` repository, which
relied on the `mstate` package), `clusteredMSM` is self-contained
(depending only on `survival`) and supports **non-monotone multistate
processes**, including illness-death with recovery and other models
with cyclic transitions.

## Installation

```r
# install.packages("devtools")
devtools::install_github("gbakoyannis/clusteredMSM")
```

After CRAN release:

```r
install.packages("clusteredMSM")
```

## A single function with a formula interface

`clusteredMSM` exposes one main function, `patp()`, modelled after
`survival::Surv()`:

```r
library(clusteredMSM)

# Synthetic clustered illness-death-with-recovery data (40 subjects,
# 8 clusters); see ?example_msm.
data(example_msm)

# Define the transition structure (illness-death with recovery)
tmat <- trans_mat(list(c(2, 3), c(1, 3), integer(0)),
                  names = c("Healthy", "Ill", "Dead"))

# One-sample analysis: P(Ill at t | Healthy at 0)
fit <- patp(msm(Tstart, Tstop, Sstart, Sstop) ~ 1,
            data = example_msm, tmat = tmat,
            id = "id", cluster = "cluster",
            h = 1, j = 2, s = 0,
            B = 1000, cband = TRUE)
fit
```

If the formula's right-hand side has a grouping variable, `patp()`
automatically estimates both group-specific curves AND tests their
equality:

```r
# Two-sample analysis (estimate + test in one call)
tt <- patp(msm(Tstart, Tstop, Sstart, Sstop) ~ treatment,
           data = example_msm, tmat = tmat,
           id = "id", cluster = "cluster",
           h = 1, j = 2, B = 1000)
tt
```

### Loading your own data

The same example is shipped as a CSV under `inst/extdata/`, so you can
mimic the typical workflow of reading a user-supplied file:

```r
f <- system.file("extdata", "example_data.csv", package = "clusteredMSM")
mydata <- read.csv(f)
head(mydata)
```

## Input data format

Each row of your data represents one **mutually-exclusive time
interval** for one subject, with columns:

| Column     | Description                                           |
|------------|-------------------------------------------------------|
| `Tstart`   | Numeric start time of the interval                    |
| `Tstop`    | Numeric end time of the interval                      |
| `Sstart`   | Integer state occupied during the interval            |
| `Sstop`    | Integer state at `Tstop` (or equal to `Sstart` if censored) |
| `id`       | Subject identifier                                    |
| `cluster`  | (optional) cluster identifier                         |
| (group)    | (optional) binary grouping variable                   |

The column names are arbitrary — `msm(...)` and the `id`/`cluster`
arguments tell the package which is which.

**Censoring** is encoded as `Sstart == Sstop` on the **final row** of a
subject's record. Subjects in absorbing states have no row after them.

Within each subject, intervals must be:
- **Temporally contiguous:** `Tstop[k] == Tstart[k+1]`
- **Spatially contiguous:** `Sstop[k] == Sstart[k+1]`

Validation is strict and informative — any violation triggers an error
with a clear message.

## Examples of valid input

**Progressive illness-death (subject who got ill, then died):**

| id | Tstart | Tstop | Sstart | Sstop |
|----|--------|-------|--------|-------|
| 1  | 0.0    | 1.5   | 1      | 2     |
| 1  | 1.5    | 3.0   | 2      | 3     |

**Subject censored healthy:**

| id | Tstart | Tstop | Sstart | Sstop |
|----|--------|-------|--------|-------|
| 2  | 0.0    | 4.0   | 1      | 1     |

**Recovery (Healthy → Ill → Healthy → censored):**

| id | Tstart | Tstop | Sstart | Sstop |
|----|--------|-------|--------|-------|
| 3  | 0.0    | 1.0   | 1      | 2     |
| 3  | 1.0    | 2.0   | 2      | 1     |
| 3  | 2.0    | 3.5   | 1      | 1     |

## Core functions

| Function | Purpose |
|---|---|
| `patp()` | The main user-facing function — formula-based estimation and testing. |
| `msm()` | Constructor for multistate intervals; used inside the formula. |
| `trans_mat()` | Build a K x K transition matrix. |
| `validate_intervals()` | Validate user data (called automatically by `patp()`; usable directly). |

## Citation

If you use `clusteredMSM` in your work, please cite:

> Bakoyannis, G. (2021). Nonparametric analysis of nonhomogeneous
> multistate processes with clustered observations. *Biometrics*,
> 77(2), 533-546. doi:10.1111/biom.13327

You can retrieve the BibTeX entry within R via `citation("clusteredMSM")`.

## License

GPL-3
