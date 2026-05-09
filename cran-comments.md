# cran-comments.md

## Test environments

* local: macOS Sonoma 14.x, R 4.4.x
* win-builder (devtools::check_win_devel()): R-devel
* R-hub (rhub::rhub_check()): Linux (Ubuntu 22.04, R-release and R-devel)
* GitHub Actions: ubuntu-latest, macos-latest, windows-latest

## R CMD check results

0 errors | 0 warnings | 0 notes

## This is a new release

This is the first submission of clusteredMSM.

The package implements methods from Bakoyannis (2021)
<doi:10.1111/biom.13327> for nonparametric analysis of clustered
multistate process data. It is intended as the canonical software
implementation accompanying that paper, and extends it to support
non-monotone (recovery) multistate processes.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
