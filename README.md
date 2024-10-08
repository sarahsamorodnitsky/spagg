
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spagg

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/spagg)](https://CRAN.R-project.org/package=spagg)
[![R-CMD-check](https://github.com/sarahsamorodnitsky/spagg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sarahsamorodnitsky/spagg/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The goal of `spagg` is to provide tools for aggregating spatial summary
statistics generated from multiple regions-of-interest (ROIs) collected 
from the same tissue sample using multiplexed spatial proteomics technologies.

The analysis of spatial proteomics often involves calculating a spatial summary
statistic, such as Ripley's K, to quantify the level of clustering, repulsion, or 
complete spatial randomness exhibited by the cells. Given multiple ROIs for each sample,
aggregating the spatial summary statistics for each ROI into a single value can facilitate
downstream association testing with clinical outcomes. 

`spagg` provides several methods for aggregating spatial summary statistics. These include
three weighted means (`diggle.avg`, `baddeley.avg`, and `landau.avg`) and two ensemble
approaches (`ensemble.avg` and`combo.weight.avg`). The weighted means aggregate
the summary statistics using a weighted mean based on the number of cells in each ROI and/or the 
area of each ROI. The ensemble approaches use random weights to construct an aggregation and combine the 
resulting p-values across many randomly-generated weights for an ensemble test. 

The weighted means can be easily incorporated into other analytical approaches for spatial proteomics,
such as the SPatial Omnibus Test (SPOT). Incorporation of the ensemble approaches in SPOT is in development. 

## Installation

You can install the development version of spagg from
[GitHub](https://github.com/) with:

``` r
# First, install devtools
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install from Github
devtools::install_github("sarahsamorodnitsky/spagg")
```

`spagg` relies on several required dependencies: `ACAT`, `dplyr`, `magrittr`, the `spatstat` family of packages, `stats`, `survival`, and `tidyselect`. The `ACAT` package is currently in development on Github. To properly install `spagg`, the latest version of `devtools` is required to ensure it can install `ACAT` from Github. 

In addition, `spagg` relies on several suggested packages to run the vignettes: `knitr`, `rmarkdown`, `ggplot2`, `tidyr`, `spatstat.utils`, `spatstat.univar`, `cowplot`, and `SPOT`. `SPOT` is also available on Github, so the latest version of `devtools` is required. 

## Vignettes

For example usage of `spagg` in single-cell spatial proteomics imaging
analysis, please see the associated vignettes. The `spagg` package contains three vignettes:

1. `Getting Started`: this vignette illustrates how to use `spagg` for univariate or bivariate colocalization analyses on simulated data.
2. `Analysis of a Non-Small Cell Lung Cancer Dataset`: this vignette illustrates how to use `spagg` to analyze a non-small cell lung cancer dataset.
3. `Accommodating Multiple Radii`: this vignette illustrates how to incorporate multiple radii into a `spagg` analysis using the `SPOT` method.

## Bugs and Improvements

Please feel free to use the `Issues` tab on [the spagg Github site](https://github.com/sarahsamorodnitsky/spagg) to note bugs or to suggest improvements. 
