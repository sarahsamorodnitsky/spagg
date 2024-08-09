
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spagg

<!-- badges: start -->
<!-- badges: end -->

The goal of `spagg` is to provide tools for aggregating spatial summary
statistics generated from multiple regions-of-interest (ROIs) imaged
using multiplexed spatial proteomics from the same tissue sample. 

`spagg` contains several methods for aggregating spatial summaries, 
including three types of weighted averages and three ensemble approaches. 
The ensemble approaches randomly generate weights to construct a weighted
average of the spatial summaries within each tissue sample. The p-values
across these replications are aggregated using the Cauchy combination test. 

## Installation

You can install the development version of spagg from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sarahsamorodnitsky/spagg")
```

## Vignettes

For example usage of `spagg` in single-cell spatial proteomics imaging
analysis, please see the associated vignettes.
