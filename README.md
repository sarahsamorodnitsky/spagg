
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spagg

<!-- badges: start -->
<!-- badges: end -->

The goal of spagg is to provide tools for aggregating spatial summary
statistics generated from multiple regions-of-interest (ROIs) imaged
using spatial proteomics from the same tissue sample.

## Installation

You can install the development version of spagg from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sarahsamorodnitsky/spagg")
```

## Example

``` r
library(spagg)
```

Let’s consider analyzing a non-small cell lung cancer dataset generated
from a multiplexed immunohistochemistry study by Johnson et al. (2021).
Johnson et al. found that CD4+ T cells and tumor cells colocalized
differently in major histocompatibility complex II (MHCII)-high tumors
than in MHCII-low tumors. We are going to test the same hypothesis using
bivariate Ripley’s K as our measure of spatial colocalization between
CD4+ T cells and tumor cells.

Since this dataset consists of multiple ROIs per tumor sample of $n=153$
tumors, we need to aggregate Ripley’s K within a sample. `spagg`
provides several ways of doing so.

First, let’s load in the data we will use.

``` r
data(lung_df_tumor)
```
