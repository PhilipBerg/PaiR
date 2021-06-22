
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pair

<!-- badges: start -->
<!-- badges: end -->

The goal of `pair` is to improve your proteomics data analysis with a
specific focus on PTM data. PTM data is usually noisier than total
proteomics and `pair` helps you normalize data, perform imputation, and
make more robust statistical decisions. It includes plots to help you
decide what normalization method to use and lets you visualize the
statistical decision. In addition, it uses a recently developed gamma
regression model to capture and normalize mean-variance trends in the
data. Further, it uses a multiple imputation pipeline that helps correct
for errors in the decision cause by imputation and generates a
probability of error in the statistical decision for imputed data. For
each individual imputation, it assumes that the mean can be correctly
estimated from the data and then uses the gamma regression to estimate
the variance given the mean.

## Installation

`pair` is still under development but will be released to CRAN shortly.

``` r
install.packages("pair")
```

In the meantime, you can download `pair` from github.

``` r
# install.packages("devtools")
devtools::install_github("PhilipBerg/pair")
```

## Example

For examples see the vignette.

``` r
library(pair)
```
