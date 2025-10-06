
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msmi

<!-- badges: start -->

<!-- badges: end -->

The goal of msmi is to use tools from multiple imputation to handle
missing data present in multistate models subject to censoring and
estimate state occupation probabilities from these data.

## Installation

You can install the development version of msmi from
[GitHub](https://github.com/) with:

``` r
if (!require(remotes)) {install.packages("remotes")}
remotes::install_github("elizabethchase/mici")
remotes::install_github("rtucker8/msmi")
```

## Usage

This is a basic example which shows you how to use the `msmi.impute`
function with simulated data `sim.data`:

``` r
library(msmi)
head(sim.data)
#>   id       t1 event1        t2 event2 sojourn01 sojourn02 sojourn12
#> 1  1 1.082190      1  2.205080      0  1.082190  9.087453  3.463207
#> 2  2 3.725816      1 11.119825      0  3.725816  5.008579 15.416609
#> 3  3 1.961939      0  1.961939      0  3.955065 19.201314  4.318890
#> 4  4 2.507155      0  2.507155      1  7.482763  2.507155        NA
#> 5  5 2.422337      1  3.928518      0  2.422337 10.997895 10.222074
#> 6  6 4.017424      0  4.017424      1  4.879177  4.017424        NA
head(msmi.impute(sim.data, M = 5, n.states = 3, prefix.states = c("event", "t"), method = "marginal")[[1]])
#>          t1 event1        t2 event2
#> 1  1.082190      1  3.758061      1
#> 2  3.725816      1 14.185732      1
#> 3  7.756142      1 10.777875      1
#> 4  2.986981      1 13.446897      1
#> 5  2.422337      1 12.882253      1
#> 6 10.631872      0 10.631872      1
```
