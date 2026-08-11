
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msmi

<!-- badges: start -->

<!-- badges: end -->

The goal of msmi is to use tools from multiple imputation to handle
missing data present in multistate models subject to censoring and
estimate state occupation probabilities from these data. Currently, msmi is
only suitable for the irreversible illness-death model.

## Installation

You can install the development version of msmi from
[GitHub](https://github.com/) with:

``` r
if (!require(remotes)) {install.packages("remotes")}
remotes::install_github("rtucker8/msmi")
```

## Usage

This package has three main functions. `msmi.impute` performs the imputation
algorithm and outputs a list of imputed datasets with no censoring. 
`msmi.tprobs` uses the imputed datasets to calculate the stat occupation
probabilities and the associated uncertainties at a set of requested times.
`msmi.plot` visualizes the returned estimates and associated confidence regions
on the probability simplex.

This is a basic example which shows you how to use the `msmi.impute`
function with simulated data `sim.data`:

``` r
library(msmi) 

head(sim.data)
#>   id         t1 event1         t2 event2  sojourn01 sojourn02 sojourn12
#> 1  1 0.84955504      0 0.84955504      0 2.97646713 34.359272  2.904316
#> 2  2 0.05423889      1 2.05011485      0 0.05423889  8.801891  2.282400
#> 3  3 2.14371939      1 4.91716494      1 2.14371939  7.834441  2.773446
#> 4  4 6.30786601      1 9.63997431      1 6.30786601 14.148366  3.332108
#> 5  5 2.50932926      1 4.49595865      1 2.50932926  8.355165  1.986629
#> 6  6 0.09055337      0 0.09055337      0 5.63360549  6.239124  2.886236

imps <- msmi.impute(sim.data, M = 5, n.states = 3, prefix.states = c("event", "t"), method = "marginal")
head(imps[[1]])
#>           t1 event1        t2 event2
#> 1 3.81009456      1  9.358730      1
#> 2 0.05423889      1  5.290110      1
#> 3 2.14371939      1  4.917165      1
#> 4 6.30786601      1  9.639974      1
#> 5 2.50932926      1  4.495959      1
#> 6 4.21332332      1 10.567040      1
```

Now we can use `msmi.tprobs` to estimate the state occupation probabilities at time 4 and 'msmi.plot` to visualize the results :

```r
tprobs <- msmi.tprobs.v2(imps, times = c(2,4))

print(trpobs$mi_estimate)

msmi.plot(tprobs, 4)
```

## References

This work is currently ongoing. A reference to the method paper will be added when available.


