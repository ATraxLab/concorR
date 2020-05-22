
<!-- README.md is generated from README.Rmd. Please edit that file -->

# concorR

<!-- badges: start -->

<!-- badges: end -->

The goal of concorR is to implement the CONCOR (CONvergence of iterated
CORrelations) algorithm for positional analysis. Positional analysis
divides a network into blocks based on the similarity of links between
actors. CONCOR uses structural equivalence—“same ties to same others”—as
its criterion for grouping nodes, and calculates this by correlating
columns in the adjacency matrix. For more details on CONCOR, see the
original description by Breiger, Boorman, and Arabie (1975), or Chapter
9 in Wasserman and Faust (1994).

## Installation

You can install the released version of concorR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("concorR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("atraxler/concorR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(concorR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!

## References

R. L. Breiger, S. A. Boorman, P. Arabie, An algorithm for clustering
relational data with applications to social network analysis and
comparison with multidimensional scaling. *J. of Mathematical
Psychology*. **12**, 328 (1975).
<http://doi.org/10.1016/0022-2496(75)90028-0>

S. Wasserman and K. Faust, *Social Network Analysis: Methods and
Applications* (Cambridge University Press, 1994).
