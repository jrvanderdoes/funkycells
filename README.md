
<!-- README.md is generated from README.Rmd. Please edit that file -->

# funkycells

<!-- badges: start -->

[![R-CMD-check](https://github.com/jrvanderdoes/funkycells/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jrvanderdoes/funkycells/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The term ${\tt funkycells}$ comes from **fun**ctional data anal**y**sis
of **K** functions+ for multiplexed images of **cells**. This package
organizes ways to analyze cell relationships based on their (empirical)
K functions, or other such functions. The approach achieves effective
analysis for many different data constructions and accounts for issues
such as overfitting. We encourage all feedback and improvements, which
can be submitted through this github repo or the package
[website](https://jrvanderdoes.github.io/funkycells/).

Please see the package
[website](https://jrvanderdoes.github.io/funkycells/) for an
introduction to ${\tt funkycells}$, given in the vignette
`vignette("funkycells")`. Additional vignettes are also present showing
applications and simulated performance. The website also contains the
change log, documenting changes for each version of the package.

## Installation

The stable version of ${\tt funkycells}$ can be found on CRAN. Install
using:

``` r
install.packages("funkycells")
```

You can also install the development version of ${\tt funkycells}$ from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("jrvanderdoes/funkycells")
```

We are actively developing the package. So while new functionality
consistently pops up in the development version, we regularly update the
CRAN version with the stable additions.

<!-- 
Don't forget to build this! devtools::build_readme() , also check out https://github.com/r-lib/actions/tree/v1/examples to setup github actions for it 
-->
