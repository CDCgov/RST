# The RST R Package
## Introduction

The Rate Stabilizing Tool (RST) package is a tool that uses a Bayesian spatiotemporal model in conjunction with C++ to help you quickly and easily generate spatially smoothed estimates for your spatial data. For the ArcGIS Plugin, visit the [RST ArcGIS Plugin webpage](https://www.cdc.gov/dhdsp/maps/gisx/rst.html).

### What is an MSTCAR model?

A Multivariate Spatiotemporal Conditional Auto-Regressive (MSTCAR) model is a Bayesian spatiotemporal model that generates spatially smoothed estimates across discrete regions. It is used to stabilize estimates in regions with low population counts, thereby bringing the inherent spatial structure of your event data to the forefront. It does this by taking information from various aspects of the data:
- The adjacency structure of your data. Neighboring regions consolidate information to stabilize estimates.
- Sociodemographic groups. The MSTCAR model also lends rate information from the sociodemographic groups included in the data.
- Time. The MSTCAR model helps to reveal trends in regions across time that may not be evident from the raw data.

RST is an R package is used to better understand binomial- or Poisson-distributed spatial count data by calculating spatially smoothed estimates from a Gibbs sampler.

### Motivation and Features

It is often very difficult to use an MSTCAR model in practice because it is 1. Very difficult to build without prior knowledge of spatiotemporal Bayesian hierarchical models, and 2. The complexity of the model makes it very computationally intensive. Because of this, There are several key ideas taken into consideration when creating the RST package:

- Designed to simply import count/population data for discrete regions and receive estimates smoothed by region, sociodemographic groups, and time periods.
- Designed to be easy-to-use and accessible, even for people without experience with Bayesian hierarchical models.
- Designed to be *fast*! Runs more than 10 times faster than equivalent models written only in R.

## Installation

To install the RST package, a few dependent packages are needed first. RST depends on the following:

##### Packages
- [Rcpp](https://cran.r-project.org/package=Rcpp): The most important package for RST. RST relies heavily on C++ code for efficiency. Rcpp is a package that allows integration of C++ code into R.
- [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo): This is an addon for Rcpp that facilitates matrix and array manipulation.
- [RcppDist](https://cran.r-project.org/package=RcppDist): This is an addon for Rcpp that allows the use of specialized distributions.
- [abind](https://cran.r-project.org/package=abind): A package designed to combine multidimensional arrays. Necessary for the concatenation of large arrays when gathering samples.
- [knitr](https://cran.r-project.org/package=knitr): Facilitates creation of package vignettes.
##### Programs
- [RTools](https://cran.r-project.org/bin/windows/Rtools/): C++ compiler for R. Necessary for Rcpp and its dependencies.

### Installation Instructions
The RST package can be easily installed with a few lines of R code. From the R console, make sure that all the necessary dependencies are installed by running the following lines:

```sh
# Install dependent packages
install.packages(c("Rcpp", "RcppArmadillo", "RcppDist", "abind", "knitr"))
# Install RTools
install.packages("installr") # optional: RTools can be downloaded and installed manually from the above URL
installr::install.Rtools()
```
Then, the package can be easily installed from GitHub:
```sh
install.packages("remotes")
remotes::install_github("CDCgov/RST", build_vignettes = TRUE)
```
If this is your first time using the RST R package, check the introductory vignette to learn how to use the package:
```sh
vignette("RST")
```

### Updates
Updating the RST package is fairly simple. All that has to be done is remove the package from R, restart R, and re-install from GitHub:
```sh
remove.packages("RST")
# Restart R to allow package installation
remotes::install_github("CDCgov/rst", build_vignettes = TRUE)
```

## Beta
This package is currently in the beta stage. Each public function and dataset has documentation to help orient you to using the package, along with an introductory vignette to walk you through a simple analysis using the example data. For debugging purposes, private functions can be accessed with RST:::problem_function().

## Thanks!
Thanks for checking out the RST package. I hope you find it useful in your spatiotemporal modeling! Please feel free to give feedback on bugs, ways to make the package more accessible to use, and features you'd like to see added to the package.
