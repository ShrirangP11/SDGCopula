# SDGCopula
<!-- badges: start -->
<!-- badges: end -->
R Package to generate synthetic data from tabular input data consisting of continuous, integer or categorical variables using copula modelling techniques.
## Installation
You can install this package with using
``` r
remotes::install_github("ShrirangP11/SDGCopula")
```
## Example
To use the function provided by this package, run the following code.
```r
library(SDGCopula)
original <- iris
synthetic <- fitCop(original, copula='normal', parametric=FALSE)
```

