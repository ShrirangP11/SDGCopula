# SDGCopula
<!-- badges: start -->
<!-- badges: end -->
R Package to generate synthetic data from tabular input data consisting of continous, interger or categorical using copula modelling techniques.
## Installation
You can install this package with using
``` r
remotes::install_github("ShrirangP11/SDGCopula")
```
## ExampleTo use the function provided by this package, run the following code.
```r
library(SDGCopula)
original <- iris
synthetic <- fitCop(original, copula='normal', parametric=FALSE)
```
## License
This package is licensed under Apache 2.0 license.
