# SDGCopula<!-- badges: start -->
<!-- badges: end -->This package links to an article about creating R packages.## InstallationYou can install this package with using``` r
remotes::install_github("ShrirangP11/SDGCopula")
```## ExampleTo use the function provided by this package, run the following code.```r
library(SDGCopula)
original <- iris
synthetic <- fitCop(original, copula='normal', parametric=FALSE)
```## LicenseThis package is licensed under Apache 2.0 license.
