pkgname <- "robcdcc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('robcdcc')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("rcpparma_hello_world")
### * rcpparma_hello_world

flush(stderr()); flush(stdout())

### Name: RcppArmadillo-Functions
### Title: Set of functions in example RcppArmadillo package
### Aliases: rcpparma_hello_world rcpparma_innerproduct
###   rcpparma_outerproduct rcpparma_bothproducts

### ** Examples

  x <- sqrt(1:4)
  rcpparma_innerproduct(x)
  rcpparma_outerproduct(x)



cleanEx()
nameEx("robcdcc-package")
### * robcdcc-package

flush(stderr()); flush(stdout())

### Name: robcdcc-package
### Title: Simulation and robust estimation for GARCH-cDCC process of order
###   (1,1)
### Aliases: robcdcc-package robcdcc
### Keywords: package

### ** Examples

~~ simple examples of the most important functions ~~



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
