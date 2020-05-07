# EDGE: Ensemble Dimensionality Reduction 

## Prerequisites
    Rcpp >= 1.0.1
    C++11
## Install
You could install the package from the source file: EDGE_1.0.tar.gz
## Examples: embedding iris dataset
    library(EDGE)
    library(RSpectra)
    custom_defs <- endr_defs
### change the number of weak learners
    custom_defs$n_wl <- 150
### change the number of dimensions to be sampled
    custom_defs$n_dm <- 3
### change the number of nearest neighbors
    custom_defs$n_neigs <- 15
### optimize the embedding
    custom_defs$opt <- TRUE
    iris.endr <- endr(iris[,1:4],custom_defs)
    plot(iris.endr,col=iris[,5])
