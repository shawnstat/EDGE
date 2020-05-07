\docType{data}
\name{endr_defs}
\alias{endr_defs}
\title{Default parameters for endr}
\format{An object of class \code{endr.config}}
\usage{
endr_defs
}
\description{
A list with parameters customizing a endr embedding
}
\details{
n_wl: integer; number of weak learners

n_dm: integer; number of dimensions to be sampled

n_neigs: integer; number of nearest neighbors

n_comps: integer; dimension of target (output) space

n_epos: integer; number of iterations performed during
layout optimization

seed: integer; seed for random number generator

alpha: numeric; initial value of "learning rate" of layout optimization

gamma: numeric; determines, together with alpha, the learning rate of layout optimization

a: numeric; contributes to gradient calculations during layout optimization.
When left at NA, a suitable value will be estimated automatically.

b: numeric; contributes to gradient calculations during layout optimization.

spread: numeric; used during automatic estimation of a/b parameters.

min_dist: numeric; determines how close points appear in the final layout

negative_sample_rate: integer; determines how many non-neighbor points are
used per point and per iteration during layout optimization

opt: logical; determines whether to perform layout optimization
}
\examples{
# display all default parameters
endr_defs
}
\keyword{datasets}
