\name{endr}
\alias{endr}
\title{Ensemble dimensionality reduction}
\usage{
endr(dat, config = endr_defs, ...)
}
\arguments{
\item{dat}{matrix, input data}

\item{config}{object of class endr.config}

\item{...}{list of settings; overwrite default values from config}
}
\value{
object of class endr, containing an embedding
}
\description{
Spectral embedding of the data
}
\examples{
library(RSpectra)
library(splatter)
# simulated data
myseed <- 2020
BC <- 500
GC <- 200
gp <- c(0.98,0.02)
params <- newSplatParams(nGenes=GC, seed=myseed, batchCells=BC, group.prob=gp)
simu_dt <- splatSimulate(params=params,method="groups")
counts_dt <- counts(simu_dt)
labs <- colData(simu_dt)[,3]
rs <- colSums(counts_dt)
rs_med <- median(rs)
dat_norm <- sweep(counts_dt,2,rs/rs_med,"/")
dat_log <- t(log2(dat_norm+1))
# embedd simulated data
custom_defs <- endr_defs
# change the number of weak learners
custom_defs$n_wl <- 1000
# change the number of dimensions to be sampled
custom_defs$n_dm <- 10
# change the number of nearest neighbors
custom_defs$n_neigs <- 15
# optimize the embedding
custom_defs$opt <- TRUE
custom_defs$seed <- 2020
simu.endr <- endr(dat_log,custom_defs)
plot(simu.endr,col=labs)
}
