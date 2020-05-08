# EDGE: Ensemble Dimensionality Reduction for Single Cell RNA-seq Data

## Prerequisites
    Rcpp >= 1.0.1
    C++11
    RSpectra
    
## Install
    library("devtools")
    install_github("shawnstat/EDGE")

## Load required packages
    library("RSpectra")
    library("splatter")
    library("EDGE")

## Simulate single cell RNA-seq data
We will use the simulator, splatter, to simulate a single cell RNA-seq dataset with 500 cells and 200 genes. There are two groups in the dataset. The minor group contains 15 cells. 

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

## Embed simulated data
We embed the normalized data into a two-dimensional subspace. We will use 1,000 weak lerners. 

    custom_defs <- endr_defs
    custom_defs$n_wl <- 1000
    custom_defs$n_dm <- 10
    custom_defs$n_neigs <- 15
    custom_defs$opt <- TRUE
    custom_defs$seed <- 2020
    simu.endr <- endr(dat_log,custom_defs)
    plot(simu.endr,col=labs)
