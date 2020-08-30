# EDGE

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
    library("umap")
    library("Rtsne")

## Simulate single cell RNA-seq data
We will use the simulator, Splatter, to simulate a single cell RNA-seq dataset with 500 cells and 200 genes. There are two groups in the dataset. The minor group contains 15 cells. 

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
We embed the normalized data into a two-dimensional subspace. Theoretically, the more weak learners and larger hash table sizes you use, the better result you will obtain. However, large memory will be required.  The optimal hash table size is 1,017,881. We reduce the hash table size to 101,107 in the example. Most of the time, this works well and uses less memory.  

    # change default parameters
    custom_defs <- endr_defs
    # the number of weak learners
    custom_defs$n_wl <- 1000
    # the number of genes selected to construct weak learners
    custom_defs$n_dm <- 10
    # the number of nearest neighbors
    custom_defs$n_neigs <- 15
    # low-dimensional embedding optimization
    custom_defs$opt <- TRUE
    # random seed
    custom_defs$seed <- 2020
    # the hash table size
    custom_defs$H <- 101107
    # EDGE
    simu.endr <- endr(dat_log,custom_defs)
    # UMAP: the number of nearest neighbors is also 15
    simu.umap <- umap(dat_log,random_state=2020)
    # t-SNE
    set.seed(2020)
    simu.tsne <- Rtsne(dat_log)
    # visualization
    op <- par(mfrow=c(1,3), bty="l",mar=c(1.5,1.5,1.5,1.5),xaxt="n",yaxt="n")
    plot(simu.endr, pch=16, cex=.8, col=labs)
    mtext(side=3,font=2,line=0,"EDGE")
    plot(simu.umap$layout, pch=16, cex=.8, col=labs)
    mtext(side=3,font=2,line=0,"UMAP")
    plot(simu.tsne$Y, pch=16, cex=.8, col=labs)
    mtext(side=3,font=2,line=0,"t-SNE")
    par(op)
