# ENGE : Ensemble Dimensionality Reduction
# Author: Xiaoxiao Sun
# 08/26/2019

#default parameters
endr_defs <- list(
  n_wl = 500,
  n_dm = 15,
  n_neigs = 20,
  n_comps = 2,
  n_epos = 200,
  seed = 5489,
  alpha = 1,
  gamma = 1.0,
  a = NA,
  b = NA,
  spread = 1,
  min_dist = 0.1,
  negative_sample_rate = 5,
  H = 101107,
  opt = FALSE
)
class(endr_defs) <- "endr.config"

# main function
endr <- function(dat,config=endr_defs,...){
  config <- check_config(config,...)
  dat <- input_check(dat)
  if (config$n_dm > ncol(dat)){
    endr_error("'n_dm' must be less than the number of columns of the data")
  }
  if (is.na(config$a)|is.na(config$b)){
    config[c("a","b")] <- ab_params(config$spread,config$min_dist)
  }
  n_wl <- config$n_wl
  n_dm <- config$n_dm
  n_comps <- config$n_comps
  seed <- config$seed
  nk <- config$n_neigs
  H <- config$H
  if (dim(dat)[1] > H) {
    H <- 1017881
  }
  # graph representation
  amat <- affinity(dat,n_wl=n_wl,n_dm=n_dm,nk=nk,H=H,seed=seed)
  amat[,3] <- amat[,3]/(2*n_wl)
  # optimize embedding
  if(config$opt){
    amat[amat[,1]==amat[,2],3] <- 0.0
    colnames(amat) <- c("from","to","value")
    coor_knn <- make_coor(amat,rownames(dat),nrow(dat))
    amat_spec <- spectral_eigen(coor_knn,n_comps)
    opt_embedding <- endr_embed(coor_knn,amat_spec,config)
  } else {
    colnames(amat) <- c("from","to","value")
    coor_knn <- make_coor(amat,rownames(dat),nrow(dat))
    amat_spec <- spectral_eigen(coor_knn,n_comps)
    opt_embedding <- amat_spec
  }
  out_embedding <- center_embed(opt_embedding)
  class(out_embedding) <- "endr"
  out_embedding
}




