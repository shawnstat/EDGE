# ENGE : Ensemble Dimensionality Reduction
# 08/26/2019

# get k nearest neighbors
endr_knn <- function(d,config){
  nk <- config$n_neigs
  nn <- sqrt(length(d))
  if (nk>nn) {
    endr.error("the number of neighbors must be smaller than number of cells")
  }
  if (nk<=1) {
    endr.error("the number of neighbors must be greater than 1")
  }
  distances <- big.matrix(ncol=nn,nrow=nn,init=0.0)
  afmat <- big.matrix(nrow=nn,ncol=nn,init=0.0)
  afmat[1:nn,1:nn] <- d
  for (i in 1:nn) {
    simi <- sort(afmat[i,],decreasing=TRUE)[1:nk]
    for(s in simi){
      j <- which(afmat[i,]%==%s)
      distances[i,j] <- afmat[i,j]
      distances[j,i] <- afmat[i,j]
    }
  }
#  rownames(distances) <- rownames(d)
  distances
}

# coordinate list
coor_list <- function(x){
  if (nrow(x)!=ncol(x)){
    stop("x must be a square matrix\n")
  }
  nx <- nrow(x)
  coor <- matrix(0,ncol=3,nrow=nx*nx)
  coor[,1] <- rep(1:nx,nx)
  coor[,2] <- rep(1:nx,each=nx)
  coor[,3] <- as.vector(x[1:nx,1:nx])
  colnames(coor) <- c("from","to","value")
  coor <- coor[order(coor[,1], coor[,2]),]
  make_coor(coor,rownames(x),nrow(x))
}

# make coordinate list
make_coor <- function(x,names,n_elems){
  x <- x[,1:3]
  colnames(x) <- c("from","to","value")
  result <- list(coor=x,names=names,n_elems=n_elems)
  class(result) <- "coor"
  result
}

#  remove zero entires
reduce_coor <- function(x){
  if (class(x)!="coor"){
    stop(paste0("expecting object of class coo!","\n"),call.=FALSE)
  }
  del_rows <- x$coor[,"value"]==0|!is.finite(x$coor[,"value"])
  x$coor <- x$coor[!del_rows,,drop=FALSE]
  rownames(x$coor) <- NULL
  x
}

# normalized Laplacian for a graph
lap_coor <- function(x){
  if (class(x)!="coor"){
    stop(paste0("expecting object of class coo ","\n"),call.=FALSE)
  }
  x <- reduce_coor(x)
  result <- x
  num_from <- length(unique(x$coor[,"from"]))
  if (num_from!=x$n_elems){
    stop(paste0("singular degrees!","\n"),call.=FALSE)
  }
  degrees <- x$coor[order(x$coor[,"from"]), ]
  degrees <- sapply(split(degrees[,"value"],degrees[,"from"]),sum)
  allelem <- x$coor
  oldvalues <- allelem[,"value"]
  allelem[,"value"] <- degrees[allelem[,"from"]] * degrees[allelem[,"to"]]
  allelem[,"value"] <- oldvalues/sqrt(allelem[,"value"])
  result$coor <- allelem
  result$coor <- result$coor[order(result$coo[,"from"],result$coor[,"to"]),]
  result
}

# get a set of k eigenvectors
spectral_eigen <- function(x,k){
  x_lap <- lap_coor(x)
  x_spa <- methods::new("dgTMatrix",
                          i=as.integer(x_lap$coor[,"from"]-1),
                          j=as.integer(x_lap$coor[,"to"]-1),
                          Dim=as.integer(rep(x_lap$n_elems,2)),
                          x=as.numeric(x_lap$coor[,"value"]))
  x_spa <- methods::as(x_spa, "dgCMatrix")
  result <- RSpectra::eigs_sym(x_spa,k,which="LA")$vectors
  rownames(result) <- x$names
  result
}

# compute a value to capture how often each item contributes to layout optimization
make_epochs <- function(w,epochs){
  result <- w
  result[1:length(w)] <- rep(-1,length(w))
  n_samps <- epochs*(w/max(w))
  n_pos <- n_samps>0
  result[n_pos] <- epochs/n_samps[n_pos]
  result
}

# estimate a/b parameters
ab_params <- function(spread,min_dist,alim=c(0,20),
                      blim=c(0,20),tolerance=1e-8){
  ## compute square error between two vectors
  sum_err <- function(x,y) {
    xy <- (x-y)
    sum(xy*xy)
  }
  ## compute y values given parameters a, b
  abcurve <- function(x,a,b) {
    xb2 <- x^(2*b)
    1.0/(1+(a*xb2))
  }
  ## create x values and target y values
  xv <- seq(0,spread*3,length=300)
  xv_high <- xv[xv>=min_dist]
  yv <- rep(0,length(xv))
  yv[xv<min_dist] <- 1
  yv[xv>min_dist] <- exp((min_dist-xv_high)/spread)
  ## internal recursive helper. Tries different combinations of a/b.
  ab_recursive <- function(alim,blim) {
    avals <- seq(alim[1],alim[2],length=10)
    bvals <- seq(blim[1],blim[2],length=10)
    ## compute square error of curve for all combinations of a/b
    ab <- expand.grid(avals,bvals)
    errors <- apply(ab,1,function(x) {
      yvals <- abcurve(xv,x[1],x[2])
      sum_err(yvals,yv)
    })
    ## identify combination with smallest error
    best <- as.numeric(ab[which.min(errors)[1],])
    ## determine if exit or keep looking in narrower interval
    mid <- c(mean(alim),mean(blim))
    if (sum(abs(best-mid))>tolerance){
      alim <- best[1] + (avals[2]-avals[1])*c(-1.5, 1.5)
      blim <- best[2] + (bvals[2]-bvals[1])*c(-1.5, 1.5)
      best <- ab_recursive(alim,blim)
    }
    best
  }
  result <- ab_recursive(alim,blim)
  names(result) <- c("a","b")
  result
}

# optimal embeddings
endr_embed <- function(g,embedding,config){
  if (config$n_epos==0) {
    return(embedding)
  }
  g <- reduce_coor(g)
  total_weight <- sum(g$coor[,"value"])
  ## simplify graph a little bit
  gmax <- max(g$coor[,"value"])
  g$coor[g$coor[,"value"]<gmax/config$n_epos,"value"] <- 0
  g <- reduce_coor(g)
  ## create an epochs-per-sample.
  eps <- cbind(g$coor,
              eps=make_epochs(g$coo[,"value"],config$n_epos))
  result <- opt_embed(embedding,config,eps)
  rownames(result) <- g$names
  result
}

# modify an existing embedding
opt_embed <- function(embedding,config,eps) {
  set.seed(config$seed)
  ## number of vertices in embedding
  V <- nrow(embedding)
  ## transpose to get observations in columns
  embedding <- t(embedding)
  ## define some vectors for book-keeping
  ## integer matrix with pairs of data
  eps_pairs <- matrix(as.integer(eps[,c("from","to")]),ncol=2)-1
  eps_val <- eps[,"eps"]
  ## epns is short for "epochs per negative sample"
  epns <- eps_val/config$negative_sample_rate
  ## eon2s is short for "epochs of next negative sample"
  eon2s <- epns
  ## eons is short for "epochs of next sample"
  eons <- eps_val
  ## nns is next negative sample
  nns <- rep(0,nrow(eps))
  ## extract some variables from config
  abg <- c(config$a,config$b,config$gamma)
  for (n in seq_len(config$n_epos)){
    ## set the learning rate for this epoch
    alpha <- config$alpha * (1 - ((n - 1)/config$n_epos))
    ## identify links in graph that require attention, then process those in loop
    adjust <- eons<=n
    ihits <- which(adjust)
    nns[ihits] <- floor((n-eon2s[ihits])/epns[ihits])
    embedding <- opt_epoch(embedding,eps_pairs,
                           as.integer(adjust),nns,abg,alpha)
    ## prepare for next epoch
    eons[ihits] <- eons[ihits] + eps_val[ihits]
    eon2s[ihits] <- eon2s[ihits] + (nns[ihits]*epns[ihits])
  }
  t(embedding)
}

# adjust a matrix so that each column is centered around zero
center_embed <- function(x){
  colcenters <- colMeans(x)
  V <- nrow(x)
  x - matrix(rep(colcenters,each=V),ncol=ncol(x),nrow=V)
}


