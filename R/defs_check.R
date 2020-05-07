# ENGE : Ensemble Dimensionality Reduction
# 08/26/2019

# check on configurations
check_config <- function(config=endr_defs,...){
  ## check on config class
  if (class(config)!="endr.config") {
    endr_error("config is absent or corrupt")
  }
  ## more arguments
  arguments <- list(...)
  for (args in names(arguments)){
    config[[args]] <- arguments[[args]]
  }
  ## check on missing arguments
  missing <- setdiff(names(endr_defs), names(config))
  if (length(missing)>0) {
    endr_error(paste0("missing arguments: ",paste(missing,collapse=", ")))
  }
  ## checks on individual parameters
  if (!is.finite(config$n_wl)|config$n_wl<0) {
    endr_error("the number of weak leaners must be positive")
  }
  if (!is.finite(config$n_dm)|config$n_dm<0) {
    endr_error("the number of selected genes must be positive")
  }
  if (config$n_neigs<2) {
    endr_error("the number of neighbors must be greater than 1")
  }
  if (!is.finite(config$n_epos)|config$n_epos<0) {
    endr_error("the number of epochs must be positive")
  }
  if (!is.finite(config$seed)|config$seed<0) {
    endr_error("the seed must be positive")
  }
  # checks on embedding control via (a, b) vs. (min_dist, spread)
  if (!identical(config$a,NA)&identical(config$b,NA)) {
    endr_warning("parameter 'a' is set but 'b' is not;\n",
                 "value of 'a' will be ignored.\n",
                 "(embedding will be controlled via 'min_dist' and 'spread')")
  }

  if (!identical(config$b,NA)&identical(config$a,NA)) {
    endr_warning("parameter 'b' is set but 'a' is not;\n",
                 "value of 'b' will be ignored.\n",
                 "(embedding will be controlled via 'min_dist' and 'spread')")
  }
  if (!identical(config$a,NA)&!identical(config$b,NA)) {
    abcontrol <- "(embedding will be controlled via 'a' and 'b')"
    if (!identical(config$min_dist,endr_defs$min_dist)) {
      endr_warning("parameters 'min_dist', 'a', 'b' are set to non-default values;\n",
                   "parameter 'min_dist' will be ignored.\n",abcontrol)
    }
    if (!identical(config$spread,endr_defs$spread)) {
      endr_warning("parameters 'spread', 'a', 'b' are set to non-default values;\n",
                   "parameter 'spread' will be ignored.\n",abcontrol)
    }
  }
  ## check on 'min_dist' and 'spread'
  if (config$min_dist>=config$spread) {
    endr_error("setting 'min_dist' must be smaller than 'spread'")
  }
  if (config$min_dist<=0) {
    endr_error("setting 'min_dist' must be > 0")
  }
  ## force some data types
  for (x in c("n_wl","n_ge","n_neigs","n_comps","n_epos","seed")){
    config[[x]] <- as.integer(config[[x]])
  }
  config
}

# check on data
input_check <- function(dat) {
  # should be matrix
  if (class(dat)=="matrix") {
    dat <- dat
  } else if (sum(class(dat)%in%c("data.frame","data.table"))) {
    dat <- as.matrix(dat)
  } else {
    endr_error("input must be a matrix or matrix-compatible\n")
  }
  dat
}

# error message
endr_error <- function(...) {
  x <- paste(...,collapse=" ")
  stop(paste0("endr: ",x,"\n"),call.=FALSE)
}

# warning message
endr_warning <- function(...) {
  x <- paste(...,collapse=" ")
  warning(paste0("endr: ",x,"\n"),call.=FALSE)
}


