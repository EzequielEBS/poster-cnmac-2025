# wrapper function for getting Stan data for sampling from NPP
library(hdbayes)

get.stan.data.npp.prior = function(
    formula,
    family,
    data.list, # list of external/historical data sets
    a0.lognc,
    lognc,
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL
) {
  res          = hdbayes:::stack.data(formula = formula, data.list = data.list)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  p            = ncol(X)
  N            = length(y)
  K            = length(end.index)
  fam.indx     = hdbayes:::get.dist.link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]

  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }

  ## Default prior on regression coefficients is N(0, 10^2)
  if ( !is.null(beta.mean) ){
    if ( !( is.vector(beta.mean) & (length(beta.mean) %in% c(1, p)) ) )
      stop("beta.mean must be a scalar or a vector of length ", p, " if beta.mean is not NULL")
  }
  beta.mean = hdbayes:::to.vector(param = beta.mean, default.value = 0, len = p)
  if ( !is.null(beta.sd) ){
    if ( !( is.vector(beta.sd) & (length(beta.sd) %in% c(1, p)) ) )
      stop("beta.sd must be a scalar or a vector of length ", p, " if beta.sd is not NULL")
  }
  beta.sd = hdbayes:::to.vector(param = beta.sd, default.value = 10, len = p)

  ## check a0.lognc and lognc
  if ( length(a0.lognc) != nrow(lognc) )
    stop('the number of rows in lognc must be the same as the length of a0.lognc')
  if ( ncol(lognc) != K )
    stop('the number of columns in lognc must be the same as the number of historical data sets')
  if ( any(is.na(a0.lognc) ) )
    stop('a0.lognc must not have missing values')
  if ( any(is.na(lognc)) )
    stop('lognc must not have missing values')
  if ( any(a0.lognc < 0) || any(a0.lognc > 1) )
    stop('each element of a0.lognc should be between 0 and 1')

  ## Default half-normal prior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp.mean) ){
    if ( !( is.vector(disp.mean) & (length(disp.mean) == 1) ) )
      stop("disp.mean must be a scalar if disp.mean is not NULL")
  }
  disp.mean = hdbayes:::to.vector(param = disp.mean, default.value = 0, len = 1)
  if ( !is.null(disp.sd) ){
    if ( !( is.vector(disp.sd) & (length(disp.sd) == 1) ) )
      stop("disp.sd must be a scalar if disp.sd is not NULL")
  }
  disp.sd = hdbayes:::to.vector(param = disp.sd, default.value = 10, len = 1)

  ## Default lower bound for each a0 is 0; default upper bound for each a0 is 1
  if ( !is.null(a0.lower) ){
    if ( !( is.vector(a0.lower) & (length(a0.lower) %in% c(1, K)) ) )
      stop("a0.lower must be a scalar or a vector of length ", K, " if a0.lower is not NULL")
  }
  a0.lower = hdbayes:::to.vector(param = a0.lower, default.value = 0, len = K)
  if ( !is.null(a0.upper) ){
    if ( !( is.vector(a0.upper) & (length(a0.upper) %in% c(1, K)) ) )
      stop("a0.upper must be a scalar or a vector of length ", K, " if a0.upper is not NULL")
  }
  a0.upper = hdbayes:::to.vector(param = a0.upper, default.value = 1, len = K)

  standat = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'mean_beta'       = beta.mean,
    'sd_beta'         = beta.sd,
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    's'               = length(a0.lognc),
    'a0_lognc'        = a0.lognc,
    'lognc'           = lognc,
    'a0_shape1'       = a0.shape1,
    'a0_shape2'       = a0.shape2,
    'a0_lower'        = a0.lower,
    'a0_upper'        = a0.upper,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )
  return(standat)
}

glm_npp_prior = cmdstan_model("code/glm_npp_prior.stan")


get.stan.data.pp.prior = function(
    formula,
    family,
    data.list, # list of external/historical data sets
    offset.list       = NULL,
    beta.mean         = NULL,
    beta.sd           = NULL,
    disp.mean         = NULL,
    disp.sd           = NULL,
    a0.lower          = NULL,
    a0.upper          = NULL,
    a0s
) {
  res          = hdbayes:::stack.data(formula = formula, data.list = data.list)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  p            = ncol(X)
  N            = length(y)
  K            = length(end.index)
  fam.indx     = hdbayes:::get.dist.link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]
  
  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }
  
  ## Default prior on regression coefficients is N(0, 10^2)
  if ( !is.null(beta.mean) ){
    if ( !( is.vector(beta.mean) & (length(beta.mean) %in% c(1, p)) ) )
      stop("beta.mean must be a scalar or a vector of length ", p, " if beta.mean is not NULL")
  }
  beta.mean = hdbayes:::to.vector(param = beta.mean, default.value = 0, len = p)
  if ( !is.null(beta.sd) ){
    if ( !( is.vector(beta.sd) & (length(beta.sd) %in% c(1, p)) ) )
      stop("beta.sd must be a scalar or a vector of length ", p, " if beta.sd is not NULL")
  }
  beta.sd = hdbayes:::to.vector(param = beta.sd, default.value = 10, len = p)
  
  # ## check a0.lognc and lognc
  # if ( length(a0.lognc) != nrow(lognc) )
  #   stop('the number of rows in lognc must be the same as the length of a0.lognc')
  # if ( ncol(lognc) != K )
  #   stop('the number of columns in lognc must be the same as the number of historical data sets')
  # if ( any(is.na(a0.lognc) ) )
  #   stop('a0.lognc must not have missing values')
  # if ( any(is.na(lognc)) )
  #   stop('lognc must not have missing values')
  # if ( any(a0.lognc < 0) || any(a0.lognc > 1) )
  #   stop('each element of a0.lognc should be between 0 and 1')
  
  ## Default half-normal prior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp.mean) ){
    if ( !( is.vector(disp.mean) & (length(disp.mean) == 1) ) )
      stop("disp.mean must be a scalar if disp.mean is not NULL")
  }
  disp.mean = hdbayes:::to.vector(param = disp.mean, default.value = 0, len = 1)
  if ( !is.null(disp.sd) ){
    if ( !( is.vector(disp.sd) & (length(disp.sd) == 1) ) )
      stop("disp.sd must be a scalar if disp.sd is not NULL")
  }
  disp.sd = hdbayes:::to.vector(param = disp.sd, default.value = 10, len = 1)
  
  ## Default lower bound for each a0 is 0; default upper bound for each a0 is 1
  if ( !is.null(a0.lower) ){
    if ( !( is.vector(a0.lower) & (length(a0.lower) %in% c(1, K)) ) )
      stop("a0.lower must be a scalar or a vector of length ", K, " if a0.lower is not NULL")
  }
  a0.lower = hdbayes:::to.vector(param = a0.lower, default.value = 0, len = K)
  if ( !is.null(a0.upper) ){
    if ( !( is.vector(a0.upper) & (length(a0.upper) %in% c(1, K)) ) )
      stop("a0.upper must be a scalar or a vector of length ", K, " if a0.upper is not NULL")
  }
  a0.upper = hdbayes:::to.vector(param = a0.upper, default.value = 1, len = K)
  
  standat = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'mean_beta'       = beta.mean,
    'sd_beta'         = beta.sd,
    'disp_mean'       = disp.mean,
    'disp_sd'         = disp.sd,
    'a0_lower'        = a0.lower,
    'a0_upper'        = a0.upper,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset,
    'a0s'             = a0s
  )
  return(standat)
}

glm_pp_prior = cmdstan_model("code/glm_pp_prior.stan")
