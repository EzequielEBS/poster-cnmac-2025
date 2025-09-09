library(mvtnorm)
library(parallel)
library(hdbayes)
library(cmdstanr)
library(dplyr)

# # load data
# data(airquality)
# hist_data <- airquality %>%
#   filter(Month <= 6)
# curr_data <- airquality %>%
#   filter(Month > 6)
# 
# # remove rows with missing values
# hist_data <- na.omit(hist_data)
# curr_data <- na.omit(curr_data)
# 
# # create response variable
# hist_data$logWind <- log(hist_data$Wind)
# curr_data$logWind <- log(curr_data$Wind)
# 
# # normalize predictor variables
# hist_data$Ozone <- (hist_data$Ozone - mean(hist_data$Ozone)) / sd(hist_data$Ozone)
# curr_data$Ozone <- (curr_data$Ozone - mean(curr_data$Ozone)) / sd(curr_data$Ozone)
# hist_data$Solar.R <- (hist_data$Solar.R - mean(hist_data$Solar.R)) / sd(hist_data$Solar.R)
# curr_data$Solar.R <- (curr_data$Solar.R - mean(curr_data$Solar.R)) / sd(curr_data$Solar.R)
# hist_data$Temp <- (hist_data$Temp - mean(hist_data$Temp)) / sd(hist_data$Temp)
# curr_data$Temp <- (curr_data$Temp - mean(curr_data$Temp)) / sd(curr_data$Temp)
# 
# # set parameters for Stan data
# formula     <- logWind ~ Ozone + Solar.R + Temp + Day
# family      <- gaussian()

load("data/sim_lm_data.RData")
formula <- y ~ X1
family <- gaussian()

sample_priorpred_pp <- function(eta) {
  res_hist          = hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
  res_curr          = hdbayes:::stack.data(formula = formula, data.list = list(curr_data))
  y0            = res_hist$y
  X0            = res_hist$X
  y = res_curr$y
  X = res_curr$X
  n <- length(y)
  p <- ncol(X)
  alpha <- 2
  beta <- 1
  mu_beta <- rep(0, p)
  S_beta <- diag(p)
  inv_S_beta <- solve(S_beta)
  S_eta <- inv_S_beta + eta * t(X0) %*% X0
  inv_S_eta <- solve(S_eta)
  mu_eta <- inv_S_eta %*% (inv_S_beta %*% mu_beta + eta * t(X0) %*% y0)
  mu <- X %*% mu_eta
  S <- beta/alpha * (diag(1, n) + X %*% inv_S_eta %*% t(X))
  
  m <- 10000
  tildey <- mvtnorm::rmvt(m, delta = mu, sigma = S, df = 2 * alpha)
  return(tildey)
}

ncores  <- 14
a0_list      <- seq(0, 1, length.out = 40)

# cl <- makeCluster(ncores)
# clusterSetRNGStream(cl, 123)
# clusterExport(cl, varlist = c('formula', 'family', 'hist_data'))
# clusterEvalQ(cl, {library("dplyr") 
#                  library("mvtnorm")})
# draws_priorpred <- parLapply(
#   cl = cl, X = a0_list, fun = sample_priorpred
# )
# stopCluster(cl)

draws_priorpred_pp <- mclapply(a0_list, FUN = sample_priorpred_pp,
                            mc.cores = ncores)

save(draws_priorpred_pp, file = 'samples_ppc/samples_priorpred_pp.RData')



#-------------------------------------------------------------------------------
# Sample NPP
#-------------------------------------------------------------------------------

sample_priorpred_npp <- function(eta) {
  res_hist          = hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
  res_curr          = hdbayes:::stack.data(formula = formula, data.list = list(curr_data))
  y0            = res_hist$y
  X0            = res_hist$X
  y = res_curr$y
  X = res_curr$X
  n <- length(y)
  p <- ncol(X)
  alpha <- 2
  beta <- 1
  mu_beta <- rep(0, p)
  S_beta <- diag(p)
  inv_S_beta <- solve(S_beta)
  S_eta <- inv_S_beta + eta * t(X0) %*% X0
  inv_S_eta <- solve(S_eta)
  mu_eta <- inv_S_eta %*% (inv_S_beta %*% mu_beta + eta * t(X0) %*% y0)
  mu <- X %*% mu_eta
  S <- beta/alpha * (diag(1, n) + X %*% inv_S_eta %*% t(X))
  tildey <- mvtnorm::rmvt(1, delta = mu, sigma = S, df = 2 * alpha)
  return(tildey)
}

m <- 10000
etas <- rbeta(m, 1, 1)

draws_priorpred_npp <- mclapply(etas, FUN = sample_priorpred_npp,
                               mc.cores = ncores)
draws_priorpred_npp <- do.call(rbind, draws_priorpred_npp)

save(draws_priorpred_npp, file = 'samples_ppc/samples_priorpred_npp.RData')
