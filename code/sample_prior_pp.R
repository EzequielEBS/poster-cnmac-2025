# example for sampling from the PP (prior)
library(cmdstanr)
library(hdbayes)
library(parallel)
library(tidyverse)
library(dplyr)
library(patchwork)

# load wrapper function for getting Stan data for sampling from PP
source("code/get_stan_data.R")

# load data
data(airquality)
hist_data <- airquality %>%
  filter(Month <= 6)
curr_data <- airquality %>%
  filter(Month > 6)

# remove rows with missing values
hist_data <- na.omit(hist_data)
curr_data <- na.omit(curr_data)

# create response variable
hist_data$logWind <- log(hist_data$Wind)
curr_data$logWind <- log(curr_data$Wind)

# normalize predictor variables
hist_data$Ozone <- (hist_data$Ozone - mean(hist_data$Ozone)) / sd(hist_data$Ozone)
curr_data$Ozone <- (curr_data$Ozone - mean(curr_data$Ozone)) / sd(curr_data$Ozone)
hist_data$Solar.R <- (hist_data$Solar.R - mean(hist_data$Solar.R)) / sd(hist_data$Solar.R)
curr_data$Solar.R <- (curr_data$Solar.R - mean(curr_data$Solar.R)) / sd(curr_data$Solar.R)
hist_data$Temp <- (hist_data$Temp - mean(hist_data$Temp)) / sd(hist_data$Temp)
curr_data$Temp <- (curr_data$Temp - mean(curr_data$Temp)) / sd(curr_data$Temp)

# set parameters for Stan data
formula     <- logWind ~ Ozone + Solar.R + Temp + Day
family      <- gaussian()

# define stuff to compute log normalising constant
ncores  <- 14
a0_list      <- seq(0.05, 1, length.out = 20)

sample_pp <- function(a0){

  standat <- get.stan.data.pp.prior(
    formula        = formula,
    family         = family,
    data.list      = list(hist_data),
    offset.list    = NULL,
    beta.mean      = 0,
    beta.sd        = 10,
    disp.mean      = NULL,
    disp.sd        = NULL,
    a0.lower       = 0,
    a0.upper       = 1,
    a0s            = c(a0)
  )
  
  fit <- glm_pp_prior$sample(
    data = standat,
    iter_warmup = 1000, iter_sampling = 2500, chains = 4, parallel_chains = 4
  )
  
  # fit$metadata()$model_params ## check the parameter names in the model
  pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]") # parameters needed for generating new data
  
  ## obtain prior samples of the parameters of interest
  d    <- fit$draws(format = 'draws_df', variables = pars) %>%
    select(all_of(pars))
  return(d)
}

cl <- makeCluster(ncores)
clusterSetRNGStream(cl, 123)
clusterExport(cl, varlist = c('formula', 'family', 'hist_data',
                              'get.stan.data.pp.prior',
                              'glm_pp_prior'))
clusterEvalQ(cl, library("dplyr"))
draws_pp <- parLapply(
  cl = cl, X = a0_list, fun = sample_pp
)
stopCluster(cl)

# save results
save(draws_pp, file = "samples_ppc/samples_pp_prior.RData")
# load("samples_ppc/samples_pp_prior.RData")