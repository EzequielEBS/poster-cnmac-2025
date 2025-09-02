# example for sampling from the NPP (prior)
library(cmdstanr)
library(hdbayes)
library(parallel)
library(tidyverse)
library(dplyr)
library(patchwork)

# load wrapper function for getting Stan data for sampling from NPP
source("code/get_stan_data.R")


# define data
current_data <- actg036
hist_data <- actg019

# normalise data
age_stats <- with(current_data,
                  c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats <- with(current_data,
                  c('mean' = mean(cd4), 'sd' = sd(cd4)))
age_stats_hist <- with(hist_data,
                       c('mean' = mean(age), 'sd' = sd(age)))
cd4_stats_hist <- with(hist_data,
                       c('mean' = mean(cd4), 'sd' = sd(cd4)))

current_data$age <- (current_data$age - age_stats['mean']) / age_stats['sd']
current_data$cd4 <- (current_data$cd4 - cd4_stats['mean']) / cd4_stats['sd']
hist_data$age <- (hist_data$age - age_stats_hist['mean']) / age_stats_hist['sd']
hist_data$cd4 <- (hist_data$cd4 - cd4_stats_hist['mean']) / cd4_stats_hist['sd']

# define formula and family
formula     <- outcome ~ age + treatment + cd4
family      <- binomial(link = "logit")

# define stuff to compute log normalising constant
ncores  <- 14
a0      <- seq(0, 1, length.out = 21) # for demonstration, change it to a large number in practice
## wrapper to obtain log normalizing constant in parallel package

logncfun <- function(a0, ...){
  hdbayes::glm.npp.lognc(
    formula = formula, family = family, a0 = a0, histdata = hist_data,
    ...
  )
}

cl <- makeCluster(ncores)
clusterSetRNGStream(cl, 123)
clusterExport(cl, varlist = c('formula', 'family', 'hist_data'))
a0.lognc <- parLapply(
  cl = cl, X = a0, fun = logncfun, iter_warmup = 1000,
  iter_sampling = 2500, chains = 4
)
stopCluster(cl)
a0.lognc <- data.frame( do.call(rbind, a0.lognc) )

## get Stan data for sampling from NPP
## similar input as for hdbayes::glm.npp()
standat1 <- get.stan.data.npp.prior(
  formula        = formula,
  family         = family,
  data.list      = list(hist_data),
  a0.lognc       = a0.lognc$a0,
  lognc          =  matrix(a0.lognc$lognc, ncol = 1),
  offset.list    = NULL,
  beta.mean      = 0,
  beta.sd        = 10,
  disp.mean      = NULL,
  disp.sd        = NULL,
  a0.shape1      = 1,
  a0.shape2      = 1,
  a0.lower       = 0,
  a0.upper       = 1
)
standat2 <- get.stan.data.npp.prior(
  formula        = formula,
  family         = family,
  data.list      = list(hist_data),
  a0.lognc       = a0.lognc$a0,
  lognc          =  matrix(a0.lognc$lognc, ncol = 1),
  offset.list    = NULL,
  beta.mean      = 0,
  beta.sd        = 10,
  disp.mean      = NULL,
  disp.sd        = NULL,
  a0.shape1      = 2,
  a0.shape2      = 2,
  a0.lower       = 0,
  a0.upper       = 1
)
standat3 <- get.stan.data.npp.prior(
  formula        = formula,
  family         = family,
  data.list      = list(hist_data),
  a0.lognc       = a0.lognc$a0,
  lognc          =  matrix(a0.lognc$lognc, ncol = 1),
  offset.list    = NULL,
  beta.mean      = 0,
  beta.sd        = 10,
  disp.mean      = NULL,
  disp.sd        = NULL,
  a0.shape1      = 1,
  a0.shape2      = 10,
  a0.lower       = 0,
  a0.upper       = 1
)
standat4 <- get.stan.data.npp.prior(
  formula        = formula,
  family         = family,
  data.list      = list(hist_data),
  a0.lognc       = a0.lognc$a0,
  lognc          =  matrix(a0.lognc$lognc, ncol = 1),
  offset.list    = NULL,
  beta.mean      = 0,
  beta.sd        = 10,
  disp.mean      = NULL,
  disp.sd        = NULL,
  a0.shape1      = 10,
  a0.shape2      = 1,
  a0.lower       = 0,
  a0.upper       = 1
)

## fit model in cmdstanr
fit1 <- glm_npp_prior$sample(
  data = standat1,
  iter_warmup = 1000, iter_sampling = 2500, chains = 4, parallel_chains = 4
)
fit2 <- glm_npp_prior$sample(
  data = standat2,
  iter_warmup = 1000, iter_sampling = 2500, chains = 4, parallel_chains = 4
)
fit3 <- glm_npp_prior$sample(
  data = standat3,
  iter_warmup = 1000, iter_sampling = 2500, chains = 4, parallel_chains = 4
)
fit4 <- glm_npp_prior$sample(
  data = standat4,
  iter_warmup = 1000, iter_sampling = 2500, chains = 4, parallel_chains = 4
)

# fit$metadata()$model_params ## check the parameter names in the model
pars <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]") # parameters needed for generating new data

## obtain prior samples of the parameters of interest
d1    <- fit1$draws(format = 'draws_df', variables = pars) %>%
  select(all_of(pars))
d2    <- fit2$draws(format = 'draws_df', variables = pars) %>%
  select(all_of(pars))
d3    <- fit3$draws(format = 'draws_df', variables = pars) %>%
  select(all_of(pars))
d4    <- fit4$draws(format = 'draws_df', variables = pars) %>%
  select(all_of(pars))

# save the draws
save(d1, d2, d3, d4, file = "samples_ppc/samples_npp_prior.RData")
# load("samples_ppc/samples_npp_prior.RData")