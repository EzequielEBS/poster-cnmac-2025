library(mvtnorm)

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

sample_priorpred <- function(a0) {
  res          = hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
  y0            = res$y
  X0            = res$X
  alpha <- 2
  beta <- 1
  sigma0 <- 1 / a0 * solve(t(X0) %*% X0)
  mu <- a0 * X0 %*% sigma0 %*% t(X0) %*% y0
  sigma <- beta / alpha * (diag(1, nrow(X0)) + X0 %*% sigma0 %*% t(X0))
  
  n <- 10000
  tildey <- rmvt(n, delta = mu, sigma = sigma, df = 2 * alpha)
  return(tildey)
}

ncores  <- 14
a0_list      <- seq(0.05, 1, length.out = 40)

cl <- makeCluster(ncores)
clusterSetRNGStream(cl, 123)
clusterExport(cl, varlist = c('formula', 'family', 'hist_data'))
clusterEvalQ(cl, {library("dplyr") 
                 library("mvtnorm")})
draws_priorpred <- parLapply(
  cl = cl, X = a0_list, fun = sample_priorpred
)
stopCluster(cl)

save(draws_priorpred, file = 'samples_ppc/samples_priorpred.RData')


hist(draws_priorpred[[1]][,1])
abline(v = hist_data$logWind[1], col = 'red')
