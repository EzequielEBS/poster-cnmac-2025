library(scoringRules)
library(mvtnorm)
library(ggplot2)

load("samples_ppc/samples_priorpred.RData")


comp_crps <- function(a0){
  res          = hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
  y0            = res$y
  X0            = res$X
  alpha <- 2
  beta <- 1
  sigma0 <- 1 / a0 * solve(t(X0) %*% X0)
  mu <- a0 * X0 %*% sigma0 %*% t(X0) %*% y0
  sigma <- beta / alpha * (diag(1, nrow(X0)) + X0 %*% sigma0 %*% t(X0))
  
  mean(crps_t(y0, 2*alpha, mu, diag(sigma)))
}

crps_list <- sapply(a0_list, comp_crps)
plot(a0_list, crps_list, type = 'b', xlab = expression(a[0]), ylab = 'CRPS')

best_a0_crps <- a0_list[which.min(crps_list)]

comp_hyvarinen <- function(a0) {
  res          = hdbayes:::stack.data(formula = formula, data.list = list(hist_data))
  y0            = res$y
  X0            = res$X
  n <- length(y0)
  alpha <- 2
  beta <- 1
  sigma0 <- 1 / a0 * solve(t(X0) %*% X0)
  mu <- a0 * X0 %*% sigma0 %*% t(X0) %*% y0
  sigma <- beta / alpha * (diag(1, nrow(X0)) + X0 %*% sigma0 %*% t(X0))
  Sinv <- solve(sigma)
  q <- as.numeric(t(y0 - mu) %*% Sinv %*% (y0 - mu))
  df <- 2*alpha
  quad  <- as.numeric(t(y0 - mu) %*% (Sinv %*% (Sinv %*% (y0 - mu))))
  - (df + n) /(df + q) * sum(diag(Sinv)) + (df + n) * (df + n + 4) / 
    (2*(df + q)^2) * quad
}
hyvarinen_list <- sapply(a0_list, comp_hyvarinen)
plot(a0_list, hyvarinen_list, type = 'b', xlab = expression(a[0]), ylab = 'Hyvarinen score')
best_a0_hyvarinen <- a0_list[which.min(hyvarinen_list)]

tildey_best_crps <- draws_priorpred[[which.min(crps_list)]]
tildey_best_hyvarinen <- draws_priorpred[[which.min(hyvarinen_list)]]

# Create legend names
name_crps <- paste0('Best a0 (CRPS): ', round(best_a0_crps, 2))
name_hyv  <- paste0('Best a0 (Hyvarinen): ', round(best_a0_hyvarinen, 2))

# Wrap vectors as data frames
df_crps <- data.frame(value = tildey_best_crps[,1], Method = name_crps)
df_hyv  <- data.frame(value = tildey_best_hyvarinen[,1], Method = name_hyv)
df_plot <- rbind(df_crps, df_hyv)

# Observed value for vertical line
obs_val <- hist_data$logWind[1]

ggplot(df_plot, aes(x = value, color = Method)) +
  geom_density(linetype = 'dashed', size = 1) +
  geom_vline(xintercept = obs_val, color = 'black', size = 1) +
  labs(x = expression(tilde(y)[1]), y = 'Density') +
  theme_minimal() +
  theme(text = element_text(size = 20), legend.position = 'top') +
  ggtitle('Prior predictive distributions for observation 30') +
  scale_color_manual(name = "", 
                     values = c("Best a0 (CRPS): 1" = 'blue', "Best a0 (Hyvarinen): 1" = 'red'))


