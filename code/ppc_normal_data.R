library(scoringRules)
library(mvtnorm)
library(ggplot2)
library(hdbayes)
library(dplyr)
library(MASS)
library(gridExtra)

load("samples_ppc/samples_priorpred_pp.RData")
load("samples_ppc/samples_priorpred_npp.RData")
load("samples_ppc/ess_eta.RData")
load("data/sim_lm_data.RData")

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

# # set parameters for Stan data
# formula     <- logWind ~ Ozone + Solar.R + Temp + Day
# family      <- gaussian()

formula <- y ~ X1
family      <- gaussian()

comp_crps <- function(eta){
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
  
  mean(crps_t(y, 2*alpha, mu, diag(S)))
}

a0_list      <- seq(0, 1, length.out = 40)
crps_list <- sapply(a0_list, comp_crps)
plot(a0_list, crps_list, type = 'b', xlab = expression(a[0]), ylab = 'CRPS')

best_a0_crps <- a0_list[which.min(crps_list)]

comp_hyvarinen <- function(eta) {
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
  Sinv <- solve(S)
  
  v <- y - mu
  q <- as.numeric(t(v) %*% Sinv %*% v)
  quad <- as.numeric(t(v) %*% Sinv %*% Sinv %*% v)                             # v^T S^{-2} v
  
  df <- 2*alpha
  # Hyvarinen score
  - (df + n)/(df + q) * sum(diag(Sinv)) + (df + n) * (df + n + 4) / (2*(df + q)^2) * quad
}
hyvarinen_list <- sapply(a0_list, comp_hyvarinen)
plot(a0_list, hyvarinen_list, type = 'b', xlab = expression(a[0]), ylab = 'Hyvarinen score')
best_a0_hyvarinen <- a0_list[which.min(hyvarinen_list)]

tildey_best_crps <- draws_priorpred_pp[[which.min(crps_list)]]
tildey_best_hyvarinen <- draws_priorpred_pp[[which.min(hyvarinen_list)]]

# Create legend names
name_crps <- "CRPS"
name_hyv  <- "Hyva"

# Wrap vectors as data frames
obs <- 3
df_crps <- data.frame(value = tildey_best_crps[,obs], Method = name_crps)
df_hyv  <- data.frame(value = tildey_best_hyvarinen[,obs], Method = name_hyv)
df_plot <- rbind(df_crps, df_hyv)

# Observed value for vertical line
obs_val <- curr_data$y[obs]

name_a01 <- "eta = 1"
name_a00 <- "eta = 0"

mytable <- cbind(
  Dist = c("p[eta^'*']", "p[eta^'**']", "pi", "p[1]", "p[0]"),  # remove "expression(...)" wrapper
  ESS = c(round(ess_eta[[which.min(crps_list)]], 2),
          round(ess_eta[[which.min(hyvarinen_list)]], 2),
          round(ess_npp, 2),
          round(ess_eta[[40]], 2),
          round(ess_eta[[1]], 2))
)

tg <- tableGrob(
  mytable,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(parse = TRUE, fontsize = 12)  # parse plotmath
    ),
    colhead = list(
      fg_params = list(parse = TRUE, fontsize = 12)
    )
  )
)


ggplot() +
  geom_density(data = df_plot, aes(x = value, color = Method), linewidth = 1) +
  geom_vline(aes(xintercept = obs_val, color = 'Obs value'), linewidth = 0.5) +
  labs(x = expression(tilde(y)), y = '') +
  theme_bw() +
  geom_density(data = data.frame(y = draws_priorpred_pp[[40]][,obs]), 
               aes(x = y, color = "eta = 1"), linetype = 'dashed',
               linewidth = 1.6) +
  geom_density(data = data.frame(y = draws_priorpred_pp[[1]][,obs]), 
               aes(x = y, color = "eta = 0"), linetype = 'dotted',
               linewidth = 1.6) +
  geom_density(data = data.frame(y = draws_priorpred_npp[,obs]), 
               aes(x = y, color = "NPP"),
               linewidth = 1) +
  scale_color_manual(
    name = NULL,
    # values = c(
    #   "CRPS"    = "#2EA1D1",
    #   "Hyva"    = "#B02ED1",
    #   "eta = 0" = "#4FD12E",
    #   "eta = 1" = "#D15E2E",
    #   "NPP" = "#D1C92E",
    #   "Obs value" = "black"
    # ),
    values = c(
      "CRPS"    = "#66A8D0",
      "Hyva"    = "#D06673",
      "eta = 0" = "#7f7f7f",
      "eta = 1" = "#7f7f7f",
      "NPP" = "#D0C366",
      "Obs value" = "black"
    ),
    labels = c(
      "CRPS"    = expression(p[eta^"*"]),
      "Hyva"    = expression(p[eta^"**"]),
      "eta = 0" = expression(p[0]),
      "eta = 1" = expression(p[1]),
      "NPP" = expression(pi),
      "Obs value" = expression(tilde(y)[obs])
    )
  ) +
  # geom_label(
  #   data = ess_labels,
  #   aes(x, y, label = label),
  #   inherit.aes = FALSE,
  #   fill = "white",   # box color
  #   color = "black",    # text color
  #   hjust = 0,
  #   size = 4,
  #   parse = TRUE
  # ) +
  annotation_custom(tg, xmin=-4.5, xmax = -3, ymin=0.27, ymax=0.5)+
  xlim(c(obs_val-3,obs_val+3)) +
  theme(text = element_text(size = 12),        # Base text size
        axis.title = element_text(size = 14),  # Axis titles
        axis.text = element_text(size = 12),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 11),
        legend.position = c(0.9, 0.7),
        legend.background = element_rect(
          fill = "white",     # background color of the legend
          color = "black",    # border color
          size = 0.3,         # border thickness
          linetype = "solid"  # border type
        ),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

ggsave("figures/ppc_npp_lm.png", width = 8, height = 5)
