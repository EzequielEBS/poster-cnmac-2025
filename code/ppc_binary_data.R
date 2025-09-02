# example for sampling from the NPP (prior)
library(hdbayes)
library(tidyverse)
# library(dplyr)
library(patchwork)
library(ggplot2)
library(ggridges)
library(RColorBrewer)

# load draws
load("samples_ppc/samples_npp_prior.RData")

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

## generate new data based on the prior samples and design matrix from actg019
Xhist       <- stats::model.matrix(formula, hist_data)
Xcurr       <- stats::model.matrix(formula, current_data)

## for demonstration, only generate 200 data sets
nsim <- 10000
d1.sub <- d1[sample(x = seq_len(nrow(d1)), size = nsim, replace = F), ]
d2.sub <- d2[sample(x = seq_len(nrow(d2)), size = nsim, replace = F), ]
d3.sub <- d3[sample(x = seq_len(nrow(d3)), size = nsim, replace = F), ]
d4.sub <- d4[sample(x = seq_len(nrow(d4)), size = nsim, replace = F), ]

## generate new probabilities p
p1curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d1.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(p)
}) # each column = one simulated data set
p2curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d2.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(p)
}) # each column = one simulated data set
p3curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d3.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(p)
}) # each column = one simulated data set
p4curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d4.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(p)
}) # each column = one simulated data set

p1hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d1.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(p)
}) # each column = one simulated data set
p2hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d2.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(p)
}) # each column = one simulated data set
p3hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d3.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(p)
}) # each column = one simulated data set
p4hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d4.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(p)
}) # each column = one simulated data set

## generate new outcome y
y1curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d1.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set
y2curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d2.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set
y3curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d3.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set
y4curr <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d4.sub[i, ])
  p <- binomial('logit')$linkinv(Xcurr %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set

y1hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d1.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set
y2hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d2.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set
y3hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d3.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set
y4hist <- sapply(seq_len(nsim), function(i){
  beta.sim <- as.numeric(d4.sub[i, ])
  p <- binomial('logit')$linkinv(Xhist %*% beta.sim)
  return(rbinom(length(p), 1, p))
}) # each column = one simulated data set

## compare the number of events (i.e, outcome = 1) in the simulated data sets to that in the current data
counts_sim1curr <- colSums(y1curr)
counts_sim2curr <- colSums(y2curr)
counts_sim3curr <- colSums(y3curr)
counts_sim4curr <- colSums(y4curr)
counts_obscurr <- sum(current_data$outcome)

counts_sim1hist <- colSums(y1hist)
counts_sim2hist <- colSums(y2hist)
counts_sim3hist <- colSums(y3hist)
counts_sim4hist <- colSums(y4hist)
counts_obshist <- sum(hist_data$outcome)

mean_p1curr <- colMeans(p1curr)
mean_p2curr <- colMeans(p2curr)
mean_p3curr <- colMeans(p3curr)
mean_p4curr <- colMeans(p4curr)
mean_p_obscurr <- mean(current_data$outcome)

mean_p1hist <- colMeans(p1hist)
mean_p2hist <- colMeans(p2hist)
mean_p3hist <- colMeans(p3hist)
mean_p4hist <- colMeans(p4hist)
mean_p_obshist <- mean(hist_data$outcome)

# Combine the simulated counts into a single data frame
df_countscurr <- data.frame(
  sim1 = counts_sim1curr,
  sim2 = counts_sim2curr,
  sim3 = counts_sim3curr,
  sim4 = counts_sim4curr
) %>%
  pivot_longer(cols = everything(), names_to = "prior", values_to = "counts") %>%
  mutate(prior_label = case_when(
    prior == "sim1" ~ "eta %~% Beta(1,1)",
    prior == "sim2" ~ "eta %~% Beta(2,2)",
    prior == "sim3" ~ "eta %~% Beta(1,10)",
    prior == "sim4" ~ "eta %~% Beta(10,1)"
  ))

df_countshist <- data.frame(
  sim1 = counts_sim1hist,
  sim2 = counts_sim2hist,
  sim3 = counts_sim3hist,
  sim4 = counts_sim4hist
) %>%
  pivot_longer(cols = everything(), names_to = "prior", values_to = "counts") %>%
  mutate(prior_label = case_when(
    prior == "sim1" ~ "eta %~% Beta(1,1)",
    prior == "sim2" ~ "eta %~% Beta(2,2)",
    prior == "sim3" ~ "eta %~% Beta(1,10)",
    prior == "sim4" ~ "eta %~% Beta(10,1)"
  ))

df_meanscurr <- data.frame(
  sim1 = mean_p1curr,
  sim2 = mean_p2curr,
  sim3 = mean_p3curr,
  sim4 = mean_p4curr
) %>%
  pivot_longer(cols = everything(), names_to = "prior", values_to = "mean_p") %>%
  mutate(prior_label = case_when(
    prior == "sim1" ~ "eta %~% Beta(1,1)",
    prior == "sim2" ~ "eta %~% Beta(2,2)",
    prior == "sim3" ~ "eta %~% Beta(1,10)",
    prior == "sim4" ~ "eta %~% Beta(10,1)"
  ))

df_meanshist <- data.frame(
  sim1 = mean_p1hist,
  sim2 = mean_p2hist,
  sim3 = mean_p3hist,
  sim4 = mean_p4hist
) %>%
  pivot_longer(cols = everything(), names_to = "prior", values_to = "mean_p") %>%
  mutate(prior_label = case_when(
    prior == "sim1" ~ "eta %~% Beta(1,1)",
    prior == "sim2" ~ "eta %~% Beta(2,2)",
    prior == "sim3" ~ "eta %~% Beta(1,10)",
    prior == "sim4" ~ "eta %~% Beta(10,1)"
  ))

df_meanscurr$dataset <- "Current data"
df_meanshist$dataset <- "Historical data"
df_meansall <- rbind(df_meanscurr, df_meanshist)

# Plot densities with colors
plot_y <- ggplot(df_counts, aes(x = counts, fill = prior)) +
  geom_histogram(bins = 80, color = "black", fill = "lightblue") +
  geom_vline(xintercept = counts_obs, color = "red", linetype = "dashed", 
             linewidth = 1) +
  facet_wrap(~ prior, scales = "free_y") +
  labs(
    x = "Number of events",
    y = "Count",
    title = "Prior predictive distribution of number of events"
  ) +
  theme_minimal()

plot_p1 <- ggplot(df_means, aes(x = mean_p, fill = prior)) +
  geom_density(color = "black", fill = "lightgreen") + 
  geom_vline(xintercept = mean_p_obs, color = "red", linetype = "dashed", 
             linewidth = 1) +
  facet_wrap(~ prior, scales = "free") +
  labs(
    x = "Mean predicted probability",
    y = "Density",
    title = "Prior predictive distribution of mean predicted probability"
  ) +
  theme_minimal()

plot_p2 <- ggplot(df_meanshist, aes(x = mean_p, y = prior_label, 
                                    fill = dataset, color = dataset)) +
  geom_density_ridges(alpha = 0.6, scale = 2) +
  theme_ridges() +
  
  # vertical lines (optional: same as dataset)
  geom_vline(aes(xintercept = ifelse(dataset == "Current data", mean_p_obscurr, mean_p_obshist),
                 linetype = dataset),
             color = "red", linewidth = 1) +
  
  scale_x_continuous(limits = c(0, 0.2)) +
  scale_fill_manual(values = c("Current data" = "#FFD580", 
                               "Historical data" = "skyblue"),
                    name = " ") +
  scale_color_manual(values = c("Current data" = "#FFD580", 
                               "Historical data" = "skyblue"),
                    name = " ") +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  scale_linetype_manual(
    name = "", 
    values = c("Current data" = "dashed", "Historical data" = "dotted")
  )  +
  theme_minimal() +
  theme(legend.position = "none",
        # legend.position = c(0.8, 0.85),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  ) +
  theme(text = element_text(size = 11),        # Base text size
        axis.title = element_text(size = 16),  # Axis titles
        axis.text = element_text(size = 16),   # Axis tick labels
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 16)) +
  labs(
    x = "Mean predicted probability",
    y = "",
  )
  

plot_p1
plot_p2

# save
ggsave("figures/ppc_npp_mean_prob.png", plot = plot_p2, width = 14, height = 10)
