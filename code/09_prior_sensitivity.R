library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
rstan_options("auto_write" = TRUE)

set.seed(1234987)


# precompile modesl
fit_model = stan_model("models/b_paretocounts_singlesample.stan")
fit_model_strong_prior = stan_model("models/b_paretocounts_singlesample_strong_prior.stan")
fit_model_strongest_prior = stan_model("models/b_paretocounts_singlesample_strongest_prior.stan")

# load data created 03_recover_parameters.R
sim_data = readRDS(file = "data/sim_data.rds")

# limit to a single sample
stan_dat <- list(x = sim_data[[3]]$x,
                 N = nrow(sim_data[[3]]),
                 counts = sim_data[[3]]$counts,
                 xmax = sim_data[[3]]$xmax,
                 xmin = sim_data[[3]]$xmin,
                 true_lambda = sim_data[[3]]$lambda)


fit_reg_prior <- sampling(object = fit_model,
                          data = stan_dat,
                          iter = 1000,
                          chains = 2,
                          open_progress = F,
                          verbose = F)

fit_strong_prior <- sampling(object = fit_model_strong_prior,
                          data = stan_dat,
                          iter = 1000,
                          chains = 2,
                          open_progress = F,
                          verbose = F)

fit_strongest_prior <- sampling(object = fit_model_strongest_prior,
                          data = stan_dat,
                          iter = 1000,
                          chains = 2,
                          open_progress = F,
                          verbose = F)

saveRDS(fit_reg_prior, file = "models/fit_reg_prior.rds")
saveRDS(fit_strong_prior, file = "models/fit_strong_prior.rds")
saveRDS(fit_strongest_prior, file = "models/fit_stronest_prior.rds")


# summarize fits ----------------------------------------------------------

fit_reg_prior = readRDS(file = "models/fit_reg_prior.rds")
fit_strong_prior = readRDS(file = "models/fit_strong_prior.rds")
fit_strongest_prior = readRDS(file = "models/fit_stronest_prior.rds")


fit_reg_post = as_draws_df(fit_reg_prior) %>% mutate(prior = "N(-1.8, 2)")
fit_strong_post = as_draws_df(fit_strong_prior) %>% mutate(prior = "N(-1.8, 0.1)")
fit_strongest_post = as_draws_df(fit_strongest_prior) %>% mutate(prior = "N(-1.8, 0.01)")

fit_all_post = bind_rows(fit_reg_post,
                          fit_strong_post,
                          fit_strongest_post)


prior_sensitivity = fit_all_post %>% 
  group_by(prior) %>% 
  summarize(mean = mean(lambda),
            sd = sd(lambda),
            lower = quantile(lambda, probs = 0.025),
            upper = quantile(lambda, probs = 0.975)) %>% 
  arrange(mean)


write_csv(prior_sensitivity, file = "tables/prior_sensitivity.csv")













