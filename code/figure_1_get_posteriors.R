library(tidyverse)
library(rstan)
library(tidybayes)
library(isdbayes)
library(brms)


# Figure 1a - separate models ---------------------------------------------------------

# summarize posteriors
sim_vary_lambdas = readRDS(file = "models/fig1a_mods.rds")

# get only the stanfits, which are stored in every three slots
sim_vary_lambdas_fits = sim_vary_lambdas[sapply(sim_vary_lambdas, function(x) is(x, "stanfit"))]

# check that there are no missing fits
which(sapply(sim_vary_lambdas_fits, is.null))

# extract posterior summaries
separate_summaries_temp = NULL

for(i in 1:length(sim_vary_lambdas_fits)){
   separate_summaries_temp[[i]] =  median_qi(rstan::extract(sim_vary_lambdas_fits[[i]])$lambda) %>% 
    mutate(iter = sim_vary_lambdas_fits[[i]]@sim$iter,
           chains = sim_vary_lambdas_fits[[i]]@sim$chains,
           model = "separate",
           model_number = i)
}

# combine posterior summaries and save
separate_summaries = bind_rows(separate_summaries_temp) %>% 
  mutate(replicate = rep(1:1000, 7),
         true_value = rep(c(-2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2), 1000)) %>% 
  rename('50%' = y,
         '2.5%' = ymin,
         '97.5%' = ymax) %>% 
  mutate(par = "lambda") %>% 
  ungroup() %>% 
  select(-.width, -.point, -.interval) %>% 
  mutate(true_value = as.factor(true_value))

# saveRDS(separate_summaries, file = "posteriors/fig1a_posterior_summaries.rds")


# Figure 1b - fixed models----------------------------------------------------

# 1) load fitted models
fixed_mods = readRDS("models/fig1b_mods.rds")

# 2) extract posterior summaries using a for loop
fixed_lambda_summaries_list = NULL

for(i in 1:length(fixed_mods)){
    fixed_lambda_summaries_list[[i]] =  fitted(fixed_mods[[1]][[1]],
                                          newdata = fixed_mods[[1]][[1]]$data %>% distinct(b_fac) %>% 
                                            mutate(xmin = 1, xmax = 1000, counts = 1)) %>% 
      as_tibble() %>% 
      mutate(true_value = unique(fixed_mods[[1]][[1]]$data$b_fac)) %>%
      mutate(iter = summary(fixed_mods[[i]][[1]])$iter,
             chains = summary(fixed_mods[[i]][[1]])$chains,
             replicate = i) %>%
      mutate(model = "fixed",
             `2.5%` = Q2.5,
             `50%` = Estimate,
             `97.5%` = Q97.5) %>% 
      select(-Est.Error, -Estimate, -Q2.5, -Q97.5)
  }

# 3) bind summaries
fixed_lambda_summaries = bind_rows(fixed_lambda_summaries_list)

# 4) save summaries
# saveRDS(fixed_lambda_summaries, file = "posteriors/fig1b_posterior_summaries.rds")


# Figure 1c - varying intercept models -------------------------------------------------
varint_mod = readRDS("models/fig1c_mods.rds")

posts_varint_list = NULL

for(i in 1:length(varint_mod)){
  posts_varint_list[[i]] = varint_mod[[i]]$data %>% 
    distinct(b_fac, xmin, xmax) %>% 
    mutate(counts = 1) %>% 
    add_epred_draws(varint_mod[[i]]) %>% 
    group_by(xmin, xmax, b_fac) %>% 
    reframe(mean = mean(.epred), 
            se_mean = sd(.epred),
            '0.1%' = quantile(.epred, probs = 0.001),
            '2.5%' = quantile(.epred, probs = 0.025),
            '25%' = quantile(.epred, probs = 0.25),
            '50%' = quantile(.epred, probs = 0.5),
            '75%' = quantile(.epred, probs = 0.75),
            '97.5%' = quantile(.epred, probs = 0.975),
            '99.9%' = quantile(.epred, probs = 0.999)) %>% 
    rename(true_value = b_fac) %>% 
    mutate(model = "varying_intercepts",
           chains = 2, 
           iter = 2000)
}

var_lambda_summaries = bind_rows(posts_varint_list) %>% 
  mutate(replicate = rep(1:1000, each = 7)) 

# saveRDS(var_lambda_summaries, file = "posteriors/fig1c_posterior_summaries.rds")

# Combine all posteriors --------------------------------------------------
separate_lambda_summaries = readRDS(file = "posteriors/fig1a_posterior_summaries.rds")
fixed_lambda_summaries = readRDS(file = "posteriors/fig1b_posterior_summaries.rds")
var_lambda_summaries = readRDS(file = "posteriors/fig1c_posterior_summaries.rds")


fig1abc_posterior_summaries = bind_rows(fixed_lambda_summaries,
                           var_lambda_summaries,
                           separate_lambda_summaries) %>% 
  as_tibble()

# saveRDS(fig1abc_posterior_summaries.rds, file = "posteriors/fig1abc_posterior_summaries.rds")

