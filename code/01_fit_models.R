library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
library(sizeSpectra)
library(janitor)
rstan_options("auto_write" = TRUE)

# This script takes > 15 hours to run everything. To avoid the wait, here is each fitted model.
recover_sims_counts = readRDS(file = "posteriors/recover_sims_counts.rds")        #(simulate single samples result)
fit_varint = readRDS(file = "models/fit_varint.rds")                              #(simulate varying intercepts result)
fit_varint_regression = readRDS(file = "models/fit_varint_regression.rds")        #(simulate regression with varying intercepts result)
fit_ibts_hierarchical = readRDS(file = "models/fit_ibts_hierarchical.rds")        #(reanalyze IBTS regression with varying intercepts result)
fit_ibts_norand_norand = readRDS(file = "models/fit_ibts_norand_norand.rds")      #(reanalyze IBTS regression without varying intercepts result)
sample_size_sims = readRDS("models/sample_size_simulations/sample_size_sims.rds") #(test sample size result)

# To re-fit everything, use the code below.
set.seed(1234987)
# SIMULATE DATA AND RECOVER PARAMETERS
# simulate single samples ----------------------------------------------------------
# Simulate 10 body size datasets and fit each one with a separate truncated pareto.

# simulate
n_sim = 1000
xmax = 1000
xmin = 1
n_bs = 10
b = seq(-2.2, -1.2, length.out = n_bs)
b_true = b

sim_b = tibble(xmax = xmax, xmin = xmin,
               b = b) %>% 
  mutate(group = as.integer(1:nrow(.)))

sim_data = sim_b %>% 
  expand_grid(individual = 1:n_sim) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(b, xmin, xmax, group) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  # mutate(group = as.factor(b)) %>% 
  group_by(group) %>% 
  group_split()

# preallocate
posts_single = list()

# precompile stan code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

# fit
for (i in 1:n_bs) {
  b = unique(sim_data[[i]]$b)
  stan_dat <- list(x = sim_data[[i]]$x,
                   N = nrow(sim_data[[i]]),
                   counts = sim_data[[i]]$counts,
                   xmax = sim_data[[i]]$xmax,
                   xmin = sim_data[[i]]$xmin)
  
  fit <- sampling(object = fit_model,
                  data = stan_dat,
                  iter = 1000,
                  chains = 2,
                  open_progress = F,
                  verbose = F)
  
  posts_single[[i]] = as_draws_df(fit) %>% 
    mutate(true_value = b)
  
}

recover_sims = bind_rows(posts_single)

saveRDS(recover_sims, file = "posteriors/recover_sims_counts.rds")

# parameter_recovery_withcounts = recover_sims %>% 
#   group_by(true_value) %>% 
#   median_qi(b_exp) %>% 
#   rename(b_modeled = b_exp)
# 
# write_csv(parameter_recovery_withcounts, file = "tables/parameter_recovery_withcounts.csv")

# simulate varying intercepts ------------------------------------------------------
# Use the same body size dataset as before (10 samples) and fit them in a single hierarchical model with each dataset as a 
# varying intercept.

mod = stan_model("models/stan_spectra_mod_ibts_temperature_randonly.stan")

varint_sims = sim_b %>% 
  mutate(mean = mean(b),
         offset = b - mean,
         sd = sd(b)) # this is the known standard devation among the 10 body size data sets

saveRDS(varint_sims, file = "data/varint_sims.rds")

sim_data_varint = bind_rows(sim_data)

# convert to a list for stan
stan_dat <- list(x = sim_data_varint$x,
                 N = nrow(sim_data_varint),
                 counts = sim_data_varint$counts,
                 # mat_s = sim_data_varint$predictor,
                 n_years = length(unique(sim_data_varint$group)),
                 year = as.integer(as.factor(sim_data_varint$group)),
                 xmax = sim_data_varint$xmax,
                 xmin = sim_data_varint$xmin)

# fit model with varying intercepts
fit_varint <- sampling(object = mod,
                       data = stan_dat,
                       iter = 1000,
                       cores = 2,
                       chains = 2,
                       open_progress = F,
                       verbose = F)

saveRDS(fit_varint, file = "models/fit_varint.rds")

# # summarize
# var_int_summary = fit_varint %>% as_draws_df() %>% 
#   pivot_longer(cols = c(sigma_year, a)) %>% 
#   group_by(name) %>% 
#   median_qi(value) %>% 
#   mutate(true_value = c(intercept, var_int_sd))
# 
# write_csv(var_int_summary, file = "tables/var_int_summary.csv")

# simulate regression with varying intercepts ------------------------------------------------------
# Simulate data and fit a hierarchical model with one fixed predictor ("beta") and varying intercepts for 
# each sample. Each sample consists of 1000 individuals. The model is simulated and re-fit 40 times to estimate
# variation in parameter recovery. (It takes > 10 hours to run!!!!)

fit_varint_regression = list()

for(i in 1:40){
  intercept = -1.5
  var_int_sd = 0.1
  beta = -0.1
  n_ind = 1000
  
  lambda_sims = tibble(group = (1:3)) %>% 
    mutate(offset = rnorm(nrow(.), 0, var_int_sd),
           group_b = intercept + offset) %>% 
    expand_grid(predictor = -seq(-2, 2, length.out = 3)) %>% 
    expand_grid(replicates = 1:2) %>% 
    mutate(b = intercept + offset + predictor*beta,
           id = 1:nrow(.))
  
  sim_data = lambda_sims %>% expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = 1, xmax = 1000,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(group, id, xmin, xmax) %>% 
    add_count(x, name = "counts") %>% 
    ungroup
  
  
  # convert to a list for stan
  stan_dat <- list(x = sim_data$x,
                   N = nrow(sim_data),
                   counts = sim_data$counts,
                   mat_s = sim_data$predictor,
                   n_years = length(unique(sim_data$group)),
                   year = as.integer(as.factor(sim_data$group)),
                   xmax = sim_data$xmax,
                   xmin = sim_data$xmin)
  
  
  mod = stan_model("models/stan_spectra_mod_varint_regression.stan")
  
  # fit model with varying intercepts
  fit_varint_regression[[i]] <- sampling(object = mod,
                                         data = stan_dat,
                                         iter = 1000,
                                         # cores = 2,
                                         chains = 2,
                                         open_progress = F,
                                         verbose = F)
}


saveRDS(fit_varint_regression, file = "models/fit_varint_regression.rds")

# fit_varint_regression = readRDS(file = "models/fit_varint_regression.rds")


# REANALYZE DATA FROM EDWARDS ET AL. 2017. -----------

# reanalyze IBTS regression with varying intercepts ---------------------
data("IBTS_data")

stan_spectra_mod_ibts_temperature = stan_model("models/stan_spectra_mod_ibts_temperature.stan")

count_sims_thin = IBTS_data %>% clean_names() %>%
  group_by(year) %>% 
  mutate(xmin = min(body_mass)) %>% 
  ungroup %>% 
  mutate(xmax = max(body_mass)) %>%
  left_join(IBTS_data  %>% clean_names() %>% distinct(year) %>% 
              mutate(mat_s = (year - mean(year))/sd(year),
                     mean_year = mean(year),
                     sd_year = sd(year))) %>% 
  mutate(group = as.integer(as.factor(year))) %>% 
  group_by(year) %>% 
  sample_frac(1)

saveRDS(count_sims_thin, file = "data/count_sims_thin.rds")

stan_data_ibts = list(N = nrow(count_sims_thin),
                      mat_s = count_sims_thin$mat_s,
                      # gpp_s = count_sims_thin$gpp_s,
                      year = as.integer(as.factor(count_sims_thin$year)),
                      n_years = length(unique(count_sims_thin$year)),
                      # n_sites = length(unique(count_sims_thin$site_id_int)),
                      # site = count_sims_thin$site_id_int,
                      counts = count_sims_thin$number,
                      x = count_sims_thin$body_mass,
                      xmin = count_sims_thin$xmin,
                      xmax = count_sims_thin$xmax)

fit_ibts = sampling(object = stan_spectra_mod_ibts_temperature, 
                    data = stan_data_ibts,
                    iter = 1000, chains = 1, cores = 4)

saveRDS(fit_ibts, file = "models/fit_ibts_hierarchical.rds")

# reanalyze IBTS regression without varying intercepts ---------------------
data("IBTS_data")

stan_spectra_mod_ibts_temperature_norand = stan_model("models/stan_spectra_mod_ibts_temperature_norand.stan")

count_sims_thin = IBTS_data %>% clean_names() %>%
  group_by(year) %>% 
  mutate(xmin = min(body_mass)) %>% 
  ungroup %>% 
  mutate(xmax = max(body_mass)) %>%
  left_join(IBTS_data  %>% clean_names() %>% distinct(year) %>% 
              mutate(mat_s = (year - mean(year))/sd(year),
                     mean_year = mean(year),
                     sd_year = sd(year))) %>% 
  mutate(group = as.integer(as.factor(year))) %>% 
  group_by(year) %>% 
  sample_frac(1)

stan_data_ibts_norand = list(N = nrow(count_sims_thin),
                             mat_s = count_sims_thin$mat_s,
                             # gpp_s = count_sims_thin$gpp_s,
                             year = as.integer(as.factor(count_sims_thin$year)),
                             n_years = length(unique(count_sims_thin$year)),
                             # n_sites = length(unique(count_sims_thin$site_id_int)),
                             # site = count_sims_thin$site_id_int,
                             counts = count_sims_thin$number,
                             x = count_sims_thin$body_mass,
                             xmin = count_sims_thin$xmin,
                             xmax = count_sims_thin$xmax)

fit_ibts_norand = sampling(object = stan_spectra_mod_ibts_temperature_norand, 
                           data = stan_data_ibts_norand,
                           iter = 2000, chains = 2, cores = 4)

saveRDS(fit_ibts_norand, file = "models/fit_ibts_norand_norand.rds")

# test sample size --------------------------------------------------------
# make empty tibble 
sample_size_sims = tibble(iter = numeric(),
                          n = numeric(),
                          b = numeric(),
                          sd = numeric())
# precompile code
fit_model = stan_model("models/old/b_paretocounts_singlesample.stan")

# iterate over 1:10 replicates for each of 11 sample sizes
for (i in 1:10) {
  for (j in 2^seq(1, 11)) {
    for(k in seq(-2, -1.2, length.out = 3)) {
      n_sim = 3000
      xmax = 215000
      xmin = 0.003
      b = k
      u = runif(n_sim, min = 0, max = 1)
      
      sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
      
      
      sim_data <- tibble(x = round(sims,3)) %>% 
        sample_n(j) %>% 
        count(x, name = "counts") %>% 
        mutate(xmin = xmin,
               xmax = xmax)
      
      stan_dat <- list(x = sim_data$x,
                       N = nrow(sim_data),
                       counts = sim_data$counts,
                       xmax = sim_data$xmax,
                       xmin = sim_data$xmin)
      
      fit <- sampling(object = fit_model,
                      data = stan_dat,
                      iter = 1000,
                      chains = 2,
                      open_progress = F,
                      verbose = F)
      
      b_exp = as.data.frame(rstan::extract(fit, pars = "b_exp")) 
      
      sample_size_sims <- sample_size_sims %>%
        add_row(iter = i,
                n = j,
                b = mean(b_exp$b_exp),
                sd = sd(b_exp$b_exp))
    }
  }
}

saveRDS(sample_size_sims, "models/sample_size_simulations/sample_size_sims.rds")


