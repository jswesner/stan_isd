library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
library(sizeSpectra)
library(janitor)
library(isdbayes)
rstan_options("auto_write" = TRUE)

theme_set(theme_default())

# 1) set data simulation parameters
set.seed = 123123
intercept = -1.5
beta = -0.1
n_ind = 300
reps = 12

# 2) make grid of lambdas
lambda_sims = tibble(group = (1:1)) %>% 
    expand_grid(predictor = -seq(-2, 2, length.out = 12)) %>% 
    # expand_grid(replicates = 1:2) %>% 
    mutate(b = intercept + predictor*beta,
           id = 1:nrow(.)) %>% 
  add_row(group = 1,
          predictor = 2.5, b = -1.1, id = 7) 
  
saveRDS(lambda_sims, file = "data/lambda_sims.rds") # save for plotting later

# 3) simulate data
sim_data = lambda_sims %>% 
  expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = 1, xmax = 1000,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(group, id, xmin, xmax, predictor) %>% 
    add_count(x, name = "counts") %>% 
  mutate(remove = case_when(predictor == 2.5 & individual > 50 ~ "remove",
                            TRUE ~ "keep")) %>% 
  filter(remove == "keep") %>%
    ungroup %>% 
  mutate(group = as.factor(predictor))


# 4) fit two-step models -------------------------------------
# Step 1: estimate lambdas separately

# split data into a list (one for each lambda)
sim_data_list = sim_data %>% group_by(group) %>% group_split()

# fit and save separate intercept only models to each lambda
brm_single_mods = brm_multiple(x | vreal(counts, xmin, xmax) ~ 1,
                               data = sim_data_list,
                               stanvars = stanvars,
                               family = paretocounts(),
                               chains = 4, iter = 2000,
                               cores = 4,
                               combine = F)

saveRDS(brm_single_mods, file = "models/brm_single_mods.rds")
brm_single_mods = readRDS("models/brm_single_mods.rds")

# extract, wrangle, and save posterior lambdas from each separate model
brm_single_fixefs = NULL

for(i in 1:length(brm_single_mods)){
  brm_single_fixefs[[i]] = as_tibble(fixef(brm_single_mods[[i]])) %>% 
    mutate(group = i, 
           reps = nrow(brm_single_mods[[i]]$data))
}

single_lambdas = bind_rows(brm_single_fixefs) %>% 
  mutate(predictor = sort(unique(sim_data$predictor)))

saveRDS(single_lambdas, file = "data/single_lambdas.rds")

# Step 2: Fit Gaussian regression between lambdas and predictor
# simple linear regression
brm_twostep_regress = brm(Estimate ~ predictor,
                          data = single_lambdas,
                          family = gaussian(),
                          iter = 2000,
                          chains = 4)

# saveRDS(brm_twostep_regress, file = "models/fig4a_mod.rds")

# weighted regression
brm_twostep_regress_mi = brm(Estimate|mi(Est.Error) ~ predictor,
                             data = single_lambdas,
                             iter = 2000,
                             chains = 4)

# saveRDS(brm_twostep_regress, file = "models/fig4b_mod.rds")


# 5) fit hierarchical models -------------------------------------------------

# fig4c: weak prior
brm_regress_weakprior_rand = brm(x | vreal(counts, xmin, xmax) ~ predictor + (1|group),
                                 data = sim_data,
                                 stanvars = stanvars,
                                 family = paretocounts(),
                                 prior = c(prior(normal(0, 0.5), class = "b"),
                                           prior(normal(-1.5, 1), class = "Intercept"),
                                           prior(exponential(1), class = "sd")),
                                 iter = 2000,
                                 chains = 4)

# saveRDS(brm_regress_weakprior_rand, file = "models/fig4c_mod.rds")

# fig4d: strong prior
brm_regress_strongprior_rand = brm(x | vreal(counts, xmin, xmax) ~ predictor + (1|group),
                             data = sim_data,
                             stanvars = stanvars,
                             family = paretocounts(),
                             prior = c(prior(normal(-0.1, 0.02), class = "b"),
                                       prior(normal(-1.5, 0.1), class = "Intercept"),
                                       prior(exponential(1), class = "sd")),
                             iter = 2000,
                             chains = 4)

# saveRDS(brm_regress_strongprior_rand, file = "models/fig4d_mod.rds")
