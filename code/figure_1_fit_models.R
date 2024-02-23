library(tidyverse)
library(rstan)
library(isdbayes)
library(brms)

n_replicates = 2 # To repeat the analysis in the paper, set n_replicates = 1000 (But this will take 5-6 days) 

# Figure 1a ----------------------------------------------------------
# Simulate 1000 body size datasets and fit each one with a separate truncated pareto.

# precompile Stan code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

# functions 

# function to simulate isd data
simulate_data = function(n_sim = 300, xmax = 1000, xmin = 1, b = -1.6){
  
  tibble(xmin = xmin) %>% 
    expand_grid(xmax, b, n_sim) %>% 
    group_by(xmax, b) %>% 
    uncount(n_sim, .remove = F) %>%
    ungroup %>% 
    mutate(u = runif(nrow(.), min = 0, max = 1),
           sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(b, xmin, xmax, n_sim) %>% 
    mutate(x = round(sims, 3)) %>% 
    count(x, name = "counts") %>% 
    group_by(n_sim, b, xmax) %>%
    mutate(id = cur_group_id()) 
}

# function to convert data to stan
make_stan_data = function(sim_data){
  list(x = sim_data$x,
       N = nrow(sim_data),
       counts = sim_data$counts,
       xmax = sim_data$xmax,
       xmin = sim_data$xmin,
       nsims = nrow(sim_data))
}

# function to fit data to models
fit_mods = function(n_sim = 300, xmax = 1000, b = NULL){
  sim_data = simulate_data(xmax = xmax, b = b)
  
  stan_data = make_stan_data(sim_data)
  
  fit = sampling(object = fit_model,
                 data = stan_data,
                 iter = 2000,
                 chains = 4,
                 open_progress = F,
                 verbose = F)
  
  return(list("fit" = fit, "data" = sim_data, "stan_data" = stan_data))
}



# fit models
# This can take ~ 2-3 hours

b = c(-2.4,
      -2.2,
      -2,
      -1.8,
      -1.6,
      -1.4,
      -1.2)

sim_vary_lambdas = replicate(mapply(fit_mods, xmax = 1000, b = b), n = n_replicates)
saveRDS(sim_vary_lambdas, file = "models/fig1a_mods.rds")


# Figure 1b ----------------------------------------------------
library(brms)
library(tidyverse)
library(isdbayes)

#1) set parameters
# lambdas
b = c(-2.4,
      -2.2,
      -2,
      -1.8,
      -1.6,
      -1.4,
      -1.2
)
xmin = 1
xmax = 1000
n_sim = 300 # number of individuals

# 2) make temporary data
temp_dat = tibble(xmin = xmin) %>% 
  expand_grid(xmax, b, n_sim) %>% 
  group_by(xmax, b) %>% 
  uncount(n_sim, .remove = F) %>%
  ungroup %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(b, xmin, xmax, n_sim) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  group_by(n_sim, b, xmax) %>%
  mutate(id = cur_group_id(),
         b_fac = as.factor(b)) %>% 
  ungroup() 

# 3) fit a toy model. We will update this with the proper chains in the loop below
fit_fixed_compiled = brm(x | vreal(counts, xmin, xmax) ~ b_fac, 
                         data = temp_dat,
                         stanvars = stanvars,
                         family = paretocounts(),
                         file = "models/fit_fixed_compiled.rds",
                         chains = 1, iter = 10)

# 4) Function to simulate data and fit model
sim_fixed = function(xmin = 1, xmax = 1000, b, n_sim = 300, chains = 1, iter = 20) {
  var_int_dat = tibble(xmin = xmin) %>% 
    expand_grid(xmax, b, n_sim) %>% 
    group_by(xmax, b) %>% 
    uncount(n_sim, .remove = F) %>%
    ungroup %>% 
    mutate(u = runif(nrow(.), min = 0, max = 1),
           sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(b, xmin, xmax, n_sim) %>% 
    mutate(x = round(sims, 3)) %>% 
    count(x, name = "counts") %>% 
    group_by(n_sim, b, xmax) %>%
    mutate(id = cur_group_id(),
           b_fac = as.factor(b)) %>% 
    ungroup() 
  
  fit_var_fixed = update(fit_fixed_compiled, newdata = var_int_dat,
                         chains = chains, iter = iter)
  
  return(list(fit_var_fixed))
}

# 5) Fit models (This is silenced. It takes multiple days)
fig1b_mods <- replicate(n = n_replicates,
                        sim_fixed(xmax = 1000, b = b, chains = 2, iter = 2000),
                        simplify = FALSE)


# saveRDS(fig1b_mods, file = "models/fig1b_mods.rds")

# Figure 1c ----------------------------------------------------
library(brms)
library(tidyverse)
library(isdbayes)

#1) Set lambdas and parameters
b = c(-2.4,
      -2.2,
      -2,
      -1.8,
      -1.6,
      -1.4,
      -1.2)
xmin = 1
xmax = 1000
n_sim = 300

# 2) Make temporary data
temp_dat = tibble(xmin = xmin) %>% 
  expand_grid(xmax, b, n_sim) %>% 
  group_by(xmax, b) %>% 
  uncount(n_sim, .remove = F) %>%
  ungroup %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(b, xmin, xmax, n_sim) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  group_by(n_sim, b, xmax) %>%
  mutate(id = cur_group_id(),
         b_fac = as.factor(b)) %>% 
  ungroup() 

# 3) Fit toy model with low iterations (to update next)
fit_varint_compiled = brm(x | vreal(counts, xmin, xmax) ~ (1|b_fac), 
                          data = temp_dat,
                          stanvars = stanvars,
                          family = paretocounts(),
                          file = "models/fit_varint_compiled.rds",
                          file_refit = "on_change",
                          prior = c(prior(normal(-1.8, 3), class = "Intercept"),
                                    prior(exponential(2), class = "sd")),
                          chains = 1, iter = 10)

# 4) Function to simulate data and fit models
sim_varint = function(xmin = 1, xmax = 1000, b, n_sim = 300){
  var_int_dat = tibble(xmin = xmin) %>% 
    expand_grid(xmax, b, n_sim) %>% 
    group_by(xmax, b) %>% 
    uncount(n_sim, .remove = F) %>%
    ungroup %>% 
    mutate(u = runif(nrow(.), min = 0, max = 1),
           sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(b, xmin, xmax, n_sim) %>% 
    mutate(x = round(sims, 3)) %>% 
    count(x, name = "counts") %>% 
    group_by(n_sim, b, xmax) %>%
    mutate(id = cur_group_id(),
           b_fac = as.factor(b)) %>% 
    ungroup() 
  
  fit_var_varint = update(fit_varint_compiled, newdata = var_int_dat,
                          chains = 2, iter = 2000)
  
  return(list(fit_var_varint))
  
}

# 5) Fit models (Silenced due to time limitations)
varint_mods = replicate(n = n_replicates, sim_varint(xmax = 1000, b = b))

# saveRDS(varint_mods, file = "models/fig1c_mods.rds")






