library(tidyverse)
library(rstan)

# rstan_options("auto_write" = TRUE)

# To re-fit everything, use the code below.
# SIMULATE DATA AND RECOVER PARAMETERS
# simulate single samples ----------------------------------------------------------
# Simulate 10 body size datasets and fit each one with a separate truncated pareto.

# precompile Stan code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")


# functions ---------------------------------------------------------------

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
                 chains = 2,
                 open_progress = F,
                 verbose = F)
  
  return(list("fit" = fit, "data" = sim_data, "stan_data" = stan_data))
}



# fit models --------------------------------------------------------------
# This can take ~ 2-3 hours

b = c(-2.4,
      -2.2,
      -2,
      -1.8,
      -1.6,
      -1.4,
      -1.2)

sim_vary_lambdas = replicate(mapply(fit_mods, xmax = 1000, b = b), n = 1000)
saveRDS(sim_vary_lambdas, file = "models/sim_vary_lambdas1000.rds")





















