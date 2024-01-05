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
fit_mods = function(n_sim = 300, xmax = 1000, b = NULL,
                    n_iter = 2000,
                    n_chains = 2){
  sim_data = simulate_data(xmax = xmax, b = b)
  
  stan_data = make_stan_data(sim_data)
  
  fit = sampling(object = fit_model,
                 data = stan_data,
                 iter = n_iter,
                 chains = n_chains,
                 open_progress = F,
                 verbose = F)
  
  return(list("fit" = fit, "data" = sim_data, "stan_data" = stan_data))
}



# separate models --------------------------------------------------------------
# This can take ~ 2-3 hours

b = c(-2.4,
      -2.2,
      -2,
      -1.8,
      -1.6,
      -1.4,
      -1.2)

sim_vary_lambdas = replicate(mapply(fit_mods, xmax = 1000, b = b), n = 2)
# saveRDS(sim_vary_lambdas, file = "models/coverage_models/sim_vary_lambdas1000.rds")
saveRDS(sim_vary_lambdas, file = "models/sim_vary_lambdasTEMP.rds")

# summarize posteriors
sim_vary_lambdas = readRDS(file = "models/coverage_models/sim_vary_lambdas1000_4chains.rds")

indices_models <- seq(1, 21000, by = 3)  # Adjust these indices accordingly

# Use lapply to summarize each element at the specified indices
# errors for "does not contain samples
summaries_list <- lapply(indices_models, function(i) {
  summary(sim_vary_lambdas[[i]])
})       

indices_b = seq(2, 21000, by = 3)

# Use lapply to summarize each element at the specified indices
b_list <- lapply(indices_b, function(i) {
  unique(sim_vary_lambdas[[i]]$b)
})   

xmin_list <- lapply(indices_b, function(i) {
  unique(sim_vary_lambdas[[i]]$xmin)
}) 

xmax_list <- lapply(indices_b, function(i) {
  unique(sim_vary_lambdas[[i]]$xmax)
}) 

separate_summaries_temp = NULL

for(i in 1:length(summaries_list)){
  separate_summaries_temp[[i]] = summaries_list[[i]]$summary %>% 
    as_tibble() %>% 
    rownames_to_column() %>% 
    filter(rowname == 1) %>% 
    mutate(true_value = b_list[[i]],
           xmin = xmin_list[[i]],
           xmax = xmax_list[[i]],
           model = "separate_models")
}

# # repeat for model with b = -2. This was fit as part of the size range experiment
# size_range_brm1000 = readRDS("models/size_range_brm1000.rds")
# 
# separate_summaries_neg2 = NULL
# 
# for(i in 1:1000){
#   separate_summaries_neg2[[i]] = summary(size_range_brm1000[[3]][, i]$fit)$summary %>% 
#     as_tibble() %>% 
#     rownames_to_column() %>% 
#     filter(rowname == 1) %>% 
#     mutate(true_value = -2,
#            xmin = 1,
#            xmax = 1000,
#            model = "separate_models")
# }

separate_summaries = bind_rows(separate_summaries_temp) %>% 
  arrange(true_value, mean) %>% 
  mutate(replicate = rep(1:1000, 7))

saveRDS(separate_summaries, file = "posteriors/separate_summaries.rds")



# fixed model----------------------------------------------------
library(brms)
library(tidyverse)
library(isdbayes)

b = c(-2.4,
      -2.2,
      -2,
      -1.8,
      -1.6,
      -1.4,
      -1.2)

xmin = 1
xmax = 1000
n_sim = 1000

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

fit_fixed_compiled = brm(x | vreal(counts, xmin, xmax) ~ b_fac, 
                         data = temp_dat,
                         stanvars = stanvars,
                         family = paretocounts(),
                         file = "models/fit_fixed_compiled.rds",
                         chains = 1, iter = 10)


sim_fixed = function(xmin = 1, xmax = 1000, b, n_sim = 1000){
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
                         chains = 2, iter = 2000)
  
  return(list(fit_var_fixed))
  
}

fixed_mods = replicate(sim_fixed(xmax = 1000, b = b), n = 2)

# saveRDS(fixed_mods, file = "models/fixed_mods1000.rds")
saveRDS(fixed_mods, file = "models/fixed_mods_TEMP.rds")

# varying intercept model -------------------------------------------------

library(brms)
library(tidyverse)
library(isdbayes)

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

fit_varint_compiled = brm(x | vreal(counts, xmin, xmax) ~ (1|b_fac), 
                          data = temp_dat,
                          stanvars = stanvars,
                          family = paretocounts(),
                          file = "models/fit_varint_compiled.rds",
                          file_refit = "on_change",
                          prior = c(prior(normal(-1.8, 3), class = "Intercept"),
                                    prior(exponential(2), class = "sd")),
                          chains = 1, iter = 10)


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

# NOTE: this was fit twice, once with 400 sims and again with 600 sims for a total of 1000. 
# The code here just shows the 600 sims
# varint_mods = replicate(sim_varint(xmax = 1000, b = b), n = 600)

# saveRDS(varint_mods, file = "models/varint_mods600.rds")

varint_mods = replicate(sim_varint(xmax = 1000, b = b), n = 2)
saveRDS(varint_mods, file = "models/varint_mods_TEMP.rds")

# Iteration Experiment --------------------------------------------------------
# Re-run a sample of the code above with 4 chains and 2000 iterations. Ensure that it matches the original. 

b = c(-2.4,
      -2.2,
      -2,
      -1.8,
      -1.6,
      -1.4,
      -1.2)

sim_vary_lambdas_iter = replicate(mapply(fit_mods, xmax = 1000, b = b, n_chains = 4, n_iter = 2000), n = 2)
sim_vary_lambdas_smalliter = replicate(mapply(fit_mods, xmax = 1000, b = b, n_chains = 2, n_iter = 2000), n = 2)
saveRDS(sim_vary_lambdas_iter, file = "models/sim_vary_lambdas1000_iter.rds")
saveRDS(sim_vary_lambdas_smalliter, file = "models/sim_vary_lambdas1000_smalliter.rds")





