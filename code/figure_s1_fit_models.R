library(tidyverse)
library(rstan)
library(isdbayes)
library(brms)

# 1) Simulate data
# simulate stream fish and inverts using approximate size ranges of fish from
# the National Ecological Observatory Network

sim_dat = function(mu = -1.2, n = 1000, fish_xmin = 1000, fish_xmax = 2e05,
                   invert_xmin = 0.003, invert_xmax = 30600){
  fish_sims = tibble(x = rparetocounts(vreal2 = fish_xmin, vreal3 = fish_xmax, mu = mu, 
                                       n = n)) %>% 
    mutate(sample_area_m2 = 500, 
           counts = 1,
           xmin = fish_xmin, 
           xmax = fish_xmax,
           mu = mu,
           taxon = "fish, counts = 1")
  
  invert_sims = tibble(x = rparetocounts(vreal2 = invert_xmin, vreal3 = invert_xmax, mu = mu,
                                         n = n)) %>% 
    mutate(sample_area_m2 = 0.09,
           counts = 1,
           xmin = invert_xmin,
           xmax = invert_xmax,
           mu = mu,
           taxon = "inverts, counts = 1") 
  
  fish_invert_sims = bind_rows(fish_sims, invert_sims) %>% 
    mutate(xmin = invert_xmin, 
           xmax = fish_xmax,
           mu = mu,
           taxon = "inverts + fish, counts = 1")
  
  fish_invert_sims_fixcounts = bind_rows(fish_sims, invert_sims) %>% 
    mutate(xmin = invert_xmin, 
           xmax = fish_xmax,
           counts = 1/sample_area_m2,
           mu = mu,
           taxon = "inverts + fish, counts = 1/sample_area")
  
  dat_sims = list(fish_sims, invert_sims, fish_invert_sims, fish_invert_sims_fixcounts)
  
}


# 2) fit dummy model
temp_mod = brm(x|vreal(counts, xmin ,xmax) ~ 1,
               data = sim_dat(mu = -1.2)[1],
               stanvars = stanvars,    # required for truncated Pareto
               family = paretocounts(),
               chains = 1, iter = 20) # required for truncated Pareto

# 3) update dummy model with data. Replicate each simulation and model fit 100 times for different lambdas
fish_mod_1.2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.2)[1], chains = 2, iter = 1000), simplify = F)
invert_mod_1.2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.2)[2], chains = 2, iter = 1000), simplify = F)
fish_invert_mod_1.2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.2)[3], chains = 2, iter = 1000), simplify = F)
fish_invert_mod_fixcounts_1.2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.2)[4], chains = 2, iter = 1000), simplify = F)

fish_invert_mods_1.2 = list(fish_mod_1.2,
                            invert_mod_1.2,
                            fish_invert_mod_1.2,
                            fish_invert_mod_fixcounts_1.2)


# saveRDS(fish_invert_mods_1.2, file = "models/figs1a_mods.rds")


fish_mod_1.6 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.6)[1], chains = 2, iter = 1000), simplify = F)
invert_mod_1.6 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.6)[2], chains = 2, iter = 1000), simplify = F)
fish_invert_mod_1.6 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.6)[3], chains = 2, iter = 1000), simplify = F)
fish_invert_mod_fixcounts_1.6 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -1.6)[4], chains = 2, iter = 1000), simplify = F)

fish_invert_mods_1.6 = list(fish_mod_1.6,
                            invert_mod_1.6,
                            fish_invert_mod_1.6,
                            fish_invert_mod_fixcounts_1.6)


# saveRDS(fish_invert_mods_1.6, file = "models/figs1b_mods.rds")


fish_mod_2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -2)[1], chains = 2, iter = 1000), simplify = F)
invert_mod_2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -2)[2], chains = 2, iter = 1000), simplify = F)
fish_invert_mod_2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -2)[3], chains = 2, iter = 1000), simplify = F)
fish_invert_mod_fixcounts_2 = replicate(100, update(temp_mod, newdata = sim_dat(mu = -2)[4], chains = 2, iter = 1000), simplify = F)


fish_invert_mods_2 = list(fish_mod_2,
                          invert_mod_2,
                          fish_invert_mod_2,
                          fish_invert_mod_fixcounts_2)


# saveRDS(fish_invert_mods_2, file = "models/figs1c_mods.rds")

