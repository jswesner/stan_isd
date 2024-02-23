
library(tidyverse)
library(rstan)
library(isdbayes)
library(brms)

# simulate stream fish and inverts using approximate size ranges of fish from
# the National Ecological Observatory Network

sim_dat = function(mu = -1.2, n = 1000, fish_xmin = 1000, fish_xmax = 2e05,
                   invert_xmin = 0.003, invert_xmax = 30600){
  fish_sims = tibble(x = rparetocounts(xmin = fish_xmin, xmax = fish_xmax, lambda = mu, 
                                       n = n)) %>% 
    mutate(sample_area_m2 = 500, 
           counts = 1,
           xmin = fish_xmin, 
           xmax = fish_xmax,
           mu = mu,
           taxon = "fish, counts = 1")
  
  invert_sims = tibble(x = rparetocounts(xmin = invert_xmin, xmax = invert_xmax, lambda = mu,
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