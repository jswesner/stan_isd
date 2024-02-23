library(tidyverse)
library(rstan)
library(isdbayes)
library(brms)

source("code/functions/sim_dat.R")

# 1) Simulate data
set.seed(1111222)
dat = sim_dat(mu = -1.6)[[4]] %>%
  mutate(counts = sample_area_m2*counts,
         sample_area_10m2 = sample_area_m2/10,
         sample_area_0.1m2 = sample_area_m2*10) %>%
  pivot_longer(cols = starts_with("sample_area")) %>%
  mutate(counts = counts/value)

# 2) Load precompiled intercept only model
temp_mod = readRDS(file = "models/precompiled/intercept_only.rds")

# 3) Fit models with different sampling areas
mod_0.1m2 = update(temp_mod, newdata = dat %>% filter(name == "sample_area_0.1m2"),
                 iter = 2000, chains = 4)

mod_m2 = update(temp_mod, newdata = dat %>% filter(name == "sample_area_m2"),
                 iter = 2000, chains = 4)

mod_10m2 = update(temp_mod, newdata = dat %>% filter(name == "sample_area_10m2"),
                 iter = 2000, chains = 4)

# 4) Save models
mods_area = list(mod_0.1m2, mod_m2, mod_10m2)
saveRDS(mods_area, file = "models/figs2_mods.rds")



