library(tidyverse)
library(rstan)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.

#1) Get fitted models
# This loads 1000-7000 fitted models. It is too big for GitHub (500-1000MB), so we have to load it from Zenodo.
# It might take 5-10 minutes
options(timeout=1000)
# reg_update = readRDS("models/coverage_models/reg_fits_fixed.rds")
reg_update = readRDS(url("https://zenodo.org/records/10689765/files/fig3_mods.rds?download=1"))

#2) Extract mean and CrI from fitted models
reg_summaries = NULL

for(i in 1:ncol(reg_update)){
  reg_summaries[[i]] = as_draws_df(reg_update[,i]$fit) %>%
    pivot_longer(starts_with("b_")) %>% 
    group_by(name) %>% 
    tidybayes::median_qi(value) %>% 
    mutate(replicate = i)
}

#3) Get true values to append
true_values = tibble(b_Intercept = -1.2,
                     b_predictor = -0.05) %>% 
  pivot_longer(cols = everything(), names_to = "name", values_to = "true_values") %>% 
  mutate(name_greek = case_when(name == "b_Intercept" ~ paste0("a) ","\u03b1"),
                                TRUE ~ paste0("b) ","\u03b2")))

#4) Append true values. Calculate coverage
reg_coverage = bind_rows(reg_summaries) %>% 
  left_join(true_values) %>% 
  mutate(contains = case_when(true_values >= .lower & true_values <= .upper ~ 1,
                              TRUE ~ 0)) %>% 
  arrange(name_greek, true_values, .lower) %>% 
  ungroup %>% 
  mutate(replicate = rep(1:1000, 2))

saveRDS(reg_coverage, file = "posteriors/fig3_posterior_summaries.rds")
