library(brms)
library(tidyverse)
library(isdbayes)


# Figure 2a) sample size -------------------------------------------------------------
# 1) load models
# This loads 1000-7000 fitted models. It is too big for GitHub (500-1000MB), so we have to load it from Zenodo.
# It might take 5-10 minutes
# fig2a_mods = readRDS("models/fig2a_mods.rds")
fig2a_mods = readRDS(url("https://zenodo.org/records/10689765/files/fig2a_mods.rds?download=1"))

# 2) Summarize posteriors and save
temp_stats <- vector("list", length = 1000)

for (i in 1:1000) {
  # Create a nested list for each i iteration
  temp_stats[[i]] <- vector("list", length = 2)
  
  for (b in 1:8) {
    temp_stats[[i]][[b]] <- fig2a_mods[[b]][14, i][[1]] %>%
      as_draws_df() %>%
      tidybayes::median_qi(b_Intercept) %>% 
      mutate(replicate = i,
             n = nrow(fig2a_mods[[b]][2, i][[1]]),
             b = b,
             true_lambda = case_when(b <= 4 ~ -2,
                                     TRUE ~ -1.6))
  }
}

sample_size_posts_df = bind_rows(temp_stats)

# saveRDS(sample_size_posts_df, file = "posteriors/fig2a_posterior_summaries.rds")

# Figure 2b) size range -------------------------------------------------------------

# 1) load models
# This loads 1000-7000 fitted models. It is too big for GitHub (500-1000MB), so we have to load it from Zenodo.
# It might take 5-10 minutes
options(timeout=1000)
# fig2b_mods = readRDS(file = "models/fig2b_mods.rds")
fig2b_mods = readRDS(url("https://zenodo.org/records/10689765/files/fig2b_mods.rds?download=1"))

# 2) Summarize posteriors and save
temp_size <- vector("list", length = 1000)

for (i in 1:1000) {
  # Create a nested list for each i iteration
  temp_size[[i]] <- vector("list", length = 2)
  
  for (b in 1:10) {
    temp_size[[i]][[b]] <- fig2b_mods[[b]][14, i][[1]] %>%
      as_draws_df() %>%
      tidybayes::median_qi(b_Intercept) %>% 
      mutate(replicate = i,
             n = nrow(fig2b_mods[[b]][2, i][[1]]),
             xmin = unique(fig2b_mods[[b]][2, i][[1]]$xmin),
             xmax = unique(fig2b_mods[[b]][2, i][[1]]$xmax),
             b = b,
             true_lambda = case_when(b <= 5 ~ -2,
                                     TRUE ~ -1.6))
  }
}

size_range_posts = bind_rows(temp_size) %>% 
  mutate(model = "separate_models", pars = "lambda",
         logxmin = log10(xmin),
         logxmax = log10(xmax),
         orders_mag = round(logxmax-logxmin), 0) %>% 
  rename(true_value = true_lambda) %>% 
  filter(true_value %in% c(-2, -1.6)) %>% 
  mutate(cov95 = case_when(true_value > .lower & true_value < .upper ~ "yes", TRUE ~ "no"))

# saveRDS(size_range_posts, file = "posteriors/fig2b_posterior_summaries.rds")
