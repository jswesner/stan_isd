library(brms)
library(tidyverse)
library(isdbayes)
library(tidybayes)

# 1) load fitted models
a1.2 = readRDS(file = "models/fit_priors_1.2.rds")
a1.5 = readRDS(file = "models/fit_priors_1.5.rds")
a1.8 = readRDS(file = "models/fit_priors_1.8.rds")
a2.0 = readRDS(file = "models/fit_priors_2.rds")

# 2) combine models
prior_fits = c(a1.2, a1.5, a1.8, a2.0)

# 3) Extract summaries
fit_prior_temp = NULL

for(i in 1:length(prior_fits)) {
  fit_prior_temp[[i]] = as_draws_df(prior_fits[[i]]) %>% 
    mutate(fit = i,
           prior = prior_fits[[i]]$prior$prior) 
}

# 4) Wrangle summaries and save
figs5_posterior_summaries = bind_rows(fit_prior_temp) %>% 
  separate(prior, into = c("prior_mean", "prior_sd"), sep = ",", remove = F) %>% 
  mutate(prior_mean = parse_number(prior_mean),
         prior_sd = parse_number(prior_sd))

saveRDS(figs5_posterior_summaries,
        file = "posteriors/figs5_posterior_summaries.rds")