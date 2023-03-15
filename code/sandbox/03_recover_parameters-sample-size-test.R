library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
rstan_options("auto_write" = TRUE)

set.seed(1234987)

# single samples ----------------------------------------------------------

# precompile code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

# simulate single samples 

# simulate
n_sim = c(1000, 10000)
xmax = c(1000, 100000)
xmin = 1
n_bs = 3
b = seq(-2.2, -1.2, length.out = n_bs)
b_true = b

sim_b = tibble(xmin = xmin,
               b = b) %>% 
  expand_grid(n_sim = n_sim, 
              xmax = xmax) %>% 
  mutate(group = as.integer(1:nrow(.))) %>% 
  group_by(group) %>%
  uncount(n_sim, .remove = F) %>% 
  mutate(id = row_number()) %>% 
  ungroup

sim_data_sample_size = sim_b %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(b, xmin, xmax, group, n_sim) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  group_by(group) %>% 
  group_split()

saveRDS(sim_data_sample_size, file = "data/sim_data-test-sample-size.rds")

# preallocate
posts_single = list()

ngroup = max(sim_b$group)

# fit
for (i in 1:length(sim_data_sample_size)) {
  b = unique(sim_data_sample_size[[i]]$b)
  n_sim = unique(sim_data_sample_size[[i]]$n_sim)
  xmax = unique(sim_data_sample_size[[i]]$xmax)
  
  stan_dat <- list(x = sim_data_sample_size[[i]]$x,
                   N = nrow(sim_data_sample_size[[i]]),
                   counts = sim_data_sample_size[[i]]$counts,
                   xmax = sim_data_sample_size[[i]]$xmax,
                   xmin = sim_data_sample_size[[i]]$xmin)
  
  fit <- sampling(object = fit_model,
                  data = stan_dat,
                  iter = 1000,
                  chains = 2,
                  open_progress = F,
                  verbose = F)
  
  posts_single[[i]] = as_draws_df(fit) %>% 
    mutate(true_value = b,
           n_sim = n_sim,
           xmax = xmax)
      
}

recover_sims_sample_size = bind_rows(posts_single)

saveRDS(recover_sims_sample_size, file = "posteriors/recover_sims-sample-size.rds")


