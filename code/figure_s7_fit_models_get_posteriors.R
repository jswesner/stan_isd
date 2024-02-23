library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
library(isdbayes)
rstan_options("auto_write" = TRUE)

set.seed(1234987)


# precompile code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

# simulate
n_sim = 1000
xmax = 1000
xmin = 1
b = c(-2, -1.6, -1.3)
b_true = b

sim_b = tibble(xmax = xmax, xmin = xmin,
               b = b) %>% 
  mutate(group = as.integer(1:nrow(.)))

sim_data = sim_b %>% 
  expand_grid(individual = 1:n_sim) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(b, xmin, xmax, group) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  # mutate(group = as.factor(b)) %>% 
  group_by(group) %>% 
  group_split()

saveRDS(sim_data, file = "data/sim_data.rds")

# preallocate
posts_single = list()

# fit
for (i in 1:length(sim_data)) {
  b = unique(sim_data[[i]]$b)
  stan_dat <- list(x = sim_data[[i]]$x,
                   N = nrow(sim_data[[i]]),
                   counts = sim_data[[i]]$counts,
                   xmax = sim_data[[i]]$xmax,
                   xmin = sim_data[[i]]$xmin)
  
  fit <- sampling(object = fit_model,
                  data = stan_dat,
                  iter = 1000,
                  chains = 2,
                  open_progress = F,
                  verbose = F)
  
  posts_single[[i]] = as_draws_df(fit) %>% 
    mutate(true_lambda = b)
  
}

recover_sims = bind_rows(posts_single)

saveRDS(recover_sims, file = "posteriors/figs7_posterior_summaries.rds")
