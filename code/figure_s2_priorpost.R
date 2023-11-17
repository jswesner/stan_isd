library(tidyverse)
library(rstan)
library(tidybayes)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.

#1) load fitted model
e = readRDS("models/brm_regress_weakprior_rand.rds")

#2) Check priors
e$prior

#3) simulate from prior predictive
iters = 100

prior_sims = tibble(beta = rnorm(iters, 0, 0.5),
       a = rnorm(iters, -1.5, 1)) %>% 
  mutate(.draw = row_number()) %>% 
  expand_grid(predictor = unique(e$data$predictor)) %>% 
  mutate(.epred = a + beta*predictor,
         model = "a) Prior")

#4) get posterior 
posterior_sims = e$data %>% distinct(predictor) %>% 
  mutate(counts = 1, 
         xmin = 1, 
         xmax = 1000) %>% 
  add_epred_draws(e, re_formula = NA, ndraws = iters) %>% 
  mutate(model = "b) Posterior")

#5) combine prior and posterior
prior_post_plot = bind_rows(prior_sims, posterior_sims) %>% 
  ggplot(aes(x = predictor, y = .epred)) + 
  geom_line(aes(group = .draw), linewidth = 0.1) + 
  facet_wrap(~model) + 
  theme(strip.text = element_text(hjust = 0)) + 
  labs(y = "\u03bb",
       x = 'Predictor')

#6) Save plot
ggview::ggview(prior_post_plot, width = 6.5, height = 3)
save_plot_and_data(prior_post_plot, file_name = "ms/figure_s2_prior_post_plot", width = 6.5, height = 3)



