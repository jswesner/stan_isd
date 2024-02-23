library(brms)
library(tidyverse)
library(isdbayes)
library(ggthemes)

# 1) load model
brm_regress_weakprior_rand = readRDS(file = "models/fig4c_mod.rds")

# 2) refit with just priors and save
brm_prior_sims = update(brm_regress_weakprior_rand,
                        sample_prior = "only")

saveRDS(brm_prior_sims, file = "models/figs6_mod.rds")

# 3) Make data to condition on

data_to_condition = tibble(predictor = seq(from = min(brm_prior_sims$data$predictor), 
                                           to = max(brm_prior_sims$data$predictor), 
                                           length.out = 20)) %>% 
  mutate(counts = 1,
         xmin = min(brm_prior_sims$data$x),
         xmax = max(brm_prior_sims$data$x))

# 4) Extract model predictions from prior and posterior
prior_lines = data_to_condition %>% 
  add_epred_draws(brm_prior_sims, re_formula = NA, ndraws = 100) %>% 
  mutate(group = "a) Prior")

posterior_lines = data_to_condition %>% 
  add_epred_draws(brm_regress_weakprior_rand, re_formula = NA, ndraws = 100) %>% 
  mutate(group = "b) Posterior")

# 5) Plot and save
figs6 = bind_rows(prior_lines, posterior_lines) %>% 
  ggplot(aes(x = predictor, y = .epred, group = .draw)) + 
  geom_line(linewidth = 0.1, alpha = 0.5) +
  facet_wrap(~group) +
  labs(y = "\u03bb",
       x = "Predictor") +
  theme(strip.text = element_text(hjust = 0))

ggview::ggview(figs6, width = 6.5, height = 3)
ggsave(figs6, file = "plots/figs6_prior_post.jpg", width = 6.5, height = 3,
       dpi = 500)
