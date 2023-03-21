library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)

# prior predictive --------------------------------------------------------------

# sample dw weighted by density
nsamples = 1000 # number of body sizes to simulate
xmin = 1
xmax = 1000
lambda_mean = -1.75
lambda_sd = 0.2
lambda_prior = rnorm(30, lambda_mean, lambda_sd)

prior_isd_sims = tibble(lambda = lambda_prior) %>% 
  expand_grid(xmin = xmin,
              xmax = xmax,
              nsamples = 1:nsamples) %>% 
  group_by(lambda) %>% 
  mutate(x.PLB = row_number()) %>% 
  ungroup %>% 
  mutate(y.PLB = (1 - (x.PLB^(lambda + 1) - (xmin^(lambda+1)))/(xmax^(lambda + 1) - (xmin^(lambda+1)))),
         model = "a) Prior")

prior_isd_sims %>% 
  ggplot(aes(x = x.PLB, y = y.PLB, group = lambda)) + 
  geom_line(alpha = 0.1) +
  scale_x_log10() + 
  scale_y_log10()

# prior vs posterior -------------------------------------------------------
recover_sims = readRDS(file = "posteriors/recover_sims_counts.rds")
#unique true lambdas to choose from for an example
unique_lambdas = unique(recover_sims$true_lambda)

posterior_isd_sims = recover_sims %>% 
  filter(true_lambda >= -1.6 & true_lambda <= -1.49) %>% 
  mutate(draw = 1:nrow(.)) %>% 
  filter(draw <= 30) %>% 
  expand_grid(xmin = xmin,
              xmax = xmax,
              nsamples = 1:nsamples) %>% 
  group_by(lambda) %>% 
  mutate(x.PLB = row_number()) %>% 
  ungroup %>% 
  mutate(y.PLB = (1 - (x.PLB^(lambda + 1) - (xmin^(lambda+1)))/(xmax^(lambda + 1) - (xmin^(lambda+1)))),
         model = "b) Posterior")


posterior_summary = posterior_isd_sims %>% 
  group_by(model) %>% 
  summarize(mean = round(mean(lambda), 2),
            sd = round(sd(lambda), 2)) %>% 
  add_row(model = "a) Prior", mean = lambda_mean, sd = lambda_sd) %>% 
  mutate(text = paste0(model, ": \u03bb = ", mean, " \u00b1 ", sd))


prior_post_plot = bind_rows(prior_isd_sims, posterior_isd_sims) %>% 
  left_join(posterior_summary) %>% 
  ggplot(aes(x = x.PLB, y = y.PLB)) + 
  geom_line(aes(group = lambda)) +
  # geom_text(data = posterior_summary, size = 3, aes(x = 50, y = 1, label = text)) +
  facet_wrap(~text) +
  scale_x_log10() + 
  scale_y_log10() + 
  coord_cartesian(ylim = c(1e-05, NA)) + 
  theme_default() + 
  theme(strip.text = element_text(hjust = 0)) +
  labs(y = "Proportion of values \u2265 x",
       x = "Individual dry mass (mg)")

ggview::ggview(prior_post_plot, width = 6, height = 3)
ggsave(prior_post_plot, width = 6, height = 3, units = "in", dpi = 600,
       file = "plots/prior_post_plot.jpg")
saveRDS(prior_post_plot, file = "plots/prior_post_plot.rds")





# prior vs posterior ibts regression --------------------------------------
# make data grid
ibts_data = readRDS(file = "data/count_sims_thin.rds")
min_x = min(ibts_data$mat_s)
max_x = max(ibts_data$mat_s)

x_preds = seq(min_x, max_x, length.out = 20)

# load posteriors
ibts_conds_hierarchical = readRDS(file = "posteriors/ibts_conds_hierarchical.rds")

# simulate priors
set.seed(98234)
nsims = 100
prior_int = rnorm(nsims, -1.5, 0.2)
prior_slope = rnorm(nsims, 0, 0.1)
prior_sigma = rexp(nsims, 2)

prior_ibts = tibble(a = prior_int, 
                    b = prior_slope,
                    sigma = prior_sigma,
                    .draw = 1:nsims) %>% 
  expand_grid(mat_s = x_preds) %>% 
  mutate(lambda = a + b*mat_s,
         model = "c) Prior slope: 0 \u00b1 0.1")

prior_post_ibts = prior_ibts %>% bind_rows(ibts_conds_hierarchical %>% 
                                             filter(.draw <= 100) %>% 
                                             select(mat_s, .draw, lambda) %>% 
                                             mutate(model = "d) Posterior slope: -0.01 \u00b1 0.01"))

prior_post_ibts = prior_post_ibts %>% 
  ggplot(aes(x = mat_s, y = lambda)) +
  geom_line(aes(group = .draw)) + 
  facet_wrap(~model) + 
  theme_default() + 
  theme(strip.text = element_text(hjust = 0)) +
  labs(y = "\u03bb",
       x = "Year (standardized)")


library(patchwork)
prior_post_4panel = prior_post_plot/prior_post_ibts

ggview::ggview(prior_post_4panel, width = 6.5, height = 7, units = "in")
saveRDS(prior_post_4panel, file = "plots/prior_post_4panel.rds")
ggsave(prior_post_4panel, file = "plots/prior_post_4panel.jpg",
       width = 6.5, height = 7, units = "in")




