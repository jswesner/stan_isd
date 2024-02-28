library(tidyverse)
library(tidybayes)
library(brms)
library(ggthemes)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.

# Figure 2a) Sample size -------------------------------------------------------------
#1) load and wrangle sample size posteriors
sample_size_summary = readRDS(file = "posteriors/fig2a_posterior_summaries.rds") %>% 
  rename(n_sim = n) %>% 
  mutate(true_value = true_lambda) %>% 
  mutate(cov95 = case_when(true_value <= .upper &
                             true_value >= .lower ~ 1,
                           TRUE ~0),
         n_sim_label = case_when(n_sim == 30 ~ "a) n = 30",
                                 n_sim == 100 ~ "b) n = 100",
                                 n_sim == 300 ~ "c) n = 300",
                                 TRUE ~ "d) n = 1000"),
         true_value_label = case_when(true_value == -2 ~ "a) \u03bb = -2",
                                      TRUE ~ "b) \u03bb = -1.6")) 

#2) calculate coverage

sample_size_coverage = sample_size_summary %>% 
  group_by(n_sim, true_value) %>% 
  reframe(cov = round(sum(cov95)/1000,2)) %>% 
  mutate(x = case_when(true_value == -1.6 ~ n_sim + 0.4*n_sim,
                       TRUE ~ n_sim - 0.19*n_sim))

# 3) make labels (for plotting next)
size_labels = 
  sample_size_summary %>% 
  distinct(true_value_label, n_sim) %>% 
  mutate(n_sim_label = paste0("n = ", n_sim))

# 4) make plot
sample_size_plot = sample_size_summary %>%
  ggplot(aes(x = n_sim, y = b_Intercept, group = interaction(n_sim, true_value))) + 
  scale_x_log10() +
  scale_color_colorblind() +
  stat_slab(position = position_dodge(width = 0.1)) +
  labs(y = "\u03bb",
       x = "Sample size",
       color = "95% CrI contains true value?",
       subtitle = "a) Sample Size") +
  theme_default() +
  geom_hline(data = . %>% ungroup %>% distinct(n_sim, true_value), aes(yintercept = true_value)) +
  theme(strip.text.y = element_blank(),
        legend.position = "top") + 
  guides(color = "none") +
  geom_text(data = sample_size_coverage, aes(y = true_value + 0.03, label = cov, x = x),
            size = 2) +
  geom_point(shape = 1, size = 0.05,
             aes(color = as.factor(cov95)), position = position_dodge(width = 0.1)) 

# 5) save plot
ggview::ggview(sample_size_plot, width = 6.5, height = 6)
save_plot_and_data(sample_size_plot, file_name = "plots/fig2a_sample_size_plot", width = 6.5, height = 6, dpi = 500)

# 6) calculate bias
sample_size_summary %>% 
  mutate(bias = b_Intercept - true_value) %>% 
  group_by(true_value, n_sim) %>% 
  reframe(mean_bias = mean(bias),
          sd_bias = sd(bias)) 

# 7) precision
sample_size_summary %>% 
  group_by(n_sim, true_value) %>% 
  mutate(range = max(b_Intercept) - min(b_Intercept)) %>% 
  distinct(n_sim, true_value, range)


# Figure 2b) Size range --------------------------------------------------------------

#1) wrangle lambda estimates 
xmax_lambdas = readRDS(file = "posteriors/fig2b_posterior_summaries.rds") 

#2) make plot

xmax_labels = 
  xmax_lambdas %>%
  filter(xmin == 1) %>% 
  distinct(xmin, xmax, true_value, orders_mag) %>% 
  mutate(xmin_label = paste0("xmin = ", xmin),
         xmax_label = paste0("xmax = ", xmax),
         true_value_label = case_when(true_value == -2 ~ "a) \u03bb = -2",
                                      TRUE ~ "b) \u03bb = -1.6"))

xmax_coverage = xmax_lambdas %>% 
  group_by(true_value, orders_mag) %>% 
  add_tally() %>% 
  reframe(cov = sum(cov95 == "yes"),
          n = n) %>% 
  distinct(true_value, n, cov, orders_mag) %>% 
  mutate(cov = round(cov/n, 2)) %>% 
  mutate(x = case_when(true_value == -1.6 ~  orders_mag - 0.12,
                       TRUE ~ orders_mag - 0.28))

range_of_sizes = xmax_lambdas %>% 
  # filter(replicate <= 1000) %>%
  ggplot(aes(x = orders_mag, y = b_Intercept, group = interaction(orders_mag, true_value))) +
  # facet_wrap(. ~ true_value, ncol = 1) +
  labs(y = "\u03bb",
       x = "Range of body sizes\n(orders of magnitude)",
       color = "95% CrI contains true value?",
       subtitle = "b) Size Range") +
  scale_color_colorblind() +
  theme_default() +
  stat_slab(position = position_dodge(width = 0.4)) +
  geom_hline(data = . %>% ungroup %>% distinct(xmin, true_value), aes(yintercept = true_value)) +
  theme(strip.text.y = element_blank(),
        legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size=3))) +
  geom_text(data = xmax_coverage, aes(y = true_value + 0.03, label = cov, x = x),
            size = 2.4) +
  geom_point(shape = 1, size = 0.05,
             aes(group = true_value, color = cov95), position = position_dodge(width = 0.4)) 

#3) save plot
ggview::ggview(range_of_sizes, width = 6.5, height = 6)
save_plot_and_data(range_of_sizes, file_name = "plots/fig2b_range_of_sizes", width = 5, height = 9, dpi = 500)

#4) calculate bias
xmax_lambdas %>% 
  mutate(bias = b_Intercept - true_value) %>% 
  group_by(orders_mag, true_value) %>% 
  reframe(mean_bias = mean(bias),
          sd_bias = sd(bias)) 

#5) precision
xmax_lambdas %>% 
  group_by(orders_mag, true_value) %>% 
  mutate(range = max(b_Intercept) - min(b_Intercept)) %>% 
  distinct(orders_mag, true_value, range)


# Combine plots -----------------------------------------------------------
sample_size_plot = readRDS(file = "plots/fig2a_sample_size_plot.rds")
range_of_sizes = readRDS(file = "plots/fig2b_range_of_sizes.rds")

library(patchwork)

sample_size_range = (sample_size_plot + ylim(-3, -1.2) +
  guides(color = "none")) / (range_of_sizes + ylim(-3, -1.2) + guides(color = "none"))

ggview::ggview(sample_size_range, width = 5.5, height = 9)
save_plot_and_data(sample_size_range, file_name = "ms/Figure2", file_type = "pdf", width = 5.5, height = 9, dpi = 500)

