library(tidyverse)
library(tidybayes)
library(brms)
library(ggthemes)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.

# Size range --------------------------------------------------------------

#1) wrangle lambda estimates 
xmax_lambdas = readRDS(file = "models/coverage_models/all_lambda_summaries.rds") %>% 
  mutate(model = "separate_models", pars = "lambda",
         logxmin = log10(xmin),
         logxmax = log10(xmax),
         orders_mag = round(logxmax-logxmin), 0) %>% 
  rename(true_value = true_lambda) %>% 
  filter(true_value %in% c(-2, -1.6)) %>% 
  mutate(cov95 = case_when(true_value > `2.5%` & true_value <`97.5%` ~ "yes", TRUE ~ "no"),
         cov50 = case_when(true_value > `25%` & true_value <`75%` ~ "yes", TRUE ~ "no"),
         cov99 = case_when(true_value > `0.5%` & true_value <`99.5%` ~ "yes", TRUE ~ "no")) 

xmax_coverage = xmax_lambdas %>% 
  group_by(true_value, orders_mag) %>% 
  add_tally() %>% 
  reframe(cov = sum(cov95 == "yes"),
          n = n) %>% 
  distinct(true_value, n, cov, orders_mag) %>% 
  mutate(cov = round(cov/n, 2)) %>% 
  mutate(x = case_when(true_value == -1.6 ~ 0.29 + orders_mag,
                       TRUE ~ orders_mag - 0.28))

#2) make plot

xmax_labels = 
  xmax_lambdas %>%
  filter(xmin == 1) %>% 
  distinct(xmin, xmax, true_value, orders_mag) %>% 
  mutate(xmin_label = paste0("xmin = ", xmin),
         xmax_label = paste0("xmax = ", xmax),
         true_value_label = case_when(true_value == -2 ~ "a) \u03bb = -2",
                                      TRUE ~ "b) \u03bb = -1.6"))

range_of_sizes = xmax_lambdas %>%
  filter(xmin == 1) %>% 
  arrange(orders_mag, true_value, `2.5%`) %>% 
  mutate(replicate = rep(1:1000, 10),
         true_value_label = case_when(true_value == -2 ~ "a) \u03bb = -2",
                                      TRUE ~ "b) \u03bb = -1.6")) %>% 
  ggplot(aes(y = replicate)) + 
  geom_linerange(linewidth = 0.3,
                 aes(x = mean, xmin = `2.5%`, xmax = `97.5%`,
                     color = as.factor(cov95))) +
  scale_color_colorblind() +
  lemon::facet_rep_grid(orders_mag ~ true_value_label,
                        scales = "free") +
  geom_vline(aes(xintercept = true_value)) +
  labs(x = "\u03bb",
       y = "Simulation",
       color = "CrI contains true value?") +
  ggthemes::scale_color_colorblind() +
  geom_text(data = xmax_labels, aes(label = xmax_label, x = true_value - 0.45, y = 620, hjust = 0),
            size = 2.4) +
  geom_text(data = xmax_labels, aes(label = xmin_label, x = true_value - 0.45, y = 700, hjust = 0),
            size = 2.4) +
  theme_default() +
  theme(legend.position = "top",
        strip.text.x = element_text(hjust = 0,
                                    vjust = 3),
        strip.text.y = element_blank()) +
  scale_y_continuous(breaks = c(1, 500, 1000)) +
  # xlim(-3, -1) +
  NULL

# range_of_sizes = xmax_lambdas %>% 
#   # filter(replicate <= 1000) %>% 
#   ggplot(aes(x = orders_mag, y = mean, group = interaction(orders_mag, true_value))) + 
#   geom_point(shape = 21, size = 0.01,
#              aes(group = true_value, color = cov95), position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.03)) +
#   # facet_wrap(. ~ true_value, ncol = 1) +
#   labs(y = "\u03bb",
#        x = "Range of body sizes\n(orders of magnitude)",
#        color = "95% CrI contains true value?",
#        subtitle = "b) Size Range") +
#   scale_color_colorblind() +
#   theme_default() +
#   geom_hline(data = . %>% ungroup %>% distinct(xmin, true_value), aes(yintercept = true_value)) +
#   theme(strip.text.y = element_blank(),
#         legend.position = "top") + 
#   guides(color = guide_legend(override.aes = list(size=3))) +
#   geom_text(data = xmax_coverage, aes(y = true_value + 0.03, label = cov, x = x),
#             size = 2) + 
#   stat_summary(
#     geom = "point",
#     fun = "mean",
#     col = "black",
#     size = 2,
#     shape = 21,
#     position = position_dodge(width = 0.2)
#   )

#3) save plot
ggview::ggview(range_of_sizes, width = 5, height = 9)
save_plot_and_data(range_of_sizes, file_name = "plots/fig2b_range_of_sizes", width = 5, height = 9, dpi = 500)

#4) calculate bias
xmax_lambdas %>% 
  mutate(bias = mean - true_value) %>% 
  group_by(orders_mag, true_value) %>% 
  reframe(mean_bias = mean(bias),
          sd_bias = sd(bias)) 

#5) precision
xmax_lambdas %>% 
  group_by(orders_mag, true_value) %>% 
  mutate(range = max(mean) - min(mean)) %>% 
  distinct(orders_mag, true_value, range)

# Sample size -------------------------------------------------------------
#1) load sample size posteriors
sample_size_posts_df = readRDS(file = "posteriors/sample_size_posts_df.rds")

#2) Summarize coverage
true_value = rep(rep(c(-2, -1.6), each = 4), 1000) # did this by hand - double check

sample_size_summary = sample_size_posts_df %>% 
  rename(n_sim = n) %>% 
  mutate(true_value = true_value) %>% 
  mutate(cov95 = case_when(true_value <= .upper &
                                true_value >= .lower ~ 1,
                              TRUE ~0),
         n_sim_label = case_when(n_sim == 30 ~ "a) n = 30",
                                 n_sim == 100 ~ "b) n = 100",
                                 n_sim == 300 ~ "c) n = 300",
                                 TRUE ~ "d) n = 1000"),
         true_value_label = case_when(true_value == -2 ~ "a) \u03bb = -2",
                                      TRUE ~ "b) \u03bb = -1.6")) 

sample_size_coverage = sample_size_summary %>% 
  group_by(n_sim, true_value) %>% 
  reframe(cov = round(sum(cov95)/1000,2)) %>% 
  mutate(x = case_when(true_value == -1.6 ~ n_sim + 0.4*n_sim,
                       TRUE ~ n_sim - 0.25*n_sim))

size_labels = 
  sample_size_summary %>% 
  distinct(true_value_label, n_sim) %>% 
  mutate(n_sim_label = paste0("n = ", n_sim))

sample_size_plot = sample_size_summary %>%
  arrange(n_sim, true_value, .lower) %>% 
  mutate(replicate = rep(1:1000, 8)) %>% 
  ggplot(aes(y = replicate)) + 
  geom_linerange(linewidth = 0.3,
                 aes(x = b_Intercept, xmin = .lower, xmax = .upper,
                     color = as.factor(cov95))) +
  scale_color_colorblind() +
  lemon::facet_rep_grid(n_sim ~ true_value_label, 
                        space = "free") +
  geom_vline(aes(xintercept = true_value)) +
  labs(x = "\u03bb",
       y = "Simulation",
       color = "CrI contains true value?") +
  ggthemes::scale_color_colorblind() +
  geom_text(data = size_labels, aes(label = n_sim_label, x = -3.2, y = 120),
            size = 3) +
  theme_default() +
  theme(legend.position = "top",
        strip.text.x = element_text(hjust = 0,
                                    vjust = 3),
        strip.text.y = element_blank()) +
  scale_y_continuous(breaks = c(1, 500, 1000),
  ) +
  NULL


#3) save plot
ggview::ggview(sample_size_plot, width = 5, height = 6)
save_plot_and_data(sample_size_plot, file_name = "plots/fig2a_sample_size_plot", width = 5, height = 6, dpi = 500)

#4) calculate bias
sample_size_summary %>% 
  mutate(bias = b_Intercept - true_value) %>% 
  group_by(true_value, n_sim) %>% 
  reframe(mean_bias = mean(bias),
          sd_bias = sd(bias)) 

#5) precision
sample_size_summary %>% 
  group_by(n_sim, true_value) %>% 
  mutate(range = max(b_Intercept) - min(b_Intercept)) %>% 
  distinct(n_sim, true_value, range)

# Combine plots -----------------------------------------------------------
sample_size_plot = readRDS(file = "plots/fig2a_sample_size_plot.rds")
range_of_sizes = readRDS(file = "plots/fig2b_range_of_sizes.rds")

library(patchwork)

sample_size_range = sample_size_plot + ylim(-3, -1.2) + guides(color = "none") + range_of_sizes + ylim(-3, -1.2) + guides(color = "none")

ggview::ggview(sample_size_range, width = 6.5, height = 4)
save_plot_and_data(sample_size_range, file_name = "ms/fig2range_of_sizes", width = 6.5, height = 4, dpi = 500)

