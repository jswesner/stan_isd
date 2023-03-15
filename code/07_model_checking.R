library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)


# single samples ----------------------------------------------------------
# compare posterior checks with 1000 or 10,000 individual samples
# get data
sim_data = readRDS(file = "data/sim_data.rds") %>% 
  bind_rows() %>% 
  rename(true_lambda = lambda,
         dw = x)


# load posteriors
recover_sims = readRDS(file = "posteriors/recover_sims_counts.rds")

posts_medians = recover_sims %>% 
  group_by(true_lambda) %>% 
  median_qi(lambda)

# pp_check ----------------------------------------------------------------

# sample dw weighted by density
nsamples = 1000

dat_sim_draw = sim_data %>% 
  right_join(recover_sims %>% filter(.draw <= 10), multiple = "all") %>% 
  group_by(true_lambda, sample_size, .draw) %>% 
  # filter(dw <= 1.01) %>% 
  sample_n(nsamples, weight = counts, replace = T) %>% 
  ungroup %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         dw = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1)))

sim_plus_raw = dat_sim_draw %>% mutate(source = "ypred") %>% 
  bind_rows(sim_data %>% mutate(source = "raw", .draw = 0))


sim_plus_raw %>%  
  mutate(.draw = as.factor(.draw),
         .draw = fct_relevel(.draw, "10", "9", "8", "7","6",
                             "5", "4", "3", "2", "1")) %>% 
  ggplot(aes(x = dw)) +
  geom_density(aes(color = source, group = .draw)) + 
  scale_x_log10() + 
  scale_color_grey() +
  theme_default() +
  facet_wrap(~true_lambda) + 
  guides(color = "none")

pp_check_plot = sim_plus_raw %>% 
  mutate(true_lambda = round(true_lambda, 1),
         source = case_when(source == "raw" ~ "y",
                            TRUE ~ "y_new")) %>%
  mutate(true_lambda = case_when(true_lambda == min(true_lambda) ~ "a) \u03bb = -2",
                                 true_lambda == max(true_lambda) ~ "c) \u03bb = -1.3",
                                 TRUE ~ "b) \u03bb = -1.6")) %>% 
  ggplot(aes(x = .draw, y = dw, color = source, fill = source)) + 
  # geom_violin(aes(group = .draw)) +
  # geom_violin(aes(group = .draw)) + 
  geom_jitter(height = 0, width = 0.2, shape = 46)  + 
  geom_boxplot(outlier.shape = NA, aes(group = .draw)) +
  scale_y_log10()  +
  scale_color_grey(end = 0.5) +
  scale_fill_grey(start = 0.4) +
  facet_wrap(~true_lambda) +
  theme_default() + 
  labs(x = "Raw data or draw from posterior (n = 10 draws)",
       y = "Individual dry mass (mg)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(hjust = 0)
        ) 

pp_check_plot

saveRDS(pp_check_plot, file = "plots/pp_check_plot.rds")
ggview::ggview(pp_check_plot, width = 6.5, height = 2, units = "in")
ggsave(pp_check_plot, file = "plots/pp_check_plot.jpg", width = 6.5, height = 2, units = "in",
       dpi = 500)







# simulate summaries of original data -----------------------------------------------------------

# sample dw weighted by density
nsamples = 1000

dat_sim_draw = sim_data %>% 
  right_join(recover_sims %>% filter(.draw <= 300), multiple = "all") %>% 
  group_by(true_lambda, sample_size, .draw) %>% 
  # filter(dw <= 1.01) %>% 
  sample_n(nsamples, weight = counts, replace = T) %>% 
  ungroup %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         ypred = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1)))


ypred_median = sim_data %>% 
  right_join(recover_sims %>% group_by(true_lambda) %>% median_qi(lambda)) %>% 
  ungroup %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         ypred = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1)),
         lower = (u*xmax^(.lower+1) +  (1-u) * xmin^(.lower+1) ) ^ (1/(.lower+1)),
         upper = (u*xmax^(.upper+1) +  (1-u) * xmin^(.upper+1) ) ^ (1/(.upper+1))) %>% 
  select(true_lambda, ypred, lower, upper) %>% 
  group_by(true_lambda) %>% 
  arrange(true_lambda, ypred) %>% 
  mutate(y_order = row_number())


raw_ypred = sim_data %>% 
  group_by(true_lambda) %>% 
  arrange(true_lambda, dw) %>% 
  mutate(y_order = row_number()) %>% 
  left_join(ypred_median)


raw_ypred %>% 
  ggplot(aes(x = dw, y = ypred)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  shape = 21, size = 0.01, linewidth = 0.001,
                  alpha = 0.2) + 
  facet_wrap(~true_lambda) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline()







sim_medians = dat_sim_draw %>% 
  group_by(true_lambda, .draw) %>% 
  summarize(t_ynew = quantile(ypred, probs = 0.5))

raw_medians = sim_data %>% 
  group_by(true_lambda) %>% 
  summarize(t_yraw = quantile(dw, probs = 0.5))

p_values = sim_medians %>% 
  left_join(raw_medians) %>% 
  mutate(diff = t_ynew - t_yraw,
         higher = case_when(diff >= 0 ~ "higher", TRUE ~ "lower"),
         ndraws = max(.draw)) %>% 
  group_by(true_lambda, higher, ndraws) %>% 
  tally() %>% 
  filter(higher == "higher") %>% 
  mutate(p_value = n/ndraws,
         p_value = round(p_value, 2))

sim_medians %>% 
  ggplot(aes(x = t_ynew)) + 
  geom_histogram(bins = 50) + 
  facet_wrap(~true_lambda) + 
  scale_x_log10() + 
  geom_vline(data = raw_medians, aes(xintercept = t_yraw)) + 
  geom_text(data = p_values, aes(x = 15, y = 60, label = p_value))


