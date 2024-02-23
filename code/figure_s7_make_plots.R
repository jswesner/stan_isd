library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)
library(patchwork)


# 1) load data
sim_data = readRDS(file = "data/sim_data.rds") %>% 
  bind_rows() %>% 
  rename(true_lambda = b,
         dw = x)

# 2) load models
recover_sims = readRDS(file = "posteriors/figs7_posterior_summaries.rds")

# Figure S7a -----------------------------------------------

# 4) resample dw weighted by density
nsamples = 1000

dat_sim_draw = sim_data %>% 
  group_by(true_lambda) %>% 
  sample_n(nsamples, weight = counts, replace = T) %>% 
  right_join(recover_sims %>% filter(.draw <= 10), relationship = "many-to-many") %>% 
  # group_by(true_lambda, sample_size, .draw) %>% 
  # filter(dw <= 1.01)  %>% 
  ungroup %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         dw = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1)))

# 5) add raw data with ".draw" = 0
sim_plus_raw = dat_sim_draw %>% 
  mutate(source = "ypred") %>% 
  bind_rows(sim_data %>% mutate(source = "raw", .draw = 0))

# 6) calculate geometric means
geom_means = sim_plus_raw %>% 
  group_by(true_lambda, .draw) %>% 
  summarize(geom_mean = exp(mean(log(dw))))

# 6) Make fig 7a
pp_check_plot = sim_plus_raw %>% 
  mutate(true_lambda = round(true_lambda, 1),
         source = case_when(source == "raw" ~ "y",
                            TRUE ~ "y_new")) %>%
  mutate(true_lambda = case_when(true_lambda == min(true_lambda) ~ "a) \u03bb = -2",
                                 true_lambda == max(true_lambda) ~ "c) \u03bb = -1.3",
                                 TRUE ~ "b) \u03bb = -1.6")) %>% 
  ggplot(aes(x = .draw, y = dw, color = source, fill = source)) + 
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

# 7) Make fig 7b
geom_mean_plot = geom_means %>% 
  filter(.draw > 0) %>% 
  ggplot(aes(x = geom_mean)) + 
  geom_histogram() + 
  facet_wrap(~true_lambda) + 
  scale_x_log10() + 
  geom_vline(data = geom_means %>% filter(.draw == 1), 
             aes(xintercept = geom_mean)) +
  theme_default() + 
  labs(x = "Geometric mean of body sizes",
       y = "Count") +
  theme(legend.title = element_blank(),
        text = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(hjust = 0)) 

#8) combine and save
pp_check_geom_mean_plot = pp_check_plot/geom_mean_plot

ggsave(pp_check_geom_mean_plot, file = "plots/figs7_ppcheck.jpg", width = 6.5, height = 4, units = "in",
       dpi = 500)



