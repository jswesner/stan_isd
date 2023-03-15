library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)

# get data
sim_data = readRDS(file = "data/sim_data.rds") %>% 
  bind_rows() %>% 
  rename(true_lambda = lambda,
         dw = x)

# load posteriors
recover_sims = readRDS(file = "posteriors/recover_sims_counts.rds")

# median lambdas
posts_medians = recover_sims %>% 
  group_by(true_lambda) %>% 
  median_qi(lambda)

# plot isd's --------------------------------------------------------------

# sample dw weighted by density
nsamples = 1000

dat_sims = sim_data %>% 
  left_join(posts_medians) %>% 
  group_by(true_lambda, xmax, sample_size) %>% 
  sample_n(nsamples, weight = counts, replace = T) %>% 
  select(dw, true_lambda, xmin, xmax, counts, lambda, .lower, .upper) 

dat_toplot = dat_sims %>% 
  # filter(sample_int %in% c(id)) %>% 
  group_by(true_lambda) %>% 
  arrange(desc(dw)) %>% 
  mutate(y_order = 1:nsamples - 1,
         true_lambda = round(true_lambda, 1)) %>% 
  ungroup() 

# estimate isd line and CrI
dat_split = dat_sims %>% 
  # filter(sample_int %in% c(id)) %>% 
  group_by(true_lambda) %>% 
  group_split

xy.PLB = NULL
for(i in 1:length(dat_split)) {
  true_lambda = unique(dat_split[[i]]$true_lambda)
  # site_id = unique(dat_split[[i]]$site_id)
  # year = unique(dat_split[[i]]$year)
  xmin = unique(dat_split[[i]]$xmin)
  xmax = unique(dat_split[[i]]$xmax)
  
  lambda = unique(dat_split[[i]]$lambda)
  .lower = unique(dat_split[[i]]$.lower)
  .upper = unique(dat_split[[i]]$.upper)
  
  x.PLB = seq(xmin, xmax,
              length=nsamples) # x values to plot PLB
  
  y.PLB = (1 - (x.PLB^(lambda + 1) - (xmin^(lambda+1)))/(xmax^(lambda + 1) - (xmin^(lambda+1))))*nsamples
  ymin.PLB = (1 - (x.PLB^(.lower + 1) - (xmin^(.lower+1)))/(xmax^(.lower + 1) - (xmin^(.lower+1))))*nsamples
  ymax.PLB = (1 - (x.PLB^(.upper + 1) - (xmin^(.upper+1)))/(xmax^(.upper + 1) - (xmin^(.upper+1))))*nsamples
  
  xy.PLB[[i]] = tibble(dw = x.PLB, 
                       y_order = y.PLB,
                       ymin = ymin.PLB,
                       ymax = ymax.PLB,
                       xmax = xmax,
                       xmin = xmin) %>%
    mutate(true_lambda = true_lambda,
           # site_id = site_id,
           # year = year,
           lambda = lambda)
}

lines_toplot = bind_rows(xy.PLB) %>% 
  mutate(true_lambda = round(true_lambda, 1)) %>%
  mutate(true_lambda = case_when(true_lambda == min(true_lambda) ~ "a) \u03bb = -2",
                                 true_lambda == max(true_lambda) ~ "c) \u03bb = -1.3",
                                 TRUE ~ "b) \u03bb = -1.3")) 

isd_single_samples_plot = dat_toplot   %>%
  mutate(true_lambda = case_when(true_lambda == min(true_lambda) ~ "a) \u03bb = -2",
                                 true_lambda == max(true_lambda) ~ "c) \u03bb = -1.3",
                                 TRUE ~ "b) \u03bb = -1.3")) %>%
  ggplot(aes(x = dw, y = y_order/1000, group = true_lambda)) + 
  geom_point(shape = 21, size = 0.3) +
  geom_line(data = lines_toplot ) +
  geom_ribbon(data = lines_toplot  , 
              aes(ymin = ymin/1000, ymax = ymax/1000), alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~true_lambda) +
  theme_default() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)) +
  labs(y = "Proportion of values \u2265 x",
       x = "Individual dry mass (mg)",
       color = "") +
  guides(color = "none") +
  coord_cartesian(ylim = c(1e-04, NA))

isd_single_samples_plot

saveRDS(isd_single_samples_plot, file = "plots/isd_single_samples_plot.rds")
ggsave(isd_single_samples_plot, file = "plots/isd_single_samples_plot.jpg", 
       width = 7, height = 7, dpi = 400)


lambda_labels = dat_toplot %>%
  mutate(letters = case_when(true_lambda == min(true_lambda) ~ "d)",
                                 true_lambda == max(true_lambda) ~ "b)",
                                 TRUE ~ "c)")) %>%
  mutate(lambda_value = case_when(true_lambda == min(true_lambda) ~ "\u03bb = -2",
                                 true_lambda == max(true_lambda) ~ "\u03bb = -1.3",
                                 TRUE ~ "\u03bb = -1.6")) %>% 
  mutate(dw = 2, y_order = 70) %>% 
  ungroup %>% 
  distinct(dw, y_order, letters, lambda_value)

isd_single_samples_plot_long = dat_toplot %>%
  mutate(letters = case_when(true_lambda == min(true_lambda) ~ "d)",
                                 true_lambda == max(true_lambda) ~ "b)",
                                 TRUE ~ "c)")) %>%
  ggplot(aes(x = dw, y = y_order/1000, group = true_lambda)) + 
  geom_point(shape = 21, size = 0.1) +
  geom_line(data = lines_toplot %>%
              mutate(letters = case_when(true_lambda == min(true_lambda) ~ "d)",
                                             true_lambda == max(true_lambda) ~ "b)",
                                             TRUE ~ "c)")),
            linewidth = 0.2) +
  geom_ribbon(data = lines_toplot %>%
                mutate(letters = case_when(true_lambda == min(true_lambda) ~ "d)",
                                               true_lambda == max(true_lambda) ~ "b)",
                                               TRUE ~ "c)"))  , 
              aes(ymin = ymin/1000, ymax = ymax/1000), alpha = 0.2) +
  geom_text(data = lambda_labels, aes(label = lambda_value),
            size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~letters, ncol = 1) +
  theme_default() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0)) +
  labs(y = "Proportion of values \u2265 x",
       x = "Individual dry mass (mg)",
       color = "") +
  guides(color = "none") +
  coord_cartesian(ylim = c(1e-04, NA))


isd_single_samples_plot_long

saveRDS(isd_single_samples_plot_long, file = "plots/isd_single_samples_plot_long.rds")
ggsave(isd_single_samples_plot_long, file = "plots/isd_single_samples_plot_long.jpg", 
       width = 3, height = 7, dpi = 400)




