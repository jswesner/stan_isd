library(brms)
library(tidyverse)
library(tidybayes)
library(janitor)


# single samples ----------------------------------------------------------

# get data
sim_data_sample_size = readRDS(file = "data/sim_data-test-sample-size.rds") %>% 
  bind_rows() %>% 
  rename(true_value = b,
         dw = x)

# load posteriors
recover_sims_sample_size = readRDS(file = "posteriors/recover_sims-sample-size.rds")

# median lambdas
posts_medians_sample_size = recover_sims_sample_size %>% 
  group_by(true_value, n_sim, xmax) %>% 
  median_qi(lambda)

# plot isd's --------------------------------------------------------------

# sample dw weighted by density
nsamples = 100000

dat_sims_sample_size = sim_data_sample_size %>% 
  left_join(posts_medians_sample_size) %>% 
  group_by(true_value, n_sim, xmax) %>% 
  sample_n(nsamples, weight = counts, replace = T) %>% 
  select(dw, true_value, xmin, xmax, counts, lambda, .lower, .upper, n_sim) 

dat_toplot_sample_size = dat_sims_sample_size %>% 
  # filter(sample_int %in% c(id)) %>% 
  group_by(true_value, n_sim, xmax) %>% 
  arrange(desc(dw)) %>% 
  mutate(y_order = 1:nsamples - 1,
         true_value = round(true_value, 1)) 

# estimate isd line and CrI
dat_split_sample_size = dat_sims_sample_size %>% 
  # filter(sample_int %in% c(id)) %>% 
  group_by(true_value, n_sim, xmax) %>% 
  group_split()

xy.PLB_sample_size = NULL
for(i in 1:length(dat_split_sample_size)) {
  true_value = unique(dat_split_sample_size[[i]]$true_value)
  xmin = unique(dat_split_sample_size[[i]]$xmin)
  xmax = unique(dat_split_sample_size[[i]]$xmax)
  n_sim = unique(dat_split_sample_size[[i]]$n_sim)
  
  lambda = unique(dat_split_sample_size[[i]]$lambda)
  .lower = unique(dat_split_sample_size[[i]]$.lower)
  .upper = unique(dat_split_sample_size[[i]]$.upper)
  
  x.PLB = seq(xmin, xmax,
              length=nsamples) # x values to plot PLB
  
  y.PLB = (1 - (x.PLB^(lambda + 1) - (xmin^(lambda+1)))/(xmax^(lambda + 1) - (xmin^(lambda+1))))*nsamples
  ymin.PLB = (1 - (x.PLB^(.lower + 1) - (xmin^(.lower+1)))/(xmax^(.lower + 1) - (xmin^(.lower+1))))*nsamples
  ymax.PLB = (1 - (x.PLB^(.upper + 1) - (xmin^(.upper+1)))/(xmax^(.upper + 1) - (xmin^(.upper+1))))*nsamples
  
  xy.PLB_sample_size[[i]] = tibble(dw = x.PLB, y_order = y.PLB,
                       ymin = ymin.PLB,
                       ymax = ymax.PLB,
                       xmax = xmax,
                       xmin = xmin) %>%
    mutate(true_value = true_value,
           lambda = lambda,
           n_sim = n_sim)
}

lines_toplot_sample_size = bind_rows(xy.PLB_sample_size) %>% 
  mutate(true_value = round(true_value, 1))

isd_sample_size = dat_toplot_sample_size %>%
  filter(true_value <= -2.1) %>% 
  ggplot(aes(x = dw, y = y_order, group = interaction(true_value, n_sim, xmax))) + 
  geom_point(shape = 21, size = 0.3) +
  geom_line(data = lines_toplot_sample_size %>%
              filter(true_value <= -2.1) ) +
  geom_ribbon(data = lines_toplot_sample_size %>%
                filter(true_value <= -2.1)  , aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(n_sim ~ xmax) +
  theme_default() +
  theme(
    strip.background = element_blank(),
    # strip.text.x = element_blank(),
    axis.text = element_blank()) +
  labs(y = "Number of values \u2265 x",
       x = "Individual dry mass (mg)",
       color = "") +
  guides(color = "none")

isd_sample_size

ggsave(isd_single_samples_plot, file = "plots/isd_single_samples_plot.jpg", 
       width = 7, height = 7, dpi = 400)

