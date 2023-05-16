library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
rstan_options("auto_write" = TRUE)

set.seed(1234987)


# single samples ----------------------------------------------------------

# precompile code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

# simulate single samples 

# simulate
n_sim = 1000
xmax = 1000
xmin = 1
# n_bs = 4
# b = seq(-2.2, -1.2, length.out = n_bs)
b = c(-2, -1.6, -1.3)
b_true = b
crops = c(0, 0.01, 0.1, 0.3, 0.5, 0.75) # get the whole data set or the lower or upper 50%, 10%, etc
rep = 1:4
n_bs = length(crops)*length(rep)*length(b_true)

sim_b = tibble(xmax = xmax, xmin = xmin) %>% 
  expand_grid(crop = crops) %>% 
  expand_grid(lambda = b) %>% 
  expand_grid(rep = rep) %>% 
  mutate(group = as.integer(1:nrow(.)),
         sample_size = n_sim) 

sim_smallest = sim_b %>% 
  expand_grid(individual = 1:n_sim) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
  group_by(lambda, xmin, xmax, group, sample_size, crop, rep) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  mutate(limit = xmax*crop,
         greater = case_when(x <= limit ~ 1,
                             crop == 0 ~ 1,
                             TRUE ~ 0),
         sizes = "smallest quantiles") %>%
  filter(greater == 1) 

sim_largest = sim_b %>% 
  expand_grid(individual = 1:n_sim) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
  group_by(lambda, xmin, xmax, group, sample_size, crop, rep) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  mutate(limit = xmax*crop,
         greater = case_when(x >= limit ~ 1,
                             crop == 0 ~ 1,
                             TRUE ~ 0),
         crop = 1-crop,
         sizes = "largest quantiles") %>%
  filter(greater == 1) 

sim_data_tibble = bind_rows(sim_largest, sim_smallest) %>% 
  group_by(crop, lambda, rep, sizes) %>% 
  mutate(group = cur_group_id())


# check
sim_data_tibble %>% 
  ggplot(aes(x = crop, y = x, color = rep)) + 
  geom_jitter() +
  facet_wrap(~sizes)

sim_data = sim_data_tibble %>% 
  group_by(crop, rep, lambda, sizes) %>% 
  mutate(xmin_x = min(x),
         xmax_x = max(x)) %>% 
  group_split()

# saveRDS(sim_data, file = "data/sim_data.rds")

# preallocate
posts_single = list()

# fit
for (i in 1:max(sim_data_tibble$group)) {
  true_lambda = unique(sim_data[[i]]$lambda)
  crop = unique(sim_data[[i]]$crop)
  rep = unique(sim_data[[i]]$rep)
  sizes = unique(sim_data[[i]]$sizes)
  cropped_sample_size = length(sim_data[[i]]$x)
  stan_dat <- list(x = sim_data[[i]]$x,
                   N = nrow(sim_data[[i]]),
                   counts = sim_data[[i]]$counts,
                   xmax = sim_data[[i]]$xmax_x,
                   xmin = sim_data[[i]]$xmin_x)
  
  fit <- sampling(object = fit_model,
                  data = stan_dat,
                  iter = 1000,
                  chains = 2,
                  open_progress = F,
                  verbose = F)
  
  posts_single[[i]] = as_draws_df(fit) %>% 
    mutate(true_lambda = true_lambda,
           crop = crop,
           rep = rep, 
           sizes = sizes,
           cropped_sample_size = cropped_sample_size)
      
}

recover_sims = bind_rows(posts_single)

# saveRDS(recover_sims, file = "posteriors/recover_sims_counts.rds")

parameter_recovery_withcounts = recover_sims %>% 
  group_by(crop, rep, true_lambda, sizes, cropped_sample_size) %>% 
  median_qi(lambda) %>% 
  rename(b_modeled = lambda) %>% 
  mutate(crop = case_when(crop == 0 & grepl("smallest", sizes) ~ 1, TRUE ~ crop)) 

saveRDS(parameter_recovery_withcounts, file = "posteriors/sandbox/parameter_recovery_withcounts.rds")

parameter_recovery_withcounts %>% 
  ggplot(aes(x = crop, y = b_modeled, group = true_lambda)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper, color = as.factor(true_lambda)),
                  position = position_jitter(width = 0.03)) + 
  geom_hline(yintercept = b, linetype = "dotted") + 
  facet_grid(true_lambda ~ sizes) +
  # guides(color = "none") +
  labs(title = "Analyzing just fish leads to bad estimates",
       subtitle = "Dotted line is the true lambda",
       x = "Proportion of sizes kept relative to xmax",
       y = "Lambda",
       color = "True Lambda") + 
  scale_x_log10() +
  scale_color_colorblind() +
  theme_default()



post_methods = bind_rows(fish_with_xminx %>% mutate(min_method = "xmin is smallest in data", crop_method = "just fish"),
          fish_with_global_xmin %>% mutate(min_method = "xmin is smallest possible", crop_method = "just fish"),
          fish_keep_smallest %>% mutate(min_method = "xmin is smallest in data", crop_method = "just inverts"),
          fish_keep_smallest_with_global_xmin %>% mutate(min_method = "xmin is smallest possible", crop_method = "just inverts"))

saveRDS(post_methods, file = "posteriors/sandbox/post_methods.rds")

(
  compare_methods = post_methods %>% 
  ggplot(aes(x = crop, y = b_modeled, color = crop_method)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper),
                  position = position_jitter(width = 0.01)) + 
  geom_hline(yintercept = b, linetype = "dotted") + 
  facet_grid(min_method ~ crop_method) +
  guides(color = "none") +
    labs(title = "Analyzing just fish leads to bad estimates",
         subtitle = "Dotted line is the true lambda",
         x = "Proportion of sizes kept relative to xmax",
         y = "Lambda") + 
    theme_default()
)

ggsave(compare_methods, file = "posteriors/sandbox/compare_methods.jpg")


# write_csv(parameter_recovery_withcounts, file = "tables/parameter_recovery_withcounts.csv")

