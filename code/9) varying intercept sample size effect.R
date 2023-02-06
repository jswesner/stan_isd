library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)

# simulate
set.seed(8787)
n_sim = 500
xmax = 1000
xmin = 1
n_bs = 10
b = seq(-2.2, -1.2, length.out = 10)
b_true = b

sim_b = tibble(xmax = xmax, xmin = xmin,
               b = b) %>% 
  mutate(group = as.integer(1:nrow(.)))

sim_data = sim_b %>% 
  expand_grid(individual = 1:n_sim) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(b, xmin, xmax, group) %>% 
  mutate(counts = 1) %>% 
  mutate(sample = "all_samples") %>% 
  group_by(group) %>% 
  add_tally()

mod = stan_model("models/stan_spectra_mod_ibts_temperature_randonly.stan")

varint_sims = sim_b %>% 
  mutate(mean = mean(b),
         offset = b - mean,
         sd = sd(b))

a = sim_data %>% 
  bind_rows(sim_data %>% 
              group_by(group) %>% 
              sample_frac(0.5) %>% 
              mutate(sample = "small_samples") %>% 
              select(-n) %>% 
              group_by(group) %>% 
              add_tally()) %>% 
  bind_rows(sim_data %>% 
              group_by(group) %>% 
              sample_frac(0.1) %>% 
              mutate(sample = "smallest_samples") %>% 
              select(-n) %>% 
              group_by(group) %>% 
              add_tally())


sim_data_varint = a %>% 
  left_join(a %>% ungroup %>% 
  distinct(group, sample) %>% 
  mutate(id = 1:nrow(.)))


# convert to a list for stan
stan_dat <- list(x = sim_data_varint$x,
                 N = nrow(sim_data_varint),
                 counts = sim_data_varint$counts,
                 # mat_s = sim_data_varint$predictor,
                 n_years = length(unique(sim_data_varint$id)),
                 year = as.integer(as.factor(sim_data_varint$id)),
                 xmax = sim_data_varint$xmax,
                 xmin = sim_data_varint$xmin)

# fit model with varying intercepts
fit_varint_samplesize <- sampling(object = mod,
                       data = stan_dat,
                       iter = 1000,
                       cores = 2,
                       chains = 2,
                       open_progress = F,
                       verbose = F)

saveRDS(fit_varint_samplesize, file = "models/fit_varint_sample_size.rds")


# refit without varying intercepts
# precompile code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

sim_data_list = sim_data %>% 
  group_by(b) %>% 
  group_split()

single_fit_posts = list()

for(i in 1:length(sim_data_list)){
  group = unique(sim_data_list[[i]]$group)
  
  stan_dat = list(x = sim_data_list[[i]]$x,
                  N = nrow(sim_data_list[[i]]),
                  counts = sim_data_list[[i]]$counts,
                  # mat_s = sim_data[[i]]$predictor,
                  # n_years = length(unique(sim_data[[i]]$id)),
                  # year = as.integer(as.factor(sim_data[[i]]$id)),
                  xmax = sim_data_list[[i]]$xmax,
                  xmin = sim_data_list[[i]]$xmin)
  
  single_fits = sampling(object = fit_model,
                           data = stan_dat,
                         chains = 2, 
                         iter = 1000)
  
  single_fit_posts[[i]] = as_draws_df(single_fits) %>% 
    mutate(group = group)
}

single_fit_posts
saveRDS(single_fit_posts, file = "models/single_fit_posts.rds")


fit_varint_samplesize = readRDS(file = "models/fit_varint_sample_size.rds")

post_varint_lambdas = as_draws_df(fit_varint_samplesize) %>% 
  pivot_longer(contains("alpha_raw")) %>% 
  mutate(id = parse_number(name)) %>% 
  left_join(sim_data_varint %>% ungroup %>% distinct(n, id, group)) %>%
  mutate(lambda = a + sigma_year*value,
         model = "partial pooling") %>% 
  filter(n == 500) %>% 
  bind_rows(bind_rows(single_fit_posts) %>% 
              mutate(model = "no pooling"))

post_varint_lambdas %>%
  ggplot(aes(x = group, group = interaction(group, model), y = lambda, fill = model)) + 
  geom_violin(position = position_dodge(width = 0.4))


plot_partial_pooling = post_varint_lambdas %>% 
  ggplot(aes(x = as.factor(group), group = id, y = lambda, fill = model)) + 
  stat_halfeye(position = position_dodge(width = 0.5), scale = 0.4 ,size = 0.6) +
  scale_fill_grey() + 
  geom_hline(aes(yintercept = mean(a))) +
  theme_default() +
  labs(x = "Sample",
       y = "\u03bb",
       fill = "Sample Size (n)")

saveRDS(plot_partial_pooling, file = "plots/plot_partial_pooling.rds")  

