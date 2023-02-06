library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
rstan_options("auto_write" = TRUE)

set.seed(1234987)


# single sample -----------------------------------------------------------
# precompile code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

n_sim = 1000
xmax = 1000
xmin = 1
n_bs = 10
b = seq(-2.2, -1.2, length.out = n_bs)

x = sizeSpectra::rPLB(1000, -1.4, xmin = xmin, xmax = xmax)

sim_data = tibble(xmax = xmax,
       xmin = xmin,
       counts = 1) %>% 
  expand_grid(x = x)


stan_dat <- list(x = sim_data$x,
                 N = nrow(sim_data),
                 counts = sim_data$counts,
                 xmax = sim_data$xmax,
                 xmin = sim_data$xmin)

fit <- sampling(object = fit_model,
                data = stan_dat,
                iter = 100,
                chains = 1,
                open_progress = F,
                verbose = F)



# single samples ----------------------------------------------------------

# precompile code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

# simulate
n_sim = 1000
xmax = 1000
xmin = 1
n_bs = 10
b = seq(-2.2, -1.2, length.out = n_bs)
b_true = b

sim_b = tibble(xmax = xmax, xmin = xmin,
               b = b) %>% 
  mutate(group = as.integer(1:nrow(.)))

# simulate single samples ----------------------------------------------------------

# precompile code
fit_model = stan_model("models/b_paretocounts_singlesample.stan")

# simulate
n_sim = 1000
xmax = 1000
xmin = 1
n_bs = 10
b = seq(-2.2, -1.2, length.out = n_bs)
b_true = b

sim_b = tibble(xmax = xmax, xmin = xmin,
               b = b) %>% 
  mutate(group = as.integer(1:nrow(.)))

sim_data = sim_b %>% 
  expand_grid(individual = 1:n_sim) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1),
         sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(b, xmin, xmax, group) %>% 
  mutate(x = round(sims, 3)) %>% 
  count(x, name = "counts") %>% 
  # mutate(group = as.factor(b)) %>% 
  group_by(group) %>% 
  group_split()

# preallocate
posts_single = list()

# fit
for (i in 1:n_bs) {
  b = unique(sim_data[[i]]$b)
  stan_dat <- list(x = sim_data[[i]]$x,
                   N = nrow(sim_data[[i]]),
                   counts = sim_data[[i]]$counts,
                   xmax = sim_data[[i]]$xmax,
                   xmin = sim_data[[i]]$xmin)
  
  fit <- sampling(object = fit_model,
                  data = stan_dat,
                  iter = 1000,
                  chains = 2,
                  open_progress = F,
                  verbose = F)
  
  posts_single[[i]] = as_draws_df(fit) %>% 
    mutate(true_value = b)
      
}

recover_sims = bind_rows(posts_single)

saveRDS(recover_sims, file = "posteriors/recover_sims_counts.rds")

parameter_recovery_withcounts = recover_sims %>% 
  group_by(true_value) %>% 
  median_qi(b_exp) %>% 
  rename(b_modeled = b_exp)

write_csv(parameter_recovery_withcounts, file = "tables/parameter_recovery_withcounts.csv")

# varying intercepts ------------------------------------------------------

mod = stan_model("models/stan_spectra_mod_ibts_temperature_randonly.stan")

varint_sims = sim_b %>% 
  mutate(mean = mean(b),
         offset = b - mean,
         sd = sd(b))

sim_data_varint = bind_rows(sim_data)

# convert to a list for stan
stan_dat <- list(x = sim_data_varint$x,
                 N = nrow(sim_data_varint),
                 counts = sim_data_varint$counts,
                 # mat_s = sim_data_varint$predictor,
                 n_years = length(unique(sim_data_varint$group)),
                 year = as.integer(as.factor(sim_data_varint$group)),
                 xmax = sim_data_varint$xmax,
                 xmin = sim_data_varint$xmin)

# fit model with varying intercepts
fit_varint <- sampling(object = mod,
                       data = stan_dat,
                       iter = 1000,
                       cores = 2,
                       chains = 2,
                       open_progress = F,
                       verbose = F)

saveRDS(fit_varint, file = "models/fit_varint.rds")

# summarize
var_int_summary = fit_varint %>% as_draws_df() %>% 
  pivot_longer(cols = c(sigma_year, a)) %>% 
  group_by(name) %>% 
  median_qi(value) %>% 
  mutate(true_value = c(intercept, var_int_sd))

write_csv(var_int_summary, file = "tables/var_int_summary.csv")

# regression with varying intercepts ------------------------------------------------------

fit_varint_regression = list()

for(i in 1:40){
  intercept = -1.5
  var_int_sd = 0.1
  beta = -0.1
  n_ind = 1000
  
  lambda_sims = tibble(group = (1:3)) %>% 
    mutate(offset = rnorm(nrow(.), 0, var_int_sd),
           group_b = intercept + offset) %>% 
    expand_grid(predictor = -seq(-2, 2, length.out = 3)) %>% 
    expand_grid(replicates = 1:2) %>% 
    mutate(b = intercept + offset + predictor*beta,
           id = 1:nrow(.))
  
  sim_data = lambda_sims %>% expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = 1, xmax = 1000,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(group, id, xmin, xmax) %>% 
    add_count(x, name = "counts") %>% 
    ungroup
  
  
  # convert to a list for stan
  stan_dat <- list(x = sim_data$x,
                   N = nrow(sim_data),
                   counts = sim_data$counts,
                   mat_s = sim_data$predictor,
                   n_years = length(unique(sim_data$group)),
                   year = as.integer(as.factor(sim_data$group)),
                   xmax = sim_data$xmax,
                   xmin = sim_data$xmin)
  
  
  mod = stan_model("models/stan_spectra_mod_varint_regression.stan")
  
  # fit model with varying intercepts
  fit_varint_regression[[i]] <- sampling(object = mod,
                                         data = stan_dat,
                                         iter = 1000,
                                         # cores = 2,
                                         chains = 2,
                                         open_progress = F,
                                         verbose = F)
}


saveRDS(fit_varint_regression, file = "models/fit_varint_regression.rds")

fit_varint_regression = readRDS(file = "models/fit_varint_regression.rds")

fit_reg_posts = list()

for(i in 1:length(fit_varint_regression)){
  fit_reg_posts[[i]] = brms::as_draws_df(fit_varint_regression[[i]]) %>% 
    mutate(model_run = i)
}


# Compare with and without varying intercepts ------------------------------------
intercept = -1.5
beta = -0.1
n_ind = 300
set.seed(234234)

lambda_sims_reg = tibble(predictor = -seq(-2, 2, length.out = 10),
                     group = 1:10) %>% 
  mutate(b = intercept + predictor*beta,
         id = 1:nrow(.),
         b = b + rgamma(nrow(.), 1, 8))


sim_data = lambda_sims_reg %>% expand_grid(individual = 1:n_ind) %>% 
  mutate(xmin = 1, xmax = 1000,
         u = runif(nrow(.), min = 0, max = 1),
         x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
  group_by(group, id, xmin, xmax) %>% 
  add_count(x, name = "counts") %>% 
  ungroup


# convert to a list for stan
stan_dat <- list(x = sim_data$x,
                 N = nrow(sim_data),
                 counts = sim_data$counts,
                 mat_s = sim_data$predictor,
                 n_years = length(unique(sim_data$group)),
                 year = as.integer(as.factor(sim_data$group)),
                 xmax = sim_data$xmax,
                 xmin = sim_data$xmin)


stan_spectra_mod_ibts_temperature_norand = stan_model("models/stan_spectra_mod_ibts_temperature_norand.stan")
mod = stan_model("models/stan_spectra_mod_varint_regression.stan")

# fit model without varying intercepts
fit_single_regression <- sampling(object = stan_spectra_mod_ibts_temperature_norand,
                                       data = stan_dat,
                                       iter = 1000,
                                       # cores = 2,
                                       chains = 2,
                                       open_progress = F,
                                       verbose = F)

saveRDS(fit_single_regression, file = "models/fit_single_regression.rds")

# fit model with varying intercepts
fit_varint_single_regression <- sampling(object = mod,
                                  data = stan_dat,
                                  iter = 1000,
                                  # cores = 2,
                                  chains = 2,
                                  open_progress = F,
                                  verbose = F)

saveRDS(fit_varint_single_regression, file = "models/fit_varint_single_regression.rds")





# bias --------------------------------------------------------------------
# single models--------------------------------------------
single_mods = readRDS(file = "posteriors/recover_sims_counts.rds") %>% 
  rename(value = b_exp) %>% 
  mutate(model = "Separate Models") %>% 
  mutate(parameter = "lambda")

single_mod_summary = single_mods %>% 
  group_by(model, true_value, parameter) %>% 
  median_qi(value) 

# varying intercepts
true_varying_int = varint_sims %>% distinct(mean, sd) %>% 
  rename(a = mean, 
         sigma_year = sd) %>% 
  pivot_longer(everything(), names_to = "parameter", values_to = "true_value")

var_intercept = readRDS("models/fit_varint.rds") %>% 
  as_draws_df() %>% 
  pivot_longer(cols = c(sigma_year, a)) %>% 
  mutate(model = "Single Model with Varying Intercepts") %>% 
  rename(parameter = name) %>% 
  left_join(true_varying_int)

parameter_recovery_table = bind_rows(single_mods, var_intercept) %>%
  group_by(model, true_value, parameter) %>% 
  median_qi(value) %>% 
  select(model, parameter, true_value, value, .lower, .upper) %>% 
  rename(q2.5 = .lower,
         q50 = value,
         q97.5 = .upper) %>% 
  select(model, parameter, true_value, q2.5, q50, q97.5) %>% 
  mutate(parameter = case_when(parameter == "lambda" ~ "\u03BB",
                               parameter == "a" ~ "\u03B1",
                               parameter == "sigma_year" ~ paste0("\u03C3", "_[group]")))
  
write_csv(parameter_recovery_table, file = "tables/parameter_recovery.csv")


var_int_groups = readRDS("models/fit_varint.rds") %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("alpha_raw")) %>% 
  mutate(value = a + sigma_year*value) %>% 
  group_by(name) %>% 
  median_qi(value) %>% 
  mutate(group = as.integer(parse_number(name)),
         model = "Varying Intercepts",
         parameter = "\u03BB") %>% 
  left_join(varint_sims %>% select(group, b) %>% rename(true_value = b))

single_and_varint_plot = bind_rows(var_int_groups, single_mod_summary) %>% 
  ggplot(aes(x = true_value, y = value, color = model, shape = model)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper),
                  position = position_dodge(width = 0.04),
                  size = 0.2) +
  geom_abline() +
  scale_color_grey(end = 0.6) + 
  theme_default() + 
  labs(y = "Modeled Value",
       x = "True Value",
       color = "Models",
       shape = "Models") + 
  ylim(-2.3, -1.1) + 
  xlim(-2.3, -1.1)

ggview(single_and_varint_plot, width = 5, height = 3, units = "in")
ggsave(single_and_varint_plot, width = 5, height = 3, units = "in",
       file = "plots/single_and_varint_plot.jpg")
saveRDS(single_and_varint_plot, file = "plots/single_and_varint_plot.rds")



# varying intercepts regression-------------------------------------------
true_reg_values = tibble(a = -1.5,
                     sigma_year = 0.1,
                     beta_mat = -0.1) %>% 
  pivot_longer(everything(), values_to = "true_value")

# check convergence
rhats = list()
for(i in 1:40){
  rhats[[i]] = summary(fit_varint_regression[[i]])$summary[,"Rhat"][c(1,2,6)]  %>% 
    as_tibble() %>% 
    mutate(model_run = i,
           name = c("a", "beta_mat", "sigma_year"))
}

rhats = bind_rows(rhats) %>% 
  group_by(model_run) %>% 
  mutate(converged = case_when(max(value) > 1.1 ~ "no",
                               TRUE ~ "yes"))

var_reg_summary = bind_rows(fit_reg_posts) %>% 
  select(-lp__, -contains("alpha_raw")) %>% 
  pivot_longer(cols = c(-.chain, -.iteration, -.draw, -model_run)) %>%
  mutate(model = "3) Varying Intercepts Regression") %>% 
  left_join(true_reg_values) %>% 
  mutate(bias = value - true_value) %>% 
  left_join(rhats %>% ungroup %>% distinct(converged, model_run)) %>% 
  filter(converged == "yes")

hdi = var_reg_summary %>% 
  group_by(model_run, name) %>% 
  median_qi(value) %>% 
  select(model_run, name, .lower, .upper, value) %>% 
  left_join(true_reg_values) %>% 
  mutate(true_inside_hdi = case_when(true_value < .lower ~ "no",
                                     true_value > .upper ~ "no",
                                     TRUE ~ "yes"),
         diff = value - true_value)

var_reg_bias = var_reg_summary %>%
  left_join(hdi %>% select(model_run, name, true_inside_hdi)) %>% 
  group_by(model_run) %>% 
  mutate(model_run_new = cur_group_id())

library(ggridges)
library(viridis)
library(ggview)

plot_linear_model_bias = var_reg_bias %>%
  mutate(greek_name = case_when(name == "a" ~ "\u03B1",
                                name == "beta_mat" ~ "\u03B2",
                                TRUE ~ paste0("\u03c3", "_[group]"))) %>% 
  ggplot(aes(x = value, y = model_run_new, group = model_run, fill = true_inside_hdi)) + 
  geom_density_ridges() +
  facet_wrap(~greek_name, scales = "free_x") + 
  geom_vline(data = . %>% ungroup %>% distinct(name, true_value, greek_name),
             aes(xintercept = true_value)) +
  theme_default() + 
  # scale_fill_viridis_d(option = "F") + 
  scale_fill_grey(start = 0, end = 0.8) +
  labs(y = "Replicate",
       x = "Parameter Value") +
  guides(fill = "none")

ggview(plot_linear_model_bias, width = 6.5, height = 3.5, units = "in")

ggsave(plot_linear_model_bias, width = 6.5, height = 3.5, units = "in",
       file = "plots/plot_linear_model_bias.jpg", dpi = 600)

saveRDS(plot_linear_model_bias, file = "plots/plot_linear_model_bias.rds")


var_reg_bias %>% 
  group_by(model_run, name, true_value) %>% 
  summarize(median = median(value)) %>% 
  mutate(bias = median - true_value) %>% 
  ggplot(aes(x = bias, y = name)) + 
  geom_point() + 
  facet_wrap(~name)


var_reg_bias %>% 
  group_by(model_run, name, true_value) %>% 
  summarize(median = median(value)) %>% 
  mutate(bias = median - true_value) %>% 
  group_by(name) %>% 
  summarize(mean = mean(bias),
            sd = sd(bias))


