library(sizeSpectra)
library(rstan)
library(janitor)
library(tidyverse)
library(brms)
library(ggthemes)
library(patchwork)
library(tidybayes)
rstan_options("auto_write" = TRUE)


# reanalyze IBTS data from Edwards et al. 2020 ----------------------------

# analyze regression with varying intercepts ---------------------
data("IBTS_data")

stan_spectra_mod_ibts_temperature = stan_model("models/stan_spectra_mod_ibts_temperature.stan")

count_sims_thin = IBTS_data %>% clean_names() %>%
  group_by(year) %>% 
  mutate(xmin = min(body_mass)) %>% 
  ungroup %>% 
  mutate(xmax = max(body_mass)) %>%
  left_join(IBTS_data  %>% clean_names() %>% distinct(year) %>% 
              mutate(mat_s = (year - mean(year))/sd(year),
                     mean_year = mean(year),
                     sd_year = sd(year))) %>% 
  mutate(group = as.integer(as.factor(year))) %>% 
  group_by(year) %>% 
  sample_frac(1)

stan_data_ibts = list(N = nrow(count_sims_thin),
                      mat_s = count_sims_thin$mat_s,
                      # gpp_s = count_sims_thin$gpp_s,
                      year = as.integer(as.factor(count_sims_thin$year)),
                      n_years = length(unique(count_sims_thin$year)),
                      # n_sites = length(unique(count_sims_thin$site_id_int)),
                      # site = count_sims_thin$site_id_int,
                      counts = count_sims_thin$number,
                      x = count_sims_thin$body_mass,
                      xmin = count_sims_thin$xmin,
                      xmax = count_sims_thin$xmax)

fit_ibts = sampling(object = stan_spectra_mod_ibts_temperature, 
                    data = stan_data_ibts,
                    iter = 1000, chains = 1, cores = 4)

saveRDS(fit_ibts, file = "models/fit_ibts_norand.rds")

fit_ibts = readRDS(file = "models/fit_ibts_norand.rds")

ibts_posts = as_draws_df(fit_ibts) 

ibts_conds = ibts_posts %>% 
  expand_grid(mat_s = seq(min(count_sims_thin$mat_s), max(count_sims_thin$mat_s), length.out = 20)) %>% 
  mutate(mean_year = unique(count_sims_thin$mean_year),
         sd_year = unique(count_sims_thin$sd_year),
         year = mat_s*sd_year + mean_year) %>% 
  mutate(lambda = a + beta_mat*mat_s)

ibts_year_posts = ibts_posts %>% 
  pivot_longer(cols = contains("alpha_raw_year")) %>% 
  mutate(group = parse_number(name)) %>% 
  left_join(count_sims_thin %>% ungroup %>% distinct(group, mat_s, year)) %>% 
  mutate(lambda = a + beta_mat*mat_s + sigma_year*value) %>% 
  group_by(group, mat_s, year) %>% 
  median_qi(lambda)

# analyze regression without varying intercepts ---------------------
data("IBTS_data")

stan_spectra_mod_ibts_temperature_norand = stan_model("models/stan_spectra_mod_ibts_temperature_norand.stan")

count_sims_thin = IBTS_data %>% clean_names() %>%
  group_by(year) %>% 
  mutate(xmin = min(body_mass)) %>% 
  ungroup %>% 
  mutate(xmax = max(body_mass)) %>%
  left_join(IBTS_data  %>% clean_names() %>% distinct(year) %>% 
              mutate(mat_s = (year - mean(year))/sd(year),
                     mean_year = mean(year),
                     sd_year = sd(year))) %>% 
  mutate(group = as.integer(as.factor(year))) %>% 
  group_by(year) %>% 
  sample_frac(1)

stan_data_ibts_norand = list(N = nrow(count_sims_thin),
                      mat_s = count_sims_thin$mat_s,
                      # gpp_s = count_sims_thin$gpp_s,
                      year = as.integer(as.factor(count_sims_thin$year)),
                      n_years = length(unique(count_sims_thin$year)),
                      # n_sites = length(unique(count_sims_thin$site_id_int)),
                      # site = count_sims_thin$site_id_int,
                      counts = count_sims_thin$number,
                      x = count_sims_thin$body_mass,
                      xmin = count_sims_thin$xmin,
                      xmax = count_sims_thin$xmax)

fit_ibts_norand = sampling(object = stan_spectra_mod_ibts_temperature_norand, 
                    data = stan_data_ibts_norand,
                    iter = 2000, chains = 2, cores = 4)

saveRDS(fit_ibts_norand, file = "models/fit_ibts_norand_norand.rds")

fit_ibts_norand = readRDS(file = "models/fit_ibts_norand_norand.rds")

ibts_posts_norand = as_draws_df(fit_ibts_norand) 

ibts_conds_norand = ibts_posts_norand %>% 
  expand_grid(mat_s = seq(min(count_sims_thin$mat_s), max(count_sims_thin$mat_s), length.out = 20)) %>% 
  mutate(mean_year = unique(count_sims_thin$mean_year),
         sd_year = unique(count_sims_thin$sd_year),
         year = mat_s*sd_year + mean_year) %>% 
  mutate(lambda = a + beta_mat*mat_s)

ibts_conds_norand_summary = ibts_conds_norand %>% 
  group_by(year) %>% 
  mean_qi(lambda) %>%
  mutate(method = "b) Bayesian - not hierarchical")

# plot ---------------------------------------
# get results of time series from Edwards et al. 2020
data("fullResults")
data("eight.results.default")

edwards_mle_estimates = fullResults %>% as_tibble() %>% filter(Method == "MLE") %>% 
  rename(lambda = b,
         .lower = confMin,
         .upper = confMax,
         method = Method,
         year = Year) %>% 
  mutate(year_c = year - mean(year)) %>% 
  select(-stdErr) %>% 
  mutate(method = "a) MLE - two steps")

edwards_reg = lm(lambda ~ year_c, data = edwards_mle_estimates)
edwards_summary = summary(edwards_reg)

mle_coefs = edwards_summary$coefficients
mle_confs = confint(edwards_reg)

# Compare

edwards_reg_slope = -0.001
# get se of slope from eq S.1 https://www.int-res.com/articles/suppl/m636p019_supp1.pdf
edwards_reg_sd = (0.0027 - (-0.0047))/(2*1.96)

mle_regression_line = tibble(slope = -0.001,
       low = -0.0047,
       high = 0.0027,
       intercept = mean(edwards_mle_estimates$lambda),
       int_upper = intercept + sd(edwards_mle_estimates$lambda),
       int_lower = intercept - sd(edwards_mle_estimates$lambda)) %>% 
  expand_grid(year = seq(1985, 2015)) %>% 
  mutate(year_c = year - mean(year),
         lambda = intercept + slope*year_c,
         .upper = int_upper + high*year_c,
         .lower = int_lower + low*year_c)

bayes_mle_regression_lines = ibts_conds %>% 
  group_by(mat_s, year) %>% 
  select(!contains("alpha_raw")) %>%
  mean_qi(lambda) %>% 
  mutate(method = "b) Bayesian - hierarchical") %>% 
  bind_rows(edwards_mle_estimates,
            ibts_conds_norand_summary,
            ibts_conds_norand_summary %>% mutate(method = "b) Bayesian - not hierarchical"))

bayes_mle_year_summaries = edwards_mle_estimates %>% 
  mutate(method = "a) MLE - two steps") %>% 
  bind_rows(ibts_year_posts %>% mutate(method = "c) Bayesian - hierarchical"))

bayes_mle_regression_plot = bayes_mle_regression_lines %>%
  ggplot(aes(x = year, y = lambda, ymin = .lower, ymax = .upper)) + 
  geom_ribbon(data = . %>% filter(grepl("Bayes", method)), alpha = 0.2) +
  geom_line(data = . %>% filter(grepl("Bayes", method)),
            color = "black", linewidth = 0.3) +
  geom_pointrange(data = bayes_mle_year_summaries,
                  size = 0.4) +
  geom_smooth(data = . %>% filter(!grepl("Bayes", method)),
              method = "lm", color = "black", size = 0.8) +
  ylim(-2, -1.3) +
  facet_wrap(~method) +
  theme_default() +
  labs(y = "\u03bb",
       x = "Year") +
  theme(strip.text.x = element_text(hjust = 0, margin=margin(l=0))) +
  NULL


saveRDS(bayes_mle_regression_plot, file = "plots/bayes_mle_regression_plot.rds")


# compare years
bayes_mle_year_summaries %>% 
  ggplot(aes(x = year, y = lambda, color = method)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper), position = position_dodge(width = 0.02))

# table -------------------------------------------------------------------
two_years = count_sims_thin %>% 
  distinct(year, mat_s, sd_year, mean_year) %>% 
  filter(year == 1998|year == 2000) 

edwards_reg_slope = -0.001
# get se of slope from eq S.1 https://www.int-res.com/articles/suppl/m636p019_supp1.pdf
edwards_reg_sd = (0.0027 - (-0.0047))

bayes_mle_regression_table = as_draws_df(fit_ibts) %>% 
  select(!contains("alpha_raw")) %>% 
  expand_grid(two_years) %>% 
  mutate(lambda = a + beta_mat*mat_s) %>% 
  select(.draw, year, sd_year, lambda) %>% 
  pivot_wider(names_from = year, values_from = lambda) %>% 
  mutate(slope = (`2000` - `1998`)/2) %>% 
  mean_qi(slope) %>% 
  mutate(model = "Bayesian - one step") %>% 
  add_row(model = "MLE - two steps",
          slope = -0.001,
          .lower = -0.0047,
          .upper = 0.0027) %>% 
  select(model, slope, .lower, .upper) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  rename(mean = slope,
         q2.5 = .lower,
         q97.5 = .upper)

write_csv(bayes_mle_regression_table, file = "tables/bayes_mle_regression_table.csv")
  
