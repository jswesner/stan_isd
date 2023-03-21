library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
library(sizeSpectra)
library(ggview)
library(ggridges)
library(patchwork)

# load models
recover_sims_counts = readRDS(file = "posteriors/recover_sims_counts.rds")  #(simulate single samples result)
fit_varint = readRDS(file = "models/fit_varint.rds")                        #(simulate varying intercepts result)
fit_varint_regression = readRDS(file = "models/fit_varint_regression.rds")  #(simulate regression with varying intercepts result)
fit_ibts_hierarchical = readRDS(file = "models/fit_ibts_hierarchical.rds")  #(reanalyze IBTS regression with varying intercepts result)
fit_ibts_norand_norand = readRDS(file = "models/fit_ibts_norand_norand.rds")#(reanalyze IBTS regression without varying intercepts result)
sample_size_sims = readRDS("models/sample_size_sims.rds")                   #(test sample size result)

# load data
varint_sims = readRDS(file = "data/varint_sims.rds")
count_sims_thin = readRDS(file = "data/count_sims_thin.rds") # for mle/bayes comparison with IBTS data


single_mods = recover_sims_counts %>% 
  rename(value = lambda,
         true_value = true_lambda) %>% 
  mutate(model = "Separate Models") %>% 
  mutate(parameter = "lambda")

single_mod_summary = single_mods %>% 
  group_by(model, true_value, parameter) %>% 
  median_qi(value) 
# Table 1 ----------------------------------------------------------------


# varying intercepts
true_varying_int = varint_sims %>% distinct(mean, sd) %>% 
  rename(a = mean, 
         sigma_year = sd) %>% 
  pivot_longer(everything(), names_to = "parameter", values_to = "true_value")

var_intercept = fit_varint %>% 
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

# Table 2 -----------------------------------------------------------------
# grid for estimating slope
two_years = count_sims_thin %>% 
  distinct(year, mat_s, sd_year, mean_year) %>% 
  filter(year == 1998|year == 2000) 

# known slope
edwards_reg_slope = -0.001
# get se of slope from eq S.1 https://www.int-res.com/articles/suppl/m636p019_supp1.pdf
edwards_reg_sd = (0.0027 - (-0.0047))

bayes_mle_regression_table = as_draws_df(fit_ibts_hierarchical) %>% 
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
# Figure 1 ----------------------------------------------------------------

var_int_groups = fit_varint %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("alpha_raw")) %>% 
  mutate(value = a + sigma_year*value) %>% 
  group_by(name) %>% 
  median_qi(value) %>% 
  mutate(group = as.integer(parse_number(name)),
         model = "Varying Intercepts",
         parameter = "\u03BB") %>% 
  left_join(varint_sims %>% select(group, b) %>% rename(true_value = b))

single_and_varint_plot = bind_rows(var_int_groups, single_mod_summary) %>% # single mode summary is from code for Table 1
  mutate(facet = "a)") %>% 
  ggplot(aes(x = true_value, y = value, color = model, shape = model)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper),
                  position = position_dodge(width = 0.04),
                  size = 0.2) +
  geom_abline() +
  scale_color_grey(end = 0.6) + 
  theme_default() + 
  facet_wrap(~facet) +
  labs(y = "Modeled \u03bb",
       x = "True \u03bb",
       color = "Models",
       shape = "Models") + 
  ylim(-2.3, -1) + 
  xlim(-2.3, -1) +
  theme(legend.title = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.position = c(0.25, 0.9),
        legend.text = element_text(size = 8))

single_and_varint_plot

ggview(single_and_varint_plot, width = 5, height = 5, units = "in")
ggsave(single_and_varint_plot, width = 5, height = 5, units = "in",
       file = "plots/single_and_varint_plot.jpg")
saveRDS(single_and_varint_plot, file = "plots/single_and_varint_plot.rds")


isd_single_samples_plot_long = readRDS(file = "plots/isd_single_samples_plot_long.rds")

fig_1_4panel = single_and_varint_plot + isd_single_samples_plot_long + 
  plot_layout(widths = c(1.8, 1))

ggview(fig_1_4panel, width = 6.5, height = 4, units = "in")
ggsave(fig_1_4panel, file = "plots/fig_1_4panel.jpg", width = 6.5, height = 4, units = "in")
ggsave(fig_1_4panel, file = "plots/fig_1_4panel.pdf", width = 6.5, height = 4, units = "in")
saveRDS(fig_1_4panel, file = "plots/fig_1_4panel.rds")

# Figure 2 ----------------------------------------------------------------
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

# sample posterior
fit_reg_posts = list()

for(i in 1:length(fit_varint_regression)){
  fit_reg_posts[[i]] = brms::as_draws_df(fit_varint_regression[[i]]) %>% 
    mutate(model_run = i)
}

var_reg_summary = bind_rows(fit_reg_posts) %>% 
  select(-lp__, -contains("alpha_raw")) %>% 
  pivot_longer(cols = c(-.chain, -.iteration, -.draw, -model_run)) %>%
  mutate(model = "3) Varying Intercepts Regression") %>% 
  left_join(true_reg_values) %>% 
  mutate(bias = value - true_value) %>% 
  left_join(rhats %>% ungroup %>% distinct(converged, model_run)) %>%  # limit to converged models only
  filter(converged == "yes")   

# summarize
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


# make the plot
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

# ggview(plot_linear_model_bias, width = 6.5, height = 3.5, units = "in")
ggsave(plot_linear_model_bias, width = 6.5, height = 3.5, units = "in",
       file = "plots/fig_2_plot_linear_model_bias.jpg", dpi = 600)
saveRDS(plot_linear_model_bias, file = "plots/plot_linear_model_bias.rds")


# Figure 3 ----------------------------------------------------------------
b_lines = tibble(b = c(-2, -1.6, -1.2),
                 b_known = c(-2, -1.6, -1.2),
                 sd = NA) %>% 
  pivot_longer(cols = -b_known)

sample_size_plot = sample_size_sims %>% 
  mutate(b_known = rep(rep(c(-2, -1.6, -1.2),11),10)) %>% 
  pivot_longer(cols = c(b, sd)) %>% 
  filter(name == "b") %>% 
  ggplot(aes(x = n, y = value)) + 
  facet_wrap(~b_known, scales = "free", ncol = 3) +
  geom_point(size = 0.1) +
  labs(x = "Number of individual body sizes in a sample",
       y = "lambda") + 
  scale_x_log10() +
  geom_hline(data = b_lines, aes(yintercept = value)) + 
  theme_default()

ggview(sample_size_plot, width = 6, height = 2, units = "in")
ggsave(sample_size_plot, width = 6, height = 2, units = "in", 
       file = "plots/fig_3_sample_size_plot.jpg")
saveRDS(sample_size_plot, file = "plots/sample_size_plot.rds")


# Figure 4 ----------------------------------------------------------------
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

# sample posteriors
ibts_posts_hierarchical = as_draws_df(fit_ibts_hierarchical) 

ibts_conds_hierarchical = ibts_posts_hierarchical %>% 
  expand_grid(mat_s = seq(min(count_sims_thin$mat_s), max(count_sims_thin$mat_s), length.out = 20)) %>% 
  mutate(mean_year = unique(count_sims_thin$mean_year),
         sd_year = unique(count_sims_thin$sd_year),
         year = mat_s*sd_year + mean_year) %>% 
  mutate(lambda = a + beta_mat*mat_s)

saveRDS(ibts_conds_hierarchical, file = "posteriors/ibts_conds_hierarchical.rds")

ibts_year_posts_hierarchical = ibts_posts_hierarchical %>% 
  pivot_longer(cols = contains("alpha_raw_year")) %>% 
  mutate(group = parse_number(name)) %>% 
  left_join(count_sims_thin %>% ungroup %>% distinct(group, mat_s, year)) %>% 
  mutate(lambda = a + beta_mat*mat_s + sigma_year*value) %>% 
  group_by(group, mat_s, year) %>% 
  median_qi(lambda)




ibts_posts_norand = as_draws_df(fit_ibts_norand_norand) 

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



# make plot
bayes_mle_regression_lines = ibts_conds_hierarchical %>% 
  group_by(mat_s, year) %>% 
  select(!contains("alpha_raw")) %>%
  mean_qi(lambda) %>% 
  mutate(method = "b) Bayesian - hierarchical") %>% 
  bind_rows(edwards_mle_estimates)

bayes_mle_year_summaries = edwards_mle_estimates %>% 
  mutate(method = "a) MLE - two steps") %>% 
  bind_rows(ibts_year_posts_hierarchical %>% mutate(method = "b) Bayesian - hierarchical"))

bayes_mle_regression_plot = bayes_mle_regression_lines %>%
  ggplot(aes(x = year, y = lambda, ymin = .lower, ymax = .upper)) + 
  geom_ribbon(data = . %>% filter(grepl("Bayes", method)), alpha = 0.2) +
  geom_line(data = . %>% filter(grepl("Bayes", method)),
            color = "black", linewidth = 0.3) +
  geom_pointrange(data = bayes_mle_year_summaries,
                  size = 0.005) +
  geom_smooth(data = . %>% filter(!grepl("Bayes", method)),
              method = "lm", color = "black", linewidth = 0.3) +
  ylim(-2, -1.3) +
  facet_wrap(~method) +
  theme_default() +
  labs(y = "\u03bb",
       x = "Year") +
  theme(strip.text.x = element_text(hjust = 0.05)) +
  NULL

saveRDS(bayes_mle_regression_plot, file = "plots/bayes_mle_regression_plot.rds")
ggview(bayes_mle_regression_plot, width = 6, height = 3, units = "in")
ggsave(bayes_mle_regression_plot, width = 6, height = 3, units = "in", file = "plots/fig_4_bayes_mle_regression_plot.jpg")

