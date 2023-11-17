library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
library(sizeSpectra)
library(janitor)
library(isdbayes)
rstan_options("auto_write" = TRUE)

theme_set(theme_default())

###
# LOAD PLOT #
hierarchical_regression_plot = readRDS(file = "plots/hierarchical_regression_plot.rds")
###


set.seed = 123123
intercept = -1.5
# var_int_sd = 0.1
beta = -0.1
n_ind = 300
reps = 12
  
lambda_sims = tibble(group = (1:1)) %>% 
    expand_grid(predictor = -seq(-2, 2, length.out = 12)) %>% 
    # expand_grid(replicates = 1:2) %>% 
    mutate(b = intercept + predictor*beta,
           id = 1:nrow(.)) %>% 
  add_row(group = 1,
          predictor = 2.5, b = -1.1, id = 7) 
  
saveRDS(lambda_sims, file = "data/lambda_sims.rds")

sim_data = lambda_sims %>% 
  expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = 1, xmax = 1000,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(group, id, xmin, xmax, predictor) %>% 
    add_count(x, name = "counts") %>% 
  mutate(remove = case_when(predictor == 2.5 & individual > 50 ~ "remove",
                            TRUE ~ "keep")) %>% 
  filter(remove == "keep") %>%
    ungroup %>% 
  mutate(group = as.factor(predictor))


# Estimate lambdas --------------------------------------------------------

sim_data_list = sim_data %>% group_by(group) %>% group_split()

brm_single_mods = brm_multiple(x | vreal(counts, xmin, xmax) ~ 1,
                               data = sim_data_list,
                               stanvars = stanvars,
                               family = paretocounts(),
                               # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
                               chains = 1, iter = 1000,
                               cores = 4,
                               combine = F)

saveRDS(brm_single_mods, file = "models/brm_single_mods.rds")
brm_single_mods = readRDS("models/brm_single_mods.rds")


brm_single_fixefs = NULL

for(i in 1:length(brm_single_mods)){
  brm_single_fixefs[[i]] = as_tibble(fixef(brm_single_mods[[i]])) %>% 
    mutate(group = i)
}

single_lambdas = bind_rows(brm_single_fixefs) %>% 
  mutate(predictor = sort(unique(sim_data$predictor)))

saveRDS(single_lambdas, file = "data/single_lambdas.rds")

# fit models ---------------------------------------------------------------
# pre load fitted model to use the update function for updating
brm_regress_prior = readRDS(file = 'models/brm_regress_prior.rds')
single_lambdas = readRDS(file = "data/single_lambdas.rds")


brm_twostep_regress = brm(Estimate ~ predictor,
                          data = single_lambdas,
                          family = gaussian(),
                          iter = 1000,
                          chains = 4,
                          file = "models/brm_twostep_regress.rds",
                          file_refit = "on_change")

brm_twostep_regress_rand = update(brm_twostep_regress, 
                                  formula = .~ predictor + (1|group),
                                  newdata = single_lambdas,
                                  file = "models/brm_twostep_regress_rand.rds",
                                  file_refit = "on_change")
# 
# brm_twostep_regress_rand = brm(Estimate ~ predictor + (1|group), 
#                           data = single_lambdas,
#                           family = gaussian(),
#                           iter = 1000, 
#                           chains = 4,
#                           file = "models/brm_twostep_regress_rand.rds",
#                           file_refit = "on_change")
# 
# brm_twostep_regress_mi = brm(Estimate|mi(Est.Error) ~ predictor,
#                              data = single_lambdas,
#                              iter = 1000, 
#                              chains = 4,
#                              file = "models/brm_twostep_regress_mi.rds",
#                              file_refit = "on_change")
# 
# brm_regress = brm(x | vreal(counts, xmin, xmax) ~ predictor, 
#                                        data = sim_data,
#                                        stanvars = stanvars,
#                                        family = paretocounts(),
#                                        # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
#                                        chains = 1, iter = 1000,
#                                        cores = 4,
#                   file = 'models/brm_regress.rds',
#                   file_refit = "on_change")
# 
# brm_regress_prior = update(brm_regress, 
#                            prior = c(prior(normal(-0.1, 0.05), class = "b")),
#                            file = 'models/brm_regress_prior.rds',
#                            file_refit = "on_change")
# 
brm_regress_weakprior_rand = update(brm_regress_prior,
                                    newdata = sim_data,
                                    prior = c(prior(normal(0, 0.5), class = "b"),
                                              prior(normal(-1.5, 1), class = "Intercept"),
                                              prior(exponential(1), class = "sd")),
                                    formula = . ~ predictor + (1|group),
                                    iter = 1000,
                                    chains = 3,
                                    file_refit = "on_change",
                                    file = 'models/brm_regress_weakprior_rand.rds')
#
brm_regress_prior_rand = update(brm_regress_prior,
                                newdata = sim_data,
                                prior = c(prior(normal(-0.1, 0.02), class = "b"),
                                          prior(normal(-1.5, 0.1), class = "Intercept"),
                                          prior(exponential(1), class = "sd")),
                                formula = . ~ predictor + (1|group),
                                iter = 1000,
                                chains = 1,
                                file_refit = "on_change",
                                file = 'models/brm_regress_prior_rand.rds')

# get posteriors ----------------------------------------------------------
brm_twostep_regress = readRDS("models/brm_twostep_regress.rds")
brm_twostep_regress_mi = readRDS("models/brm_twostep_regress_mi.rds")
brm_regress = readRDS("models/brm_regress.rds")
brm_regress_prior = readRDS("models/brm_regress_prior.rds")
brm_regress_weakprior_rand = readRDS("models/brm_regress_weakprior_rand.rds")
brm_regress_prior_rand = readRDS("models/brm_regress_prior_rand.rds")

a = brm_twostep_regress
b = brm_twostep_regress_mi
c = brm_regress
d = brm_regress_prior
e = brm_regress_weakprior_rand
f = brm_regress_prior_rand

mod_list = list(a,b,c,d,e,f)

mod_draws = NULL

for(i in 1:length(mod_list)) {
  mod_draws[[i]] = tibble(predictor = seq(min(mod_list[[i]]$data$predictor),
                                     max(mod_list[[i]]$data$predictor),
                                     length.out = 20)) %>% 
    mutate(counts = 1, # placeholders
           xmin = 1, 
           xmax = 1000,
           Est.Error = 1,
           group = 1) %>% 
    add_epred_draws(mod_list[[i]], re_formula = NA) %>% 
    mutate(model = i)
}

weak_rand_samples = mod_list[[5]]$data %>% 
  distinct(predictor, group) %>% 
  mutate(counts = 1, xmin = 1, xmax = 1000) %>% 
  add_epred_draws(mod_list[[5]], re_formula = NULL) %>% 
  mutate(model = 5,
         model_name = "b) Bayesian - One Step\nHierarchical + Weak Priors")

strong_rand_samples = mod_list[[6]]$data %>% 
  distinct(predictor, group) %>% 
  mutate(counts = 1, xmin = 1, xmax = 1000) %>% 
  add_epred_draws(mod_list[[6]], re_formula = NULL) %>% 
  mutate(model = 6,
         model_name = "c) Bayesian - One Step\nHierarchical + Strong Priors")

rand_samples = bind_rows(weak_rand_samples, strong_rand_samples)

samples_all = rand_samples %>% 
  group_by(model, model_name, predictor) %>% 
  median_qi(.epred) 

# compare slopes ----------------------------------------------------------

pars = lapply(mod_list, fixef)
pars_id = NULL
for(i in 1:length(pars)) {
  pars_id[[i]] = pars[[i]] %>% as_tibble() %>% mutate(model = i,
                                                      par = c("Intercept", "Slope"))
  pars_id_out = bind_rows(pars_id)
}

pars_fixefs = pars_id_out %>% left_join(tibble(par = c("Intercept", "Slope"),
                                               true_value = c(mean(single_lambdas$Estimate), beta))) %>% 
  left_join(all_posts %>% ungroup %>% distinct(model, model_name)) %>% 
  filter(!is.na(model_name)) %>% 
  mutate(order = as.numeric(as.factor(str_sub(model_name, 1, 1))))

pars_fixefs %>% 
  mutate(par = case_when(par == "Intercept" ~ "a) Intercept", 
                         TRUE ~ "b) Slope")) %>% 
  ggplot(aes(y = reorder(model_name, -order), x = Estimate, xmin = Q2.5, xmax = Q97.5)) + 
  geom_pointrange() +
  facet_wrap(~par, scales = "free_x") +
  geom_vline(aes(xintercept = true_value), linetype = "dashed") +
  labs(y = "",
       x = "Parameter Estimate") +
  theme(strip.text.x = element_text(hjust = 0))

# plot posteriors ---------------------------------------------------------

all_posts = bind_rows(mod_draws) %>% 
  mutate(model_name = case_when(model == 1 ~ "a) MLE - Two Steps\nSimple Linear Regression",
                                # model == 2 ~ "b) Two Steps\nErrors in Variables Regression",
                                # model == 3 ~ "c) One Step\nWeak Priors",
                                # model == 4 ~ "d) One Step\nStrong Priors",
                                model == 5 ~ "b) Bayesian - One Step\nHierarchical + Weak Priors",
                                model == 6 ~ "c) Bayesian - One Step\nHierarchical + Strong Priors")) %>%
  filter(model %in% c(1,5,6)) 

all_samples = samples_all %>% 
  mutate(dots = "Hierarchical Estimated \u03BB's") %>% 
  bind_rows(single_lambdas %>% 
              # expand_grid(model = 1:6) %>% 
              mutate(model = 1) %>% 
              rename(.epred = Estimate) %>% 
              mutate(dots = "Individually Estimated \u03BB's") %>% 
              left_join(all_posts %>% ungroup %>% distinct(model, model_name))) %>% 
  filter(!is.na(model_name)) %>% 
  mutate(dot_size = case_when(predictor == 2.5 ~ 20,
                              TRUE ~ 400),
         outlier = case_when(predictor == 2.5 ~ "Outlier", 
                             TRUE ~ "Not outlier")) 


slope_text = pars_fixefs  %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  mutate(text = paste0(par, ": ", Estimate, " (", Q2.5, " to ", Q97.5, ")")) %>% 
  filter(par == "Slope")


hierarchical_regression_plot = all_posts %>% 
  ggplot(aes(x = predictor, y = .epred)) + 
  stat_lineribbon(.width = c(0.5, 0.95), color = "black",
                  fill = "grey20", alpha = 0.4, 
                  linewidth = 0.2) + 
  facet_wrap(~model_name, nrow = 1) +
  geom_pointrange(data = all_samples, aes(y = .epred, x = predictor,
                                          ymin = .lower,
                                          ymax = .upper,
                                          shape = outlier,
                                          alpha = dots),
                  size = 0.2,
                  linewidth = 0.2) +
  facet_wrap(~model_name) +
  scale_alpha_manual(values = c(0.4, .8)) +
  theme(strip.text.x = element_text(angle = 0, hjust = 0),
        legend.title = element_blank(),
        legend.position = c(0.15, 0.9),
        legend.text = element_text(size = 9),
        text = element_text(size = 9)) + 
  labs(y = "\u03BB",
       x = "Predictor") +
  geom_abline(intercept = mean(lambda_sims$b), slope = -0.1) +
  coord_cartesian(ylim = c(-1.8, -1)) +
  guides(shape = "none",
         alpha = "none") + 
  geom_text(data = slope_text, aes(label = text),
            x = -0.5, y = -1.7, size = 2)


ggview::ggview(hierarchical_regression_plot, 
               width = 6.5, height = 2.5, units = "in")

saveRDS(hierarchical_regression_plot, file = "plots/hierarchical_regression_plot.rds")
ggsave(hierarchical_regression_plot, file = "plots/hierarchical_regression_plot.jpg",
       width = 6.5, height = 3, units = "in", dpi = 500)



# posterior processing ----------------------------------------------------

as_draws_df(a) %>% 
  reframe(slope_test = sum(b_predictor < 0)/nrow(.))
as_draws_df(e) %>% 
  reframe(slope_test = sum(b_predictor < 0)/nrow(.))
as_draws_df(f) %>% 
  reframe(slope_test = sum(b_predictor < 0)/nrow(.))
