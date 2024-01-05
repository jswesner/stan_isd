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


# simulate data
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
    mutate(group = i, 
           reps = nrow(brm_single_mods[[i]]$data))
}

single_lambdas = bind_rows(brm_single_fixefs) %>% 
  mutate(predictor = sort(unique(sim_data$predictor)))

saveRDS(single_lambdas, file = "data/single_lambdas.rds")

# fit models ---------------------------------------------------------------
# pre load fitted model to use the update function for updating
brm_regress_prior = readRDS(file = 'models/brm_regress_prior.rds')
single_lambdas = readRDS(file = "data/single_lambdas.rds")
brm_twostep_regress = readRDS("models/brm_twostep_regress.rds")

brm_twostep_regress = brm(Estimate ~ predictor,
                          data = single_lambdas,
                          family = gaussian(),
                          iter = 1000,
                          chains = 4,
                          file = "models/brm_twostep_regress.rds",
                          file_refit = "on_change")

brm_twostep_regress = update(brm_twostep_regress, iter = 2000, chains = 4)

saveRDS(brm_twostep_regress, file = "models/brm_twostep_regress.rds")

brm_twostep_regress_rand = update(brm_twostep_regress, 
                                  formula = .~ predictor + (1|group),
                                  newdata = single_lambdas,
                                  file = "models/brm_twostep_regress_rand.rds",
                                  file_refit = "on_change", 
                                  iter = 2000, 
                                  chains = 4)
# 
# brm_twostep_regress_rand = brm(Estimate ~ predictor + (1|group), 
#                           data = single_lambdas,
#                           family = gaussian(),
#                           iter = 1000, 
#                           chains = 4,
#                           file = "models/brm_twostep_regress_rand.rds",
#                           file_refit = "on_change")
# 
brm_twostep_regress_mi = brm(Estimate|mi(Est.Error) ~ predictor,
                             data = single_lambdas,
                             iter = 2000,
                             chains = 4,
                             file = "models/brm_twostep_regress_mi.rds",
                             file_refit = "on_change")
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
                                    iter = 2000,
                                    chains = 4,
                                    file_refit = "on_change",
                                    file = 'models/brm_regress_weakprior_rand.rds')
#
brm_regress_prior_rand = update(brm_regress_prior,
                                newdata = sim_data,
                                prior = c(prior(normal(-0.1, 0.02), class = "b"),
                                          prior(normal(-1.5, 0.1), class = "Intercept"),
                                          prior(exponential(1), class = "sd")),
                                formula = . ~ predictor + (1|group),
                                iter = 2000,
                                chains = 4,
                                file_refit = "on_change",
                                file = 'models/brm_regress_prior_rand.rds')


# posterior processing ----------------------------------------------------

as_draws_df(a) %>% 
  reframe(slope_test = sum(b_predictor < 0)/nrow(.))
as_draws_df(e) %>% 
  reframe(slope_test = sum(b_predictor < 0)/nrow(.))
as_draws_df(f) %>% 
  reframe(slope_test = sum(b_predictor < 0)/nrow(.))


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


# weighted regression  ----------------------------------------------------


