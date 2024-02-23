library(brms)
library(tidyverse)
library(isdbayes)
library(tidybayes)
library(scales)


# Figure S5 Prior Sensitivity

# Fit models --------------------------------------------------------------

#1) Function to simulate data
make_data = function(n_ind = 300,
                     xmin = 1, 
                     xmax = 1000,
                     b = -1.6){
  
  tibble(xmin = xmin,
         xmax = xmax,
         b = b) %>% 
    expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = xmin, 
           xmax = xmax,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    mutate(counts = 1) %>% 
    ungroup
}

#2) Precompile dummy model (for updating later)
reg_dat = make_data()

# precompile dummy model
reg_fit = brm(x | vreal(counts, xmin, xmax) ~ 1, 
              data = reg_dat,
              stanvars = stanvars,
              family = paretocounts(),
              # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
              chains = 1, iter = 10,
              cores = 4)



#3) Update dummy model with different priors

fit_prior_a_temp = update(reg_fit, prior = c(prior(normal(-1.5, 5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_a = update(fit_prior_a_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_b_temp = update(reg_fit, prior = c(prior(normal(-1.5, 2), class = "Intercept")), chains = 1, iter = 10)
fit_prior_b = update(fit_prior_b_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_c_temp = update(reg_fit, prior = c(prior(normal(-1.5, 1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_c = update(fit_prior_c_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_d_temp = update(reg_fit, prior = c(prior(normal(-1.5, 0.5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_d = update(fit_prior_d_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_e_temp = update(reg_fit, prior = c(prior(normal(-1.5, 0.1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_e = update(fit_prior_e_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_f_temp = update(reg_fit, prior = c(prior(normal(-1.5, 0.01), class = "Intercept")), chains = 1, iter = 10)
fit_prior_f = update(fit_prior_f_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_g_temp = update(reg_fit, prior = c(prior(normal(-1.5, 0.001), class = "Intercept")), chains = 1, iter = 10)
fit_prior_g = update(fit_prior_g_temp, newdata = make_data(), chains = 2, iter = 1000)


#4) Combine models
fit_priors = list(fit_prior_a,
                  fit_prior_b,
                  fit_prior_c,
                  fit_prior_d,
                  fit_prior_e,
                  fit_prior_f,
                  fit_prior_g)

saveRDS(fit_priors, file = "models/figs5_mods_1.5.rds")

priors_to_add = tibble(fit = 1:length(fit_priors),
                       prior_mean = -1.5,
                       prior_sd = c(10, 5, 2, 1, 0.1, 0.01, 0.001))

fit_prior_temp = NULL

for(i in 1:length(fit_priors)) {
  fit_prior_temp[[i]] = as_draws_df(fit_priors[[i]]) %>% 
    mutate(fit = i) %>% 
    left_join(priors_to_add)
}

#5) repeat with different prior means

# -2
fit_prior_a_temp = update(reg_fit, prior = c(prior(normal(-2, 5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_a = update(fit_prior_a_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_b_temp = update(reg_fit, prior = c(prior(normal(-2, 2), class = "Intercept")), chains = 1, iter = 10)
fit_prior_b = update(fit_prior_b_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_c_temp = update(reg_fit, prior = c(prior(normal(-2, 1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_c = update(fit_prior_c_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_d_temp = update(reg_fit, prior = c(prior(normal(-2, 0.5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_d = update(fit_prior_d_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_e_temp = update(reg_fit, prior = c(prior(normal(-2, 0.1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_e = update(fit_prior_e_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_f_temp = update(reg_fit, prior = c(prior(normal(-2, 0.01), class = "Intercept")), chains = 1, iter = 10)
fit_prior_f = update(fit_prior_f_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_g_temp = update(reg_fit, prior = c(prior(normal(-2, 0.001), class = "Intercept")), chains = 1, iter = 10)
fit_prior_g = update(fit_prior_g_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_priors = list(fit_prior_a,
                  fit_prior_b,
                  fit_prior_c,
                  fit_prior_d,
                  fit_prior_e,
                  fit_prior_f,
                  fit_prior_g)

saveRDS(fit_priors, file = "models/figs5_mods_2.rds")

# -1.8
fit_prior_a_temp = update(reg_fit, prior = c(prior(normal(-1.8, 5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_a = update(fit_prior_a_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_b_temp = update(reg_fit, prior = c(prior(normal(-1.8, 2), class = "Intercept")), chains = 1, iter = 10)
fit_prior_b = update(fit_prior_b_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_c_temp = update(reg_fit, prior = c(prior(normal(-1.8, 1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_c = update(fit_prior_c_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_d_temp = update(reg_fit, prior = c(prior(normal(-1.8, 0.5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_d = update(fit_prior_d_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_e_temp = update(reg_fit, prior = c(prior(normal(-1.8, 0.1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_e = update(fit_prior_e_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_f_temp = update(reg_fit, prior = c(prior(normal(-1.8, 0.01), class = "Intercept")), chains = 1, iter = 10)
fit_prior_f = update(fit_prior_f_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_g_temp = update(reg_fit, prior = c(prior(normal(-1.8, 0.001), class = "Intercept")), chains = 1, iter = 10)
fit_prior_g = update(fit_prior_g_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_priors = list(fit_prior_a,
                  fit_prior_b,
                  fit_prior_c,
                  fit_prior_d,
                  fit_prior_e,
                  fit_prior_f,
                  fit_prior_g)

saveRDS(fit_priors, file = "models/figs5_mods_1.8.rds")


# -1.2
fit_prior_a_temp = update(reg_fit, prior = c(prior(normal(-1.2, 5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_a = update(fit_prior_a_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_b_temp = update(reg_fit, prior = c(prior(normal(-1.2, 2), class = "Intercept")), chains = 1, iter = 10)
fit_prior_b = update(fit_prior_b_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_c_temp = update(reg_fit, prior = c(prior(normal(-1.2, 1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_c = update(fit_prior_c_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_d_temp = update(reg_fit, prior = c(prior(normal(-1.2, 0.5), class = "Intercept")), chains = 1, iter = 10)
fit_prior_d = update(fit_prior_d_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_e_temp = update(reg_fit, prior = c(prior(normal(-1.2, 0.1), class = "Intercept")), chains = 1, iter = 10)
fit_prior_e = update(fit_prior_e_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_f_temp = update(reg_fit, prior = c(prior(normal(-1.2, 0.01), class = "Intercept")), chains = 1, iter = 10)
fit_prior_f = update(fit_prior_f_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_prior_g_temp = update(reg_fit, prior = c(prior(normal(-1.2, 0.001), class = "Intercept")), chains = 1, iter = 10)
fit_prior_g = update(fit_prior_g_temp, newdata = make_data(), chains = 2, iter = 1000)

fit_priors = list(fit_prior_a,
                  fit_prior_b,
                  fit_prior_c,
                  fit_prior_d,
                  fit_prior_e,
                  fit_prior_f,
                  fit_prior_g)

saveRDS(fit_priors, file = "models/figs5_mods_1.2.rds")




