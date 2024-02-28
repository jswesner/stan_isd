library(tidyverse)
library(rstan)
library(brms)
library(scales)
library(tidybayes)
library(ggthemes)
library(sizeSpectra)
library(janitor)
library(isdbayes)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.

theme_set(theme_default())

#1) load fitted models 
a = readRDS("models/fig4a_mod.rds")
b = readRDS("models/fig4b_mod.rds")
c = readRDS("models/fig4c_mod.rds")
d = readRDS("models/fig4d_mod.rds")
single_lambdas = readRDS(file = "data/single_lambdas.rds") 
lambda_sims = readRDS(file = "data/lambda_sims.rds") 

#2) Extract posteriors of the regression lines from fitted models
mod_list = list(a,b,c,d)

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

#3) Extract posteriors of individual lambdas from the varying intercepts models 
weak_rand_samples = mod_list[[3]]$data %>% 
  distinct(predictor, group) %>% 
  mutate(counts = 1, xmin = 1, xmax = 1000) %>% 
  add_epred_draws(mod_list[[3]], re_formula = NULL) %>% 
  mutate(model = 3,
         model_name = "b) Bayesian - One Step\nHierarchical + Weak Priors")

strong_rand_samples = mod_list[[4]]$data %>% 
  distinct(predictor, group) %>% 
  mutate(counts = 1, xmin = 1, xmax = 1000) %>% 
  add_epred_draws(mod_list[[4]], re_formula = NULL) %>% 
  mutate(model = 4,
         model_name = "c) Bayesian - One Step\nHierarchical + Strong Priors")


rand_samples = bind_rows(weak_rand_samples, strong_rand_samples)

#4) summarize lambdas from the varying intercepts models 
samples_all = rand_samples %>% 
  group_by(model, model_name, predictor) %>% 
  median_qi(.epred) 

#5) Summarize fixed effects (to add text to the final plot)
pars = lapply(mod_list, fixef)
pars_id = NULL
for(i in 1:length(pars)) {
  pars_id[[i]] = pars[[i]] %>% as_tibble() %>% mutate(model = i,
                                                      par = c("Intercept", "Slope"))
  pars_id_out = bind_rows(pars_id)
}


all_posts = bind_rows(mod_draws) %>% 
  mutate(model_name = case_when(model == 1 ~ "a) MLE - Two Steps\nSimple Linear Regression",
                                model == 2 ~ "b) Two Steps\nRegression with measurement error",
                                # model == 3 ~ "c) One Step\nWeak Priors",
                                # model == 4 ~ "d) One Step\nStrong Priors",
                                model == 3 ~ "c) Bayesian - One Step\nHierarchical + Weak Priors",
                                model == 4 ~ "d) Bayesian - One Step\nHierarchical + Strong Priors")) 


pars_fixefs = pars_id_out %>% left_join(tibble(par = c("Intercept", "Slope"),
                                               true_value = c(mean(single_lambdas$Estimate), beta))) %>% 
  left_join(all_posts %>% ungroup %>% distinct(model, model_name)) %>% 
  filter(!is.na(model_name)) %>% 
  mutate(order = as.numeric(as.factor(str_sub(model_name, 1, 1))))


#6) Wrangle data for the final plot

all_samples = samples_all %>% 
  mutate(dots = "Hierarchical Estimated \u03BB's") %>% 
  bind_rows(single_lambdas %>% 
              expand_grid(model = 1:2) %>%
              rename(.epred = Estimate) %>% 
              mutate(dots = "Individually Estimated \u03BB's") %>% 
              left_join(all_posts %>% ungroup %>% distinct(model, model_name))) %>% 
  filter(!is.na(model_name)) %>% 
  mutate(dot_size = case_when(predictor == 2.5 ~ 20,
                              TRUE ~ 400),
         outlier = case_when(predictor == 2.5 ~ "Outlier", 
                             TRUE ~ "Not outlier")) %>% 
  mutate(model_name = case_when(model == 1 ~ "a) MLE - Two Steps\nSimple Linear Regression",
                                model == 2 ~ "b) Two Steps\nRegression with measurement error",
                                model == 3 ~ "c) Bayesian - One Step\nHierarchical + Weak Priors",
                                model == 4 ~ "d) Bayesian - One Step\nHierarchical + Strong Priors")) %>%
  mutate(.lower = case_when(model == 2 ~ Q2.5,
                            TRUE ~ .lower),
         .upper = case_when(model == 2 ~ Q97.5,
                            TRUE ~ .upper))

slope_text = pars_fixefs  %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  mutate(text = paste0(par, ": ", Estimate, " (", Q2.5, " to ", Q97.5, ")")) %>% 
  filter(par == "Slope")

#7) Make plot
hierarchical_regression_plot = all_posts %>% 
  ggplot(aes(x = predictor, y = .epred)) + 
  stat_lineribbon(.width = c(0.5, 0.95), color = "gray30",
                  fill = "gray30",
                  alpha = 0.5,
                  linewidth = 0.2) + 
  facet_wrap(~model_name, nrow = 2) +
  geom_pointrange(data = all_samples, aes(y = .epred, x = predictor,
                                          ymin = .lower,
                                          ymax = .upper,
                                          shape = outlier),
                  size = 0.2,
                  linewidth = 0.2,
                  color = "gray30") +
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
               width = 6.5, height = 6.5, units = "in")

save_plot_and_data(hierarchical_regression_plot, file_name = "plots/fig4hierarchical_regression_plot",
       width = 6.5, height = 6.5, dpi = 500)
save_plot_and_data(hierarchical_regression_plot, file_name = "ms/Figure4", file_type = "pdf",
                   width = 6.5, height = 6.5, dpi = 500)
