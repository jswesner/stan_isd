library(tidyverse)
library(tidybayes)
library(brms)
library(ggthemes)
source("code/sandbox/save_plot_and_data.R")

#1) read in model summaries
separate_lambda_summaries = readRDS(file = "models/coverage_models/all_lambda_summaries.rds") %>% 
  mutate(model = "separate_models", pars = "lambda") %>% rename(true_value = true_lambda)
var_lambda_summaries = readRDS("posteriors/var_summaries.rds") %>% mutate(pars = "lambda", model = "varying_intercepts")
fixed_lambda_summaries = readRDS("posteriors/fixed_summaries.rds") %>% mutate(pars = "lambda")

all_lambda_summaries = bind_rows(separate_lambda_summaries,
                                 var_lambda_summaries,
                                 fixed_lambda_summaries) 


#2) wrangle model summaries
plot_data = all_lambda_summaries %>% 
  filter(xmax == 1000) %>% # limit to xmax = 1000
  mutate(cov95 = case_when(true_value > `2.5%` & true_value <`97.5%` ~ "yes", TRUE ~ "no"),
         cov50 = case_when(true_value > `25%` & true_value <`75%` ~ "yes", TRUE ~ "no"),
         cov99 = case_when(true_value > `0.5%` & true_value <`99.5%` ~ "yes", TRUE ~ "no")) %>%
  group_by(xmax, model, pars, true_value) %>% 
  # filter(replicate <= 100) %>% 
  mutate(model = case_when(model == "separate_models" ~ "a) Separate Models",
                           model == "fixed_predictors" ~ "b) Fixed Predictor",
                           TRUE ~ "c) Varying Intercepts"))

#3) plot model summaries
coverage_plot = plot_data %>% 
  ggplot(aes(x = replicate, y = mean, ymin = `2.5%`, ymax = `97.5%`,
             color = as.factor(cov95))) + 
  geom_linerange(linewidth = 0.3) + 
  facet_grid(true_value ~ model) +
  labs(y = "\u03bb",
       x = "Simulation",
       color = "CrI contains true value?") +
  ggthemes::scale_color_colorblind() +
  theme_default() +
  geom_hline(data = plot_data %>% ungroup %>% distinct(model, cov95, true_value),
             aes(yintercept = true_value)) +
  theme(strip.text.y = element_blank(),
        legend.position = "top")

# coverage_plot = readRDS("plots/coverage_plot.rds")
ggview::ggview(coverage_plot, width = 6, height = 7)

save_plot_and_data(coverage_plot, file_name = "plots/coverage_plot", width = 6, height = 7)
save_plot_and_data(coverage_plot, file_name = "ms/coverage_plot", width = 6, height = 7)

