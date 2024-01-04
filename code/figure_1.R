library(tidyverse)
library(tidybayes)
library(brms)
library(ggthemes)
source("code/sandbox/save_plot_and_data.R")

#1) read in model summaries
# separate_lambda_summaries = readRDS(file = "models/coverage_models/all_lambda_summaries.rds") %>% 
separate_lambda_summaries = readRDS(file = "posteriors/separate_summaries.rds") 
var_lambda_summaries = readRDS("posteriors/var_summaries.rds") %>% mutate(pars = "lambda", model = "varying_intercepts")
fixed_lambda_summaries = readRDS("posteriors/fixed_summaries.rds") %>% mutate(pars = "lambda")

all_lambda_summaries = bind_rows(separate_lambda_summaries,
                                 var_lambda_summaries,
                                 fixed_lambda_summaries) 


#2) wrangle model summaries
plot_data = all_lambda_summaries %>% 
  # filter(xmax == 1000) %>% # limit to xmax = 1000
  mutate(cov95 = case_when(true_value > `2.5%` & true_value <`97.5%` ~ "yes", TRUE ~ "no"),
         cov50 = case_when(true_value > `25%` & true_value <`75%` ~ "yes", TRUE ~ "no"),
         cov99 = case_when(true_value > `0.5%` & true_value <`99.5%` ~ "yes", TRUE ~ "no")) %>%
  group_by(xmax, model, pars, true_value) %>% 
  # filter(replicate <= 100) %>% 
  mutate(model = case_when(model == "separate_models" ~ "a) Separate Models",
                           model == "fixed_predictors" ~ "b) Fixed Predictor",
                           TRUE ~ "c) Varying Intercepts"))

#3) plot model summaries
# make labels
labels = plot_data %>% 
  ungroup %>% 
  distinct(model, true_value)

# make plot
coverage_plot = plot_data %>% 
  arrange(model, true_value, `2.5%`) %>% 
  ungroup %>% 
  mutate(replicate = rep(1:1000, 21)) %>% 
  ggplot(aes(y = replicate)) + 
  geom_linerange(linewidth = 0.3,
                 aes(x = mean, xmin = `2.5%`, xmax = `97.5%`,
                     color = as.factor(cov95))) + 
  # facet_grid(true_value ~ model, scales = "free") +
  lemon::facet_rep_grid(true_value ~ model, scales = "free",
                        space = "free") +
  labs(x = "\u03bb",
       y = "Simulation",
       color = "CrI contains true value?") +
  ggthemes::scale_color_colorblind() +
  geom_text(data = labels, aes(label = true_value, x = true_value - 0.2, y = 1100),
            size = 3) +
  theme_default() +
  geom_vline(data = plot_data %>% ungroup %>% distinct(model, cov95, true_value),
             aes(xintercept = true_value)) +
  theme(legend.position = "top",
        strip.text.x = element_text(hjust = 0,
                                    vjust = 3),
        strip.text.y = element_blank()) +
  scale_y_continuous(breaks = c(1, 500, 1000),
  ) +
  coord_cartesian(clip = "off",
                  ylim = c(1, 1100)) +
  NULL


# coverage_plot = readRDS("plots/coverage_plot.rds")
ggview::ggview(coverage_plot, width = 6.5, height = 9)

save_plot_and_data(coverage_plot, file_name = "plots/coverage_plot", width = 6.5, height = 9)
save_plot_and_data(coverage_plot, file_name = "ms/coverage_plot", width = 6.5, height = 9)

