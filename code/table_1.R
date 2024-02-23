library(tidyverse)
library(tidybayes)
source("code/sandbox/save_plot_and_data.R")

#1) read in posterior summaries 
all_lambda_summaries = readRDS(file = "posteriors/fig1abc_posterior_summaries.rds") %>% 
  mutate(true_value  = parse_number(as.character(true_value))) %>% 
  mutate(pars = "lambda")


#2) calculate coverage and bias
coverage = all_lambda_summaries %>% 
  # filter(xmax == 1000) %>% # limit to xmax = 1000
  mutate(cov95 = case_when(true_value > `2.5%` & true_value <`97.5%` ~ "yes", TRUE ~ "no")) %>%
  group_by(true_value, model, cov95) %>% 
  tally() %>%
  mutate(total = sum(n)) %>% 
  filter(cov95 == "yes") %>% 
  mutate(within95 = n/total) %>% 
  select(-n, -total, -cov95)

bias = all_lambda_summaries %>% 
  mutate(bias = `50%` - true_value) %>% 
  group_by(model, true_value) %>% 
  reframe(mean_bias = mean(bias),
          sd_bias = sd(bias))


coverage_and_bias = left_join(coverage,bias) %>% 
  mutate(model = case_when(model == "fixed" ~ "fixed_predictors",
                           model == "separate" ~ "separate_models",
                           TRUE ~ model))

#3) format table and save
table_1 = coverage_and_bias %>%
  mutate(within95 = round(within95, 2),
         mean_bias = round(mean_bias, 3),
         sd_bias = round(sd_bias,2)) %>% 
  mutate(within95 = as.character(within95)) %>% 
  rename(true_lambda = true_value,
         Coverage = within95) %>% 
  select(model, true_lambda, Coverage,
         mean_bias, sd_bias) %>% 
  mutate(Bias = paste0(mean_bias, " (", sd_bias, ")"),
         model = case_when(model == "fixed_predictors" ~ "Fixed Predictors",
                           model == "separate_models" ~ "Separate Models",
                           TRUE ~ "Varying Intercepts"),
         model = as.factor(model),
         model = fct_relevel(model, "Separate Models", "Fixed Predictors", "Varying Intercepts")) %>% 
  arrange(model) %>% 
  select(-mean_bias, -sd_bias) %>% 
  pivot_longer(cols = c(-model, -true_lambda)) %>% 
  pivot_wider(names_from = model, values_from = value) %>% 
  arrange(desc(name))

write_csv(table_1, file = "tables/table_1.csv")
