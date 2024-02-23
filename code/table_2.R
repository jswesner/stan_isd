library(tidyverse)

# 1) load posteriors
fig3_posterior_summaries = readRDS(file = "posteriors/fig3_posterior_summaries.rds")

# 2) calculate coverage and bias
coverage = fig3_posterior_summaries %>% 
  mutate(cov95 = case_when(true_values > .lower & true_values < .upper ~ "yes", TRUE ~ "no")) %>%
  group_by(name, cov95) %>% 
  tally() %>%
  mutate(total = sum(n)) %>% 
  filter(cov95 == "yes") %>% 
  mutate(Coverage = round(n/total, 2)) %>% 
  select(-n, -total, -cov95)


bias = fig3_posterior_summaries %>% 
  mutate(bias = value - true_values) %>% 
  group_by(name) %>% 
  reframe(mean_bias = round(mean(bias), 3),
          sd_bias = round(sd(bias), 2)) %>% 
  mutate(`Bias (mean, sd)` = paste0(mean_bias, " (", sd_bias, ")")) %>% 
  select(-mean_bias, -sd_bias)

table_2 = left_join(coverage, bias) %>% 
  mutate(name = case_when(name == "b_Intercept" ~ "Intercept",
                           name == "b_predictor" ~ "Slope")) %>% 
  rename(Parameter = name) %>% 
  select(Parameter, `Bias (mean, sd)`, Coverage)

write_csv(table_2, file = "tables/table_2.csv")
