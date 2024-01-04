library(tidyverse)
library(rstan)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.


#1) Get fitted models
reg_update = readRDS("models/coverage_models/reg_fits_fixed.rds")

#2) Extract mean and CrI from fitted models
reg_summaries = NULL

for(i in 1:ncol(reg_update)){
  reg_summaries[[i]] = as_draws_df(reg_update[,i]$fit) %>%
    pivot_longer(starts_with("b_")) %>% 
    group_by(name) %>% 
    tidybayes::median_qi(value) %>% 
    mutate(replicate = i)
}

#3) Get true values to append
true_values = tibble(b_Intercept = -1.2,
                     b_predictor = -0.05) %>% 
  pivot_longer(cols = everything(), names_to = "name", values_to = "true_values") %>% 
  mutate(name_greek = case_when(name == "b_Intercept" ~ paste0("a) ","\u03b1"),
                                TRUE ~ paste0("b) ","\u03b2")))

#4) Append true values. Calculate coverage
reg_coverage = bind_rows(reg_summaries) %>% 
  left_join(true_values) %>% 
  mutate(contains = case_when(true_values >= .lower & true_values <= .upper ~ 1,
                              TRUE ~ 0))

#5) Make plot
reg_coverage_plot = reg_coverage %>% 
  arrange(name_greek, true_values, .lower) %>% 
  ungroup %>% 
  mutate(replicate = rep(1:1000, 2)) %>% 
  ggplot(aes(x = replicate, y = value, color = as.factor(contains),
             shape = as.factor(contains))) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper), linewidth = 0.62) +
  facet_wrap(~name_greek, scales = "free", ncol = 2) +
  geom_hline(data = . %>% distinct(name_greek, true_values), 
             aes(yintercept = true_values)) +
  scale_color_colorblind() +
  theme_default() +
  labs(y = "Parameter values",
       x = "Simulation") +
  guides(color = "none") +
  theme(strip.text = element_text(hjust = 0)) +
  coord_flip()

#6) Save plot
ggview::ggview(reg_coverage_plot, height = 3, width = 6.5)
save_plot_and_data(reg_coverage_plot, file_name = "plots/reg_coverage_plot",  height = 3, width = 6.5)
