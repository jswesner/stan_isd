library(tidyverse)
library(tidybayes)
library(brms)
library(ggthemes)
library(ggh4x)
source("code/sandbox/save_plot_and_data.R")

#1) read in posterior summaries 
all_lambda_summaries = readRDS(file = "posteriors/fig1abc_posterior_summaries.rds") %>% 
  mutate(true_value  = parse_number(as.character(true_value))) %>% 
  mutate(pars = "lambda")


#2) wrangle posterior summaries
plot_data = all_lambda_summaries %>% 
  # filter(xmax == 1000) %>% # limit to xmax = 1000
  mutate(cov95 = case_when(true_value > `2.5%` & true_value <`97.5%` ~ "yes", TRUE ~ "no")) %>%
  group_by(xmax, model, pars, true_value) %>% 
  # filter(replicate <= 100) %>% 
  mutate(model = case_when(model == "separate" ~ "a) Separate Models",
                           model == "fixed" ~ "b) Fixed Predictor",
                           TRUE ~ "c) Varying Intercepts"))

#3) plot posterior summaries
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
  facet_grid2(true_value ~ model, scales = "free",
              independent = "x") +
  labs(x = "\u03bb",
       y = "Model Run",
       color = "CrI contains true value?") +
  ggthemes::scale_color_colorblind() +
  # geom_text(data = labels, aes(label = true_value, x = true_value - 0.2, y = 1100),
  #           size = 3) +
  theme_default() +
  geom_vline(data = plot_data %>% ungroup %>% distinct(model, cov95, true_value),
             aes(xintercept = true_value)) +
  theme(legend.position = "top",
        strip.text.x = element_text(hjust = 0,
                                    vjust = 3),
        axis.text.x = element_text(size = 9)) +
  scale_y_continuous(breaks = c(1, 500, 1000),
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  coord_cartesian(clip = "off",
                  ylim = c(1, 1100)) +
  NULL

# 4) View and save plot
ggview::ggview(coverage_plot, width = 6.5, height = 9)

save_plot_and_data(coverage_plot, file_name = "plots/fig1coverage_plot", width = 6.5, height = 9.2)
save_plot_and_data(coverage_plot, file_name = "ms/fig1coverage_plot", width = 6.5, height = 9.2,
                   file_type = "pdf")

# 5) Summarize coverage

plot_data %>% 
  group_by(model, true_value, cov95) %>% 
  tally() %>% 
  pivot_wider(names_from = cov95, values_from = n) %>% 
  mutate(mean_coverage = 1 - no/yes) %>% 
  print(n = Inf)
