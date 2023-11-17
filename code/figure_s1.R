library(brms)
library(tidyverse)
library(isdbayes)
library(tidybayes)
library(scales)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.


all_fits = readRDS(file = "models/all_fits.rds")

mean_text = tibble(y = c(-1.2, -1.5, -1.8, -2)) %>% 
  mutate(x = 0.001)

prior_sens_plot = all_fits %>% 
  group_by(prior_sd, prior_mean) %>% 
  median_qi(b_Intercept) %>% 
  ggplot(aes(x = prior_sd, y = b_Intercept)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper, color = as.factor(prior_mean),
                      group = prior_mean),
                  position = position_dodge(width = 0.05)) +
  geom_line(aes(group = prior_mean, color = as.factor(prior_mean)),
            position = position_dodge(width = 0.05)) +
  scale_x_log10(labels = comma) +
  # guides(color = "none") +
  scale_y_continuous(breaks = c(-2, -1.8, -1.5, -1.2)) + 
  geom_hline(yintercept = -1.6, linetype = "dotted") +
  labs(x = "Prior standard deviation",
       y = "\u03bb",
       color = "Prior mean") +
  scale_color_colorblind() + 
  NULL

ggview::ggview(prior_sens_plot, width = 6.5, height = 4)
save_plot_and_data(prior_sens_plot, file_name = "ms/figS1", width = 6.5, height = 4, dpi = 500)


