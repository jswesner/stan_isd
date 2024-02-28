library(tidyverse)
library(brms)
source("code/sandbox/save_plot_and_data.R") # custom function that saves image files with ggsave and also makes an .rds copy of that file.

#1) load posteriors

fig3_posterior_summaries = readRDS(file = "posteriors/fig3_posterior_summaries.rds")

#2) make plot
reg_coverage_plot = fig3_posterior_summaries %>% 
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
save_plot_and_data(reg_coverage_plot, file_name = "ms/Figure3", file_type = "pdf",  height = 3, width = 6.5)
