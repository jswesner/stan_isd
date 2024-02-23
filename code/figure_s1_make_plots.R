library(tidyverse)
library(rstan)
library(isdbayes)
library(brms)
library(tidybayes)

# 1) load posteriors
all_posts = readRDS(file = "posteriors/figs1_posterior_summaries.rds")

# 2) plot and save
fig_counts_fishinverts = all_posts %>% 
  # filter(replicate <= 20) %>%
  mutate(facet = case_when(true_lambda == -2 ~ "a)",
                           true_lambda == -1.6 ~ "b)",
                           TRUE ~ "c)")) %>% 
  ggplot(aes(x = b_Intercept, y = reorder(data_structure, desc(data_structure)))) +
  stat_slab(aes(group = interaction(replicate, data_structure),
                color = as.factor(true_lambda)),
            # color = "black", 
            fill = NA, linewidth = 0.01, alpha = 0.8) + 
  facet_wrap(~facet) + 
  geom_vline(aes(xintercept = true_lambda), linetype = "dashed") +
  labs(y = "",
       x = "\u03bb") + 
  theme_bw() +
  ggthemes::scale_color_colorblind() +
  guides(color = "none") +
  theme(strip.text = element_text(hjust = 0),
        panel.grid = element_blank(),
        strip.background = element_blank())

ggview::ggview(fig_counts_fishinverts, width = 6.5, height = 3)
ggsave(fig_counts_fishinverts, file = "plots/figs1_fish_counts_fishinverts.jpg",
       width = 6.5, height = 3)
saveRDS(fig_counts_fishinverts, file = "plots/figs1_fish_counts_fishinverts.rds")
ggsave(fig_counts_fishinverts, file = "ms/figs1_fish_counts_fishinverts.jpg",
       width = 6.5, height = 3)


