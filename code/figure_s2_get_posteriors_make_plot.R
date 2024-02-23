library(tidyverse)
library(rstan)
library(isdbayes)
library(brms)

mods_area = readRDS(file = "models/figs2_mods.rds")

# get posts and plot

mods_posts = NULL
for(i in 1:length(mods_area)){
  mods_posts[[i]] = as_draws_df(mods_area[[i]]) %>% 
    mutate(model = i)
}

plot_sampling_area = bind_rows(mods_posts) %>% 
  mutate(area = case_when(model == 3 ~ "c) 10m2",
                          model == 2 ~ "b) m2",
                          model == 1 ~ "a) 0.1m2")) %>% 
  ggplot(aes(x = b_Intercept, fill = area)) + 
  viridis::scale_fill_viridis(discrete = T) +
  stat_slab(alpha = 0.8) +
  # geom_vline(xintercept = -1.6, linetype = "dashed") +
  theme_default() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.8)) + 
  labs(fill = "Sampling Area",
       x = "\u03bb")


saveRDS(plot_sampling_area, file = "plots/figure_s2_sampling_area.rds")
ggsave(plot_sampling_area, file = "plots/figure_s2_sampling_area.jpg",
       width = 6, height = 4, dpi = 500)
ggsave(plot_sampling_area, file = "ms/figure_s2_sampling_area.jpg",
       width = 6, height = 4, dpi = 500)
