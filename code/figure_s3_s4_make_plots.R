library(tidyverse)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
source("code/sandbox/save_plot_and_data.R")

# In figure 1, the fixed and varying intercept models are 2 chains with 2000 iterations. 
# This code compares a subset of those models with 4 chains to test the effect of chains

# The models were fit on a cluster. They are the same models as described in the code for figure_1_fit_models.R
# The only difference is that we fit a subset of those models with chains = 4 and iter = 2000. For simplicity, 
# we have not added that code to the repo. Instead, we just load the wrangled posteriors here and make the plot.

# 1) load posteriors for summary statistics and rhats
posts_varint_summaries = readRDS(file = "posteriors/posts_varint_summaries.rds")
posts_varint_rhats = readRDS(file = "posteriors/posts_varint_rhats.rds")
posts_fixed_summaries = readRDS(file = "posteriors/posts_fixed_summaries.rds")
posts_fixed_rhats = readRDS(file = "posteriors/posts_fixed_rhats.rds")

# 2) Figure S3
plot_compare_chains_lambdas = bind_rows(posts_varint_summaries, 
                                                      posts_fixed_summaries) %>% 
  arrange(model, chains, true_value, `2.5%`) %>% 
  filter(chains > 1) %>% 
  ungroup %>% 
  mutate(chains = paste(model, ", Chains =", chains),
         replicate = rep(1:50, 28)) %>% 
  mutate(model = case_when(chains == "fixed , Chains = 2" ~ "a) Fixed Effects\nChains = 2",
                           chains == "fixed , Chains = 4" ~ "b) Fixed EFfects\nChains = 4",
                           chains == "varying_intercepts , Chains = 2" ~ "c) Varying Intercepts\nChains = 2",
                           TRUE ~ "d) Varying Intercepts\nChains = 4")) %>% 
  ggplot(aes(y = replicate)) + 
  geom_linerange(linewidth = 0.3,
                 aes(x = mean, xmin = `2.5%`, xmax = `97.5%`,
                     color = as.factor(chains)),
                 alpha = 0.7) +
  facet_grid(true_value ~ model, scales = "free") +
  ggthemes::scale_color_colorblind()+
  labs(color = "Number of chains",
       x = "rhat",
       y = "Model Run") +
  brms::theme_default() +
  geom_vline(aes(xintercept = true_value), linetype = 'dashed',
             linewidth = 1) +
  scale_x_continuous(breaks = c(-2.5, -1.75, -1)) +
  guides(color = "none") +
  scale_y_continuous(breaks = c(1, 25, 50))

         
ggview::ggview(plot_compare_chains_lambdas, width = 6.5, height = 6)
save_plot_and_data(plot_compare_chains_lambdas, width = 6.5, height = 6,
                   file_name = "plots/figs3_compare_chains")

save_plot_and_data(plot_compare_chains_lambdas, width = 6, height = 6,
                   file_name = "ms/figs3_compare_chains")        


# 3) Figure S4
plot_compare_chains_rhats = bind_rows(posts_varint_rhats, posts_fixed_rhats) %>% 
  filter(!parameters %in% c("lprior", "lp__")) %>% 
  filter(chains > 1) %>% 
  mutate(model = case_when(model == "fixed" ~ "a) Fixed Effects",
                           TRUE ~ "b) Varying Intercepts")) %>% 
  ggplot(aes(x = parameters, y = rhats, color = as.factor(chains))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size = 0.3) +
  ylim(0.999, 1.1) +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  coord_flip() +
  brms::theme_default() +
  labs(color = "Number of chains",
       x = "Parameters",
       y = "rhat") +
  ggthemes::scale_color_colorblind() + 
  facet_wrap(model ~., scales = "free", ncol = 1) +
  theme(legend.position = c(0.65, 0.9),
        strip.text = element_text(hjust = 0))

ggview::ggview(plot_compare_chains_rhats, width = 6, height = 9)
save_plot_and_data(plot_compare_chains_rhats, width = 6, height = 9,
                   file_name = "plots/figs4_compare_chains_rhats")

save_plot_and_data(plot_compare_chains_rhats, width = 6, height = 9,
                   file_name = "ms/figs4_compare_chains_rhats")

# Summarize ---------------------------------------------------------------

compare_lambdas = readRDS("ms/figs3_compare_chains.rds")
compare_rhats = readRDS("ms/figs4_compare_chains_rhats.rds")


compare_lambdas$data %>% 
  group_by(model) %>% 
  mutate(bias = `50%` - true_value) %>% 
  reframe(mean = mean(bias),
          sd = sd(bias))


compare_rhats$data %>% 
  group_by(model, chains) %>% 
  reframe(mean = mean(rhats),
          sd = sd(rhats),
          max = max(rhats))

