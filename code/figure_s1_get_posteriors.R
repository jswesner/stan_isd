library(tidyverse)
library(rstan)
library(isdbayes)
library(brms)

# 1) load models
# these are too big for GitHub (500-1000MB), so we have to load from Zenodo.
# fish_invert_mods_2 = unlist(readRDS(file = "models/figs1a_mods.rds"), recursive = F)
# fish_invert_mods_1.6 = unlist(readRDS(file = "models/figs1b_mods.rds"), recursive = F)
# fish_invert_mods_1.2 = unlist(readRDS(file = "models/test/figs1c_mods.rds"), recursive = F)

fish_invert_mods_2 = unlist(readRDS(url("https://zenodo.org/records/10698127/files/figs1a_mods.rds?download=1")), recursive = F)
fish_invert_mods_1.6 = unlist(readRDS(url("https://zenodo.org/records/10698127/files/figs1b_mods.rds?download=1")), recursive = F)
fish_invert_mods_1.2 = unlist(readRDS(url("https://zenodo.org/records/10698127/files/figs1c_mods.rds?download=1")), recursive = F)

# 2) Extract posterior summaries
fish_invert_posts_2_list = NULL

for(i in 1:length(fish_invert_mods_2)){
  fish_invert_posts_2_list[[i]] = as_draws_df(fish_invert_mods_2[[i]]) %>% 
    mutate(true_lambda = -2, 
           replicate = i) %>% 
    mutate(data_structure = case_when(i <= 100 ~ "1) Fish\nCounts = 1",
                                      i <= 200 ~ "2) Inverts\nCounts = 1",
                                      i <= 300 ~ "3) Fish + Inverts\nCounts = 1",
                                      TRUE ~ "4) Fish + Inverts\nCounts = 1/sample_area"),
           replicate = case_when(i <= 100 ~ replicate,
                                 i <= 200 ~ replicate - 100,
                                 i <= 300 ~ replicate - 200,
                                 TRUE ~ replicate - 300))
}

fish_invert_posts_1.6_list = NULL

for(i in 1:length(fish_invert_mods_1.6)){
  fish_invert_posts_1.6_list[[i]] = as_draws_df(fish_invert_mods_1.6[[i]]) %>% 
    mutate(true_lambda = -1.6, 
           replicate = i) %>% 
    mutate(data_structure = case_when(i <= 100 ~ "1) Fish\nCounts = 1",
                                      i <= 200 ~ "2) Inverts\nCounts = 1",
                                      i <= 300 ~ "3) Fish + Inverts\nCounts = 1",
                                      TRUE ~ "4) Fish + Inverts\nCounts = 1/sample_area"),
           replicate = case_when(i <= 100 ~ replicate,
                                 i <= 200 ~ replicate - 100,
                                 i <= 300 ~ replicate - 200,
                                 TRUE ~ replicate - 300))
}

fish_invert_posts_1.2_list = NULL

for(i in 1:length(fish_invert_mods_1.2)){
  fish_invert_posts_1.2_list[[i]] = as_draws_df(fish_invert_mods_1.2[[i]]) %>% 
    mutate(true_lambda = -1.2, 
           replicate = i) %>% 
    mutate(data_structure = case_when(i <= 100 ~ "1) Fish\nCounts = 1",
                                      i <= 200 ~ "2) Inverts\nCounts = 1",
                                      i <= 300 ~ "3) Fish + Inverts\nCounts = 1",
                                      TRUE ~ "4) Fish + Inverts\nCounts = 1/sample_area"),
           replicate = case_when(i <= 100 ~ replicate,
                                 i <= 200 ~ replicate - 100,
                                 i <= 300 ~ replicate - 200,
                                 TRUE ~ replicate - 300))
}

# 3) bind and save
all_posts = bind_rows(fish_invert_posts_1.2_list,
                      fish_invert_posts_1.6_list,
                      fish_invert_posts_2_list)

saveRDS(all_posts, file = "posteriors/figs1_posterior_summaries.rds")





# try it with sizeSpectra code

mlecounts = NULL

dat_2 = replicate(sim_dat(mu = -2)[4], n = 2)

dat_2

for(i in 1:length(dat_2)){
x = dat_2[[i]]$x
c = dat_2[[i]]$counts
sumclogx = sum(c*log(x))


mlecounts[[i]] <- sizeSpectra::calcLike(negLL.fn = negLL.PLB.counts,
                                       p = -2,
                                       x = x,
                                       c = c,
                                       K = length(c),
                                       xmin = min(x),
                                       xmax = max(x),
                                       sumclogx = sumclogx)

}


mles = NULL
for(i in 301:400){
mles[[i]] = tibble(b_mle = mlecounts[[i]]$MLE)
}

cis = NULL
for(i in 301:400){
  cis[[i]] = tibble(low_mle = mlecounts[[i]]$conf[1],
                     high_mle = mlecounts[[i]]$conf[2])
}

mean_ci = bind_rows(mles) %>% 
  mutate(lower = bind_rows(cis) %>% select(low_mle),
         upper = bind_rows(cis) %>% select(high_mle))


hist(mean_ci$b_mle)
