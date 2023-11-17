library(brms)
library(tidyverse)
library(isdbayes)


# Sample Size -------------------------------------------------------------

#1) Function to simulate data
make_data = function(n_ind = 300,
                     xmin = 1, 
                     xmax = 1000,
                     b = -1.6){
  
  tibble(xmin = xmin,
         xmax = xmax,
         b = b) %>% 
    expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = xmin, 
           xmax = xmax,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    mutate(counts = 1) %>% 
    ungroup
}

#2) Precompile dummy model (for updating later)
reg_dat = make_data()

# precompile dummy model
reg_fit = brm(x | vreal(counts, xmin, xmax) ~ 1, 
              data = reg_dat,
              stanvars = stanvars,
              family = paretocounts(),
              # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
              chains = 1, iter = 10,
              cores = 4)



#3) Update dummy model with different lambdas (b_values) and sample sizes (n_ind_values)
# NOTE!!!!!! This requires a cluster/supercomputer to fit. To run it, unsilence the for loop.
b_values <- c(-2, -1.6)
n_ind_values <- c(30, 100, 300, 1000)

all_reg_fits <- list()

# for (b_val in b_values) {
#   for (n_val in n_ind_values) {
#     # Update the model and replicate as needed
#     reg_fits <- replicate(
#       update(reg_fit, newdata = make_data(n_ind = n_val, b = b_val), iter = 1000),
#       n = 1000
#     )
#     
#     # Store the results in the list
#     all_reg_fits[[paste("b", b_val, "n_ind", n_val, sep = "_")]] <- reg_fits
#   }
# }

#4) Save the model
# saveRDS(all_reg_fits, file = "models/sample_size1000.rds")

# load fitted model
sample_size_brm = readRDS(file = "models/sample_size1000.rds")

#5) Summarize posteriors and save

temp_stats <- vector("list", length = 1000)

for (i in 1:1000) {
  # Create a nested list for each i iteration
  temp_stats[[i]] <- vector("list", length = 2)
  
  for (b in 1:8) {
    temp_stats[[i]][[b]] <- sample_size_brm[[b]][14, i][[1]] %>%
      as_draws_df() %>%
      tidybayes::median_qi(b_Intercept) %>% 
      mutate(replicate = i,
             n = nrow(sample_size_brm[[b]][2, i][[1]]))
  }
}


sample_size_posts_df = bind_rows(temp_stats)

saveRDS(sample_size_posts_df, file = "posteriors/sample_size_posts_df.rds")

# Size Ranges -------------------------------------------------------------
library(brms)
library(tidyverse)
library(isdbayes)

#1) Function to simulate data
make_data = function(n_ind = 300,
                     xmin = 1, 
                     xmax = 1000,
                     b = -1.6){
  
  tibble(xmin = xmin,
         xmax = xmax,
         b = b) %>% 
    expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = xmin, 
           xmax = xmax,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    mutate(counts = 1) %>% 
    ungroup
}

#2) Precompile dummy model (for updating later)
reg_dat = make_data()

# precompile dummy model
reg_fit = brm(x | vreal(counts, xmin, xmax) ~ 1, 
              data = reg_dat,
              stanvars = stanvars,
              family = paretocounts(),
              # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
              chains = 1, iter = 10,
              cores = 4)

#3) Update dummy model with different lambdas (b_values) and sample sizes (n_ind_values)
# NOTE!!!!!! This requires a cluster/supercomputer to fit. To run it, unsilence the for loop.
b_values <- c(-2, -1.6)
x_max_values <- c(10, 100, 1000, 10000, 100000)

all_reg_fits <- list()

for (b_val in b_values) {
  for (x_max_val in x_max_values) {
    # Update the model and replicate as needed
    reg_fits <- replicate(
      update(reg_fit, newdata = make_data(xmax = x_max_val, b = b_val), iter = 1000),
      n = 1
    )

    # Store the results in the list
    all_reg_fits[[paste("b", b_val, "n_ind", x_max_val, sep = "_")]] <- reg_fits
  }
}

# 4) Save the model
saveRDS(all_reg_fits, file = "models/size_range_brm1000.rds")

# load fitted model
sample_size_brm = readRDS(file = "models/size_range_brm1000.rds")


#5) Summarize posteriors and save

dat = tibble(xmin = xmin,
             xmax = xmax) %>% 
  expand_grid(b = c(-2, -1.6)) %>% 
  expand_grid(n_sim = c(30, 100, 300, 1000)) %>% 
  expand_grid(replicate = 1:1000) %>% 
  group_by(b, n_sim, replicate) %>% 
  group_split()

sample_size_posts = NULL

for(i in 1:length(sample_size_brm)){
  sample_size_posts[[i]] = as_draws_df(sample_size_brm[[i]]) %>% 
    mutate(sample_size = nrow(sample_size_brm[[i]]$data),
           true_value = sample_size_brm[[i]]$data2$b,
           n_sim = sample_size_brm[[i]]$data2$n_sim,
           replicate = sample_size_brm[[i]]$data2$replicate)
}

sample_size_posts_df = bind_rows(sample_size_posts)
saveRDS(sample_size_posts_df, file = "posteriors/sample_size_posts_df.rds")