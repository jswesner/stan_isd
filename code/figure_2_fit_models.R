library(brms)
library(tidyverse)
library(isdbayes)

# Figure 2a) sample_size -------------------------------------------------------------
# 1) Function to simulate data
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

# 2) Simulate some data
reg_dat = make_data()

# 3) fit dummy model (to update next)
reg_fit = brm(x | vreal(counts, xmin, xmax) ~ 1, 
              data = reg_dat,
              stanvars = stanvars,
              family = paretocounts(),
              # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
              chains = 4, iter = 10,
              cores = 4)

# 4) update the dummy model with different lambdas (b_values) and sample sizes (n_ind_values)
# NOTE!!!!!! This requires a cluster/supercomputer to fit. It will take > 15 hours.
b_values <- c(-2, -1.6)
n_ind_values <- c(30, 100, 300, 1000)

# 5) create empty list to store results
all_reg_fits <- list()

# 6) for loop to iterate through b_values and n_ind_values and fit models 1000 times each
for (b_val in b_values) {
  for (n_val in n_ind_values) {
    # Update the model and replicate as needed
    reg_fits <- replicate(
      update(reg_fit, newdata = make_data(n_ind = n_val, b = b_val), iter = 2000),
      n = 1000
    )
    
    # Store the results in the list
    all_reg_fits[[paste("b", b_val, "n_ind", n_val, sep = "_")]] <- reg_fits
  }
}

# 7) Save model
# saveRDS(all_reg_fits, file = "models/fig2a_mods.rds")


# Figure 2b) size_range --------------------------------------------------------------

# 1) Function to simulate data
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

# 2) simulate data
reg_dat = make_data()

# 3) fit dummy model (to update next)
reg_fit = brm(x | vreal(counts, xmin, xmax) ~ 1, 
              data = reg_dat,
              stanvars = stanvars,
              family = paretocounts(),
              # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
              chains = 4, iter = 10,
              cores = 4)



# 4) fit dummy model (to update next)
# NOTE!!!!!! This requires a cluster/supercomputer to fit. It will take > 15 hours.
b_values <- c(-2, -1.6)
x_max_values <- c(10, 100, 1000, 10000, 100000)

# 5) create empty list to store results
all_reg_fits <- list()

# 6) for loop to iterate through b_values and x_max_values and fit models 1000 times each
for (b_val in b_values) {
  for (x_max_val in x_max_values) {
    # Update the model and replicate as needed
    reg_fits <- replicate(
      update(reg_fit, newdata = make_data(xmax = x_max_val, b = b_val), iter = 2000),
      n = 1000
    )
    
    # Store the results in the list
    all_reg_fits[[paste("b", b_val, "n_ind", x_max_val, sep = "_")]] <- reg_fits
  }
}

# 4) Save the model
# saveRDS(all_reg_fits, file = "models/fig2b_mods.rds")