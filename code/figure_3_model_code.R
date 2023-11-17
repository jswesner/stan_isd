library(brms)
library(tidyverse)
library(isdbayes)
library(ggthemes)
source("code/sandbox/save_plot_and_data.R")

reps = 4
xvals = 3
n_ind = 300
intercept = -1.2
var_int_sd = 0.01
beta = -0.05
n_groups = 1
xmin = 1
xmax = 1000

set.seed(23434)

data_grid = tibble(group = (1:n_groups)) %>% 
  expand_grid(predictor = -seq(-1, 1, length.out = 5)) %>% 
  expand_grid(replicates = 1:reps) %>% 
  mutate(offset = rnorm(nrow(.), 0, var_int_sd ),
         group_b = intercept) %>% 
  mutate(b = intercept + offset + predictor*beta,
         id = 1:nrow(.)) %>% 
  mutate(par = letters[id],
         r_true = b - group_b)

data_grid %>% 
  ggplot(aes(x = predictor, y = b)) +
  geom_point()

make_regression_data = function(n_ind = 300,
                                xmin = 1, 
                                xmax = 1000){
  
  data_grid %>% 
    expand_grid(individual = 1:n_ind) %>% 
    mutate(xmin = xmin, 
           xmax = xmax,
           u = runif(nrow(.), min = 0, max = 1),
           x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% 
    group_by(group, id, xmin, xmax) %>% 
    add_count(x, name = "counts") %>% 
    ungroup
}

reg_dat = make_regression_data()

# fixed effects model
reg_fit = brm(x | vreal(counts, xmin, xmax) ~ predictor + (1|id), 
              data = reg_dat,
              stanvars = stanvars,
              family = paretocounts(),
              # prior = c(prior(normal(-1.5, 0.01), class = "Intercept")),
              chains = 1, iter = 10,
              cores = 4)


reg_fits = replicate(update(reg_fit, newdata = make_regression_data(), iter = 1000), n = 2)

# saveRDS(reg_fits, file = "models/coverage_models/reg_fits.rds")