library(tidyverse)
library(brms)
source("code/pareto_brm.R")


sim_data = tibble(dw = rparetocounts(300, -1.2, 1, 1000)) %>% 
  mutate(sample = 1,
         xmin = 1, 
         xmax = 1000,
         counts = 1)

fit_brm = brm(dw | vreal(counts, xmin, xmax) ~ 1, 
              data = sim_data,
              family = paretocounts, 
              stanvars = stanvars)


fit_brm


