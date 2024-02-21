functions {
  real negLL_PLB_binned(real b,
                        vector w,
                        int[] d,
                        int J,
                        real xmin,
                        real xmax) {
    // Check parameter bounds
    if (xmin <= 0 || xmin >= xmax || size(d) != J || size(w) != J + 1 ||
        d[1] == 0 || d[J] == 0 || min(d) < 0) {
      reject("Parameters out of bounds in negLL_PLB_binned");
    }
    
    int n = sum(d);
    real neglogLL;
    
    if (b != -1) {
      neglogLL = n * log(fabs(w[J + 1]^(b + 1) - w[1]^(b + 1))) -
                 sum(to_vector(d) .* log(fabs(w[1:J] .^ (b + 1) - w[2:(J + 1)] .^ (b + 1))));
    } else {
      neglogLL = n * log(log(w[J + 1]) - log(w[1])) -
                 sum(to_vector(d) .* log(log(w[1:J]) - log(w[2:(J + 1)])));
    }
    
    return neglogLL;
  }
}

data {
  int<lower = 0> J; // Number of bins
  vector[J + 1] w; // Bin breaks
  int<lower = 0> d[J]; // Bin counts
  real<lower = 0> xmin; // Minimum value of bins
  real<lower = 0> xmax; // Maximum value of bins
}

parameters {
  real<lower = 0> b; // Parameter b
}

model {
  // Priors for parameters can be added here if needed

  // Likelihood
  target += negLL_PLB_binned(b, w, d, J, xmin, xmax);
}
