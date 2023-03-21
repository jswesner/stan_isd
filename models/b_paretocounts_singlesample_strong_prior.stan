functions{
  real paretocounts_lpdf(real x, real lambda, real xmin, real xmax, real counts){
    if(lambda != -1)
    return(counts*(log((lambda+1) / ( xmax^(lambda+1) - xmin^(lambda+1))) + lambda*log(x)));
    else
    return(counts*(-log(log(xmax) - log(xmin)) - log(x)));
  }
}
    
data {
	int<lower=0> N;
	vector <lower = 0>[N] x;
	vector <lower = 0>[N] xmin;
	vector <lower = 0>[N] xmax;
	vector <lower = 0>[N] counts;
}

parameters {
	real lambda;
}

model {
	lambda ~ normal(-1.8, 0.2);
	for (i in 1:N){
	  x[i] ~ paretocounts(lambda, xmin[i], xmax[i], counts[i]);
	  }
}



