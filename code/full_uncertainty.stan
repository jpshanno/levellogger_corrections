// Ideas from 
// https://mbjoseph.github.io/posts/2018-12-25-errors-in-variables-models-in-stan/
// Stan manual on measurement error
// rstan 'getting started' school example

data {
  int<lower=0> N;       // number of cases
  vector[N] x_obs;     // measurement of x
  vector<lower=0>[N] tau_x;  // measurement noise
  vector[N] y_obs;          // outcome (variate)
  vector<lower=0>[N] tau_y;
}

parameters {
  vector[N] z_x;
  vector[N] z_y; 
  real alpha;           // intercept
  real beta;            // slope
  vector<lower=0>[N] sigma;  // outcome noise
  
  // real mu_x;          // prior location
  // real sigma_x;       // prior scale
}

transformed parameters {
  vector[N] x;
  vector[N] y;
  
  x = x_obs + tau_x .* z_x;
  y = y_obs + tau_y .* z_y;
}

model {
  target += std_normal_lpdf(z_x);
  target += std_normal_lpdf(z_y);
  target += normal_lpdf(y | alpha + beta * x, sigma); // log-likelihood
  target += normal_lpdf(alpha | 0, 10);
  target += normal_lpdf(beta | 0, 10);
  target += cauchy_lpdf(sigma | 0, 0.5);
}
