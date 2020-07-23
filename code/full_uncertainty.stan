// Ideas from 
// https://mbjoseph.github.io/posts/2018-12-25-errors-in-variables-models-in-stan/
// Stan manual on measurement error
// rstan 'getting started' school example

data {
  int<lower=0> N;       // number of cases
  vector[N] x_meas;     // measurement of x
  vector<lower=0>[N] tau;  // measurement noise
  vector[N] y;          // outcome (variate)
  vector<lower=0>[N] tau_y;
}

parameters {
  vector[N] eta;          // unknown true value
  // real mu_x;          // prior location
  // real sigma_x;       // prior scale
  real alpha;           // intercept
  real beta;            // slope
  vector<lower=0>[N] sigma_y;  // outcome noise
}

transformed parameters {
  vector[N] x;
  vector[N] sigma;
  x = x_meas + tau .* eta;
  sigma = tau_y + sigma_y;
}

model {
  target += std_normal_lpdf(eta);       // prior log-density
  target += normal_lpdf(y | alpha + beta * x, sigma); // log-likelihood
  target += normal_lpdf(alpha | 0, 10);
  target += normal_lpdf(beta | 0, 10);
  target += cauchy_lpdf(sigma_y | 0, 0.5);
// 
// model {
//   x ~ normal(mu_x, sigma_x);  // prior
//   x_meas ~ normal(x, tau);    // measurement model
//   y ~ normal(alpha + beta * x[i], sigma);
//   alpha ~ normal(0, 10);
//   beta ~ normal(0, 10);
//   sigma ~ cauchy(0, 5);
}
