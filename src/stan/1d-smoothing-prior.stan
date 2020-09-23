data {
  int<lower=0> N;
  int<lower=1> L;
  real y[N];
  int<lower=1,upper=L> group_index[N];
  vector[N] x;
}
parameters {
  real intercept; // intercept
  real beta; // global slope
  real<lower=0> sigma_err;
  real<lower=0> sigma_b;
  real<lower=0> sigma_a;
  real<lower=0> s_b;
  real<lower=0> s_a;
  vector[L] alpha_l; // slope
  vector[L] beta_l; // slope
  real<lower=0,upper=1> rho_a;
  real<lower=0,upper=1> rho_b;
}
transformed parameters {
  vector[N] yhat;
  real<lower=-1,upper=1> rho_transformed_a;
  real<lower=-1,upper=1> rho_transformed_b;
  vector[L] b_transformed;
  vector[L] a_transformed;

  rho_transformed_a = (rho_a * 2) - 1; // the autoregressive coefficient
  rho_transformed_b = (rho_b * 2) - 1; // the autoregressive coefficient

  for(l in 1:L) {
    b_transformed[l] = sigma_b * beta_l[l];
    a_transformed[l] = sigma_a * alpha_l[l];
  }
  for (i in 1:N) {
    yhat[i] = intercept + x[i] * beta + a_transformed[group_index[i]] + x[i] * b_transformed[group_index[i]];
  }
}
model {
  // priors on parameter variance
  sigma_b ~ normal(0, 1);
  sigma_a ~ normal(0, 1);

  // priors on parameters
  intercept ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma_err ~ cauchy(0, 2);
  rho_a ~ beta(2, 2);
  rho_b ~ beta(2, 2);

  s_a ~ cauchy(0,1);
  s_b ~ cauchy(0,1);

  alpha_l[1] ~ normal(0, 1/sqrt(1-rho_transformed_a^2)); // before it was normal(0, 1) but this is wrong
  beta_l[1] ~ normal(0, 1/sqrt(1-rho_transformed_b^2)); // before it was normal(0, 1) but this is wrong

  for (j in 2:L) {
    alpha_l[j] ~ normal(rho_transformed_a * alpha_l[j-1],1);
    beta_l[j] ~ normal(rho_transformed_b * beta_l[j-1],1);
  }

  for (n in 1:N)
    y[n] ~ normal(yhat[n], sigma_err);
}
generated quantities {
  vector[L] tau;
  for (l in 1:L) {
    tau[l] = beta + b_transformed[l];
  }
}
