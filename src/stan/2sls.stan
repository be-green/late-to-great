data {
  int<lower=0> N; // rows
  int<lower=0> K_x; // num exog x vars
  int<lower=0> K_z; // num instruments
  matrix[N, K_z] z; // matrix of instrumental variables
  vector[N] x_end; // single endogenous x
  matrix[N, K_x] x_exog; // exogenous outcomes
  vector[N] y; // outcome variable
}
parameters {
  real alpha_first;
  vector[K_z + K_x] beta_first;
  real alpha_second;
  vector[K_x + 1] beta_second;
  real<lower=0> sigma_x;
  real<lower=0> sigma_y;
}
model {
  vector[N] x_pred;
  matrix[N, K_x + K_z] all_exog;
  matrix[N, K_x + 1] scnd_stage;

  // first stage, project on exogenous variables
  all_exog = append_col(z, x_exog);
  x_end ~ normal(alpha_first + all_exog * beta_first, sigma_x);

  // second stage, use exogenous component of end. treatment
  x_pred = alpha_first + all_exog * beta_first;
  scnd_stage = append_col(x_exog, x_pred);

  // final outcome
  y ~ normal(alpha_second + scnd_stage * beta_second, sigma_y);
}
