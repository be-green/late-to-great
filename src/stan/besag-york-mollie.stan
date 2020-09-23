functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2])
      + normal_lpdf(sum(phi) | 0, 0.00 * N);
  }
}
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure
  int<lower=1> K;                 // num covariates
  matrix[N, K] x;                 // design matrix

  real<lower=0> scaling_factor; // scales the variance of the spatial effects
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0; // intercept
  vector[K] betas; // covariates

  real <lower=0> sigma; // overall standard deviation
  real logit_rho; // proportion unstructured vs. spatially structured variance

  vector[N] theta; // heterogeneous effects
  vector[N] phi; // spatial effects
}
transformed parameters {
  real<lower=0, upper=1> rho = inv_logit(logit_rho);
  // variance of each component should be approximately equal to 1
  vector[N] convolved_re = sqrt(1 - rho) * theta
                            + sqrt(rho / scaling_factor) * phi;
}
model {
  y ~ poisson_log(log_E + beta0 + x * betas + convolved_re * sigma);  // co-variates

  // setting up our priors and model structure
  beta0 ~ normal(0.0, 1.0);
  betas ~ normal(0.0, 1.0);
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0, 1.0);
  rho ~ beta(0.5, 0.5);
  phi ~ icar_normal_lpdf(N, node1, node2);
}
generated quantities {
  real logit_rho = log(rho / (1.0 - rho));
  vector[N] eta = log_E + beta0 + x * betas + convolved_re * sigma; // co-variates
  vector[N] mu = exp(eta);
}
