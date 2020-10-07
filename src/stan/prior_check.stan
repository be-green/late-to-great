functions {
  // @param phi vector of spatial effects
  // @param N number of nodes, used for scaling soft constraint
  // @param node1 first vector in edge-set
  // @param node2 second vector in edge-set
  // @param s scaling factor for tau
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2])
    + normal_lpdf(sum(phi) | 0, 0.001 * N);
  }
}
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=0> N_nodes;
  int<lower=1> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1> node2[N_edges];  // and node1[i] < node2[i]
  vector[N] x_mu;
  vector[N] x_sigma;
}
parameters {
  real alpha0; // intercept
  real beta; // global slope

  real <lower=0> sigma; // residual variance
  real <lower=0> sigma_alpha; // variance for intercept RE
  real <lower=0> sigma_beta; // variance for slope RE

  vector[N_nodes] beta_i; // varying effects slopes
  vector[N_nodes] alpha_i; // spatial effects intercept
  vector[N] x;
  vector[N] y_rep;
}
transformed parameters {
  vector[N_nodes] b_transformed;
  vector[N_nodes] a_transformed;
  vector[N] mu;

  // non-centered parameterization
  b_transformed = beta_i * sigma_beta;
  a_transformed = alpha_i * sigma_alpha;

  // loop here for sampling efficiency
  for(i in 1:N) {
    mu[i] = alpha0 + a_transformed[group_index[i]] +
      x[i] * beta + x[i] * b_transformed[group_index[i]];
  }
}
model {

  // setting up our priors and model structure
  alpha0 ~ normal(0.0, 2);
  beta ~ normal(0.0, 2);

  // a_z ~ normal(0, 2);
  // b_z ~ normal(0, 2);

  sigma ~ cauchy(0, 2);
  // sigma_d ~ cauchy(0, 1);
  sigma_beta ~ normal(0, 1);
  sigma_alpha ~ normal(0, 1);

  x ~ normal(x_mu, x_sigma);
  // spatial priors on the random variation
  beta_i ~ icar_normal_lpdf(N_nodes, node1, node2);
  alpha_i ~ icar_normal_lpdf(N_nodes, node1, node2);
  // d ~ normal(x, sigma_d);
  y_rep ~ normal(mu, sigma);  // co-variates

}
generated quantities {
  vector[N_nodes] tau;

  for (i in 1:N_nodes) {
    tau[i] = beta + b_transformed[i];
  }

}
