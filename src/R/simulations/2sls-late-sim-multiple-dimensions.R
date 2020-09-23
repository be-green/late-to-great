library(data.table)
library(magrittr)
library(rstan)
library(estimatr)
library(ggplot2)

# load helpers that construct lattice
# for treatment effect estimate
source("src/R/helpers/make_grid.R")
source("src/R/helpers/visualize_effects.R")


# one core per chain
options(mc.cores = 4)

# compile stan model
model <- stan_model("src/stan/full-causal-model.stan")


# simulation exercise
# n = nobs
# k = ndim for covariate grid
# ngroups = number of groups per x variable, fine-ness of grid
set.seed(100)
n = 1000
k = 4
ngroups = 10

# construct random covariance matrix
z = abs(rnorm(k, 0, 2))
s = tcrossprod(z)

# simulate covariates
X = MASS::mvrnorm(n = n, Sigma = s,
                  mu = rnorm(k, 0, 1))

# construct random betas
beta = rnorm(k, 0, 10)

# true functional form for treatment effect
# linear function of these covariates
# thus treatment effects vary by individual
# can make non-linear but for now I'll leave it
tau = sin(X) %*% beta
tau = as.vector(tau)

# treatment status is a linear function of both effectiveness
# and the random component
natural_uptake = fifelse(rank(tau)/length(tau) > 0.9, 1, 0)

# strong instrument w/ lots of variance
instr = rnorm(n, mean = 0, sd = 10)
rand_assign = rbinom(n, 1, prob = plogis(instr))

# initially natural uptake, then people induced to treatment
actual_treatment = natural_uptake
actual_treatment[which(natural_uptake == 0)] <- rand_assign[which(natural_uptake == 0)]

# actual outcomes
outcomes = tau * actual_treatment + rnorm(n)

# predicted mean and standard deviation for each part of second stage
# run w/ stata robust standard errors
# use mean and sd as sufficient statistics to expand the error bars
# to adjust for stage 1 uncertainty
pred <- predict(lm_robust(actual_treatment ~ instr, se_type = "stata"),
                se.fit = T)

# construct data for stan sampler
standata <- format_data(
  X, nbins = ngroups, D = actual_treatment,
  p = pred, Y = as.numeric(outcomes)
)

# perform MCMC sampling
# 5k iterations is probably overdoing it, but
# better to have good fidelity
fit <- sampling(model, data = standata,
                iter = 5000,
                control = list(max_treedepth = 15))

# stan file already constructs pre-computed treatment effect tau
est_tau <- extract(fit, "tau")

# grab group index for treatment effect bucket
group_index <- standata$group_index

# visualize sorted effects by group vs. actual observations
plot_effects(est_tau, tau, actual_treatment, natural_uptake, rand_assign,
             standata$group_index)
