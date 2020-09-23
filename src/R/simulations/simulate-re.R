# Simulating varying treatment effects

# number of groups
J = 30
N = 100

# average effect of treatment
avg_effect <- 0.1

# average of response variable overall
y_mean <- 0

q = 0.2
r = 0.9
s = 1.5

# covariance matrix for slope and intercept terms
covmat <- matrix(c(q^2, r * q * s, r * q * s, s^2), nrow = 2,
                     byrow = TRUE)
re <- rmvnorm(J, mean = c(0, 0), sigma = covmat)

# group level averages randomly vary
alphas <- y_mean + re[, 1]

# group_level response randomly varying as well
betas <- avg_effect + re[, 2]

# treatment uptake by group (e.g. does the lottery induce uptake of treatment)
# for now the same across all groups
uptake_prob <- rep(0.5, J)

sim_group_data <- function(J, N, alphas, betas, uptake_prob) {

  group_data <- data.table(J = 1:30, alphas, betas, uptake_prob)
  group_data <- group_data[rep(1:.N,N)]

  setorder(group_data, J)

  group_data[,N := rep(1:100, 30)]

  # half of subjects get the lottery
  group_data[,Lottery := runif(nrow(group_data))]
  group_data[,Treatment := fifelse(Lottery > 0.5, 1, 0)]
  group_data[,UptakeSampling := runif(nrow(group_data))]
  group_data[Treatment == 1 & UptakeSampling < uptake_prob, Uptake := 1]
  group_data[is.na(Uptake) & UptakeSampling < uptake_prob/2, Uptake := 1]
  group_data[is.na(Uptake), Uptake := 0]

  err_scale <- 0.2

  group_data[,Y := alphas + betas * Uptake + rnorm(nrow(group_data), 0, err_scale)]
  group_data[]
}

library(estimatr)

group_data <- sim_group_data(J, N, alphas, betas, uptake_prob)

# Adjust for group level differences with group regressor
# 2sls basically nails the result
iv_robust(Y ~ Uptake + factor(J) | Treatment + factor(J), data = group_data)

# overestimated by linear regression
lm(Y ~ Uptake, data= group_data)

# actual treatment effect
group_data[,.(mean(betas)), by = Uptake]

# what if we have sorting into the uptake?
# then we get a LATE, but it's hard to extrapolate

alphas = sort(alphas)
betas = sort(betas)
uptake_prob = seq(1/J, 1, by = 1/J)

group_data <- sim_group_data(J, N, alphas, betas, uptake_prob)

# Adjust for group level differences with group regressor
# 2sls basically nails the result
# where Uptake is a weighted average of all people selected into the program
iv_robust(Y ~ Uptake + factor(J) | Treatment + factor(J), data = group_data)
group_data[,.(mean(betas)), by = Uptake]

# But what if we want to extrapolate? Then we need a model for the variation
# in the actual treatment effect
iv_robust(Y ~ Uptake + factor(J) | Treatment + factor(J), data = group_data)

first_stage <- glmer(Uptake ~ Treatment + (1 + Treatment | J),
                          data = group_data,
                     control = glmerControl(optimizer = "Nelder_Mead"),
                   family = binomial)

pred = predict(first_stage, type = "response")
group_data[,Pred := pred]

fit <- lmer(Y ~ Treatment + (1 + Treatment | J), data = group_data )

