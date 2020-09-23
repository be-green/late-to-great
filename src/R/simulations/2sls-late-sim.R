library(data.table)
library(magrittr)
library(rstan)
library(data.table)

options(mc.cores = 4)

model <- stan_model("src/models/stan/1d-smoothing-prior.stan")

n = 1000

# simulation exercise
set.seed(100)

# covariance of X variables
s = matrix(c(1,0.3,0.4,
             0.2,1,0.3,
             0.4,0.2,1), ncol = 3)

# X covariates
X = rnorm(n = n, sd = 1, mean = 0)

# betas
beta = c(0.4)

# true functional form for treatment effect
# linear function of these covariates
# thus treatment effects vary by individual
tau = X * beta
tau = as.vector(tau)
# tau = as.numeric(scale(tau))
# treatment status is a linear function of both effectiveness
# and the random component
natural_uptake = fifelse(rank(tau)/length(tau) > 0.9, 1, 0)
instr = rnorm(n, mean = 0, sd = 10)
rand_assign = rbinom(n, 1, prob = plogis(instr))

# initially natural uptake, then people induced to treatment
actual_treatment = natural_uptake
actual_treatment[which(natural_uptake == 0)] <- rand_assign[which(natural_uptake == 0)]

# proportion of population induced by instrument
prop_compliers <- length(rand_assign[which(natural_uptake == 0 & rand_assign == 1)])/n

outcomes = tau * actual_treatment + rnorm(n)

scaled_outcomes = scale(outcomes)
# instr = scale(instr)
estimatr::iv_robust(outcomes ~ actual_treatment | instr)

# true effect is 0 without treatment
# we get identification on compliers, which
# is when you wouldn't have done it already
# but now you do because of the instrument
data.table(actual_treatment, natural_uptake,
           rand_assign, tau) %>%
  .[,.(Effect = mean(tau)), by = c("natural_uptake", "rand_assign")]

# bin observable characteristics for our model
bin <- function(vec, groups) {
  findInterval(vec, quantile(vec, probs = seq(0, 1, by = 1/groups)),
               all.inside = T)
}

group_list = c(10, 30, 50, 100, 500)

for(i in 1:length(group_list)) {
  ngroups = group_list[i]

  # let's create a data.frame of our groups
  X_groups <- purrr::map_dfc(as.data.frame(X), bin, groups = ngroups)
  colnames(X_groups) <- paste0(colnames(X_groups), "_group")

  # outcomes <- scale(outcomes)

  full_data <- data.table(X_groups, X, outcomes = as.numeric(outcomes),
                          actual_treatment, instrument = as.numeric(instr))


  pred_treat <- predict(lm(actual_treatment ~ instr))

  full_data <- data.table(full_data, pred_treat = as.numeric(pred_treat))

  sdata <- list(
    D = 2,
    N = n,
    L = ngroups,
    y = as.numeric(full_data$outcomes),
    x = full_data$pred_treat,
    group_index = full_data$X_group
  )

  test <- sampling(model, data = sdata, iter = 4000, control = list(max_treedepth = 20))
  betas <- extract(test, "tau")


  # test2 <- sampling(model2, data = sdata, iter = 2000, control = list(max_treedepth = 20))
  # betas2 <- extract(test2, "b_transformed")


  b_reg <- data.table(X_group = 1:ngroups, SlightlyInformed = colMeans(betas$tau))
  # b2_reg <- data.table(X_group = 1:ngroups, RandomWalk = colMeans(betas2$b_transformed))

  ind_reg <- data.table(X_groups, pred_treat = as.numeric(pred_treat),
                        outcomes = as.numeric(outcomes))[
                          ,.(CompletelyIndependent = coef(lm(outcomes ~ pred_treat))[2]),
                          by = X_group] %>%
    setorder(X_group)

  data.table(actual_treatment, natural_uptake, X_groups,
             rand_assign, tau) %>%
    .[natural_uptake == 0 & rand_assign == 1] %>%
    .[,.(TrueEffect = mean(tau), N = length(tau)),
      by = c("X_group")] %>%
    merge(b_reg) %>%
    # merge(b2_reg) %>%
    merge(ind_reg) %>%
    setorder(X_group) %>%
    .[,.(X_group, N, TrueEffect, `AR(1) Prior` = SlightlyInformed,
         `Independent 2SLS` = CompletelyIndependent
         # ,
         # `Random Walk Prior` = RandomWalk
    )] %>%
    melt(1:2, variable.name = 'Model') %>%
    ggplot(aes(x = X_group, y = value, color = Model)) +
    geom_point() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(color = "") +
    xlab("Index of X group") +
    ylab("Treatment Effect")

  ggsave(paste0("outputs/img/posterior-means-vs-ind-2sls-n-",
                n, "groups-", ngroups, "model-tausinxbeta.png"),
         dpi = 300,
         height = 5, width = 8)

  data.table(actual_treatment, natural_uptake, X_groups,
             rand_assign, tau) %>%
    .[natural_uptake == 0 & rand_assign == 1] %>%
    .[,.(TrueEffect = mean(tau), N = length(tau)),
      by = c("X_group")] %>%
    merge(b_reg) %>%
    # merge(b2_reg) %>%
    merge(ind_reg) %>%
    setorder(X_group) %>%
    .[,.(X_group, N, TrueEffect, `AR(1) Prior` = SlightlyInformed,
         `Completely Independent` = CompletelyIndependent
    )] %>%
    melt(1:3, variable.name = 'Model') %>%
    .[,.(MSE = mean((value - TrueEffect)^2),
         RMSE = mean(sqrt((value - TrueEffect)^2)),
         R2 = cor(value, TrueEffect)^2,
         Bias = mean(value - TrueEffect)), by = "Model"] %>%
  fwrite(paste0("outputs/tbl/posterior-mean-summary/posterior-means-vs-ind-2sls-n-",
                n, "groups-", ngroups, "model-tausinxbeta.png"))

}

