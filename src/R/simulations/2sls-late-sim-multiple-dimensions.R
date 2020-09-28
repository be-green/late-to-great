library(data.table)
library(magrittr)
library(rstan)
library(estimatr)
library(ggplot2)

# load helpers that construct lattice
# for treatment effect estimate
source("src/R/helpers/make_grid.R")
source("src/R/helpers/visualize_effects.R")
source("src/R/helpers/summarize_effects.R")


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

all_k <- c(1,2,3,4,5)

squared <- function(x) x^2
cubed <- function(x) x^3



all_fun <- c("I", "squared", "sin")

# simulation exercise
# n = nobs
# k = ndim for covariate grid
# ngroups = number of groups per x variable, fine-ness of grid
set.seed(100)
n = 1000

for(q in 2:length(all_k)) {

  # correct dimension for loop
  k = all_k[q]

  # construct random covariance matrix
  s = tcrossprod(rnorm(k))
  diag(s) = exp(rnorm(k))

  # simulate covariates
  X = MASS::mvrnorm(n = n, Sigma = s,
                    mu = rnorm(k, 0, 1))

  # construct random betas
  beta = rnorm(k, 0, 2)

  for(j in 1:length(all_fun)) {

    # true functional form for treatment effect
    # linear function of these covariates
    # thus treatment effects vary by individual
    # can make non-linear but for now I'll leave it
    fun = all_fun[j]
    tau = as.vector(do.call(fun, args = list(X)) %*% beta)
    tau = as.numeric(scale(tau))

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

    group_list <- c(3, 5, 10)

    for(i in 1:length(group_list)) {

      ngroups = group_list[i]

      dname <- paste0("output/model-comparisons/",
                      paste0("n-", n, ",ngroups-", ngroups,",k-",k,",fun-",fun))

      if(!dir.exists(dname)) {
        dir.create(dname)
      }

      # construct data for stan sampler
      standata <- format_data(
        X, nbins = ngroups, D = actual_treatment,
        p = pred, Y = as.numeric(outcomes)
      )

      # perform MCMC sampling
      # 5k iterations is probably overdoing it, but
      # better to have good fidelity
      fit <- sampling(model, data = standata,
                      iter = 2000,
                      control = list(max_treedepth = 15))

      # stan file already constructs pre-computed treatment effect tau
      est_tau <- extract(fit, "tau")

      # grab group index for treatment effect bucket
      group_index <- standata$group_index

      # visualize sorted effects by group vs. actual observations
      plot_effects(est_tau, tau, actual_treatment, natural_uptake, rand_assign,
                   standata$group_index)

      ggsave(filename = paste0(dname,"/posterior_vs_avg.png"), dpi = 300, width = 9, height = 6)

      plot_pred_v_raw(est_tau, tau, actual_treatment, natural_uptake, rand_assign,
                   standata$group_index)

      ggsave(filename = paste0(dname,"/posterior_interval_vs_raw.png"), dpi = 300, width = 9, height = 6)

      plot_pred_mean_v_raw(est_tau, tau, actual_treatment, natural_uptake, rand_assign,
                      standata$group_index) +
        geom_abline(slope = 1)

      ggsave(filename = paste0(dname,"/posterior_mean_vs_raw.png"), dpi = 300, width = 9, height = 6)

      # compare to point estimates from 2sls for same structure
      data.table(actual_treatment, pred = pred$fit, outcomes, group_index) %>%
        split(by = "group_index") %>%
        lapply(function(x) data.table(Estimate = coef(lm(outcomes ~ pred, data = x))[2])) %>%
        rbindlist(idcol = T) %>% setnames(".id", "group_index") %>%
        .[,group_index := as.integer(group_index)] %>%
        merge(data.table(tau, group_index, rand_assign, natural_uptake) %>%
                .[natural_uptake == 0 & rand_assign == 1] %>%
                .[,.(tau = mean(tau)), by = group_index]) %>%
        ggplot(aes(x = Estimate, y = tau, group = Estimate)) +
        geom_point() +
        theme_minimal() +
        geom_abline(slope = 1, intercept = 0)

      ggsave(filename = paste0(dname,"/2sls_vs_mean.png"), dpi = 300, width = 9, height = 6)


      fit_table <- make_fit_table(actual_treatment, pred,
                                  outcomes, group_index,
                                  est_tau, tau)


      fit_table[, k := k]
      fit_table[, fun := fun]
      fit_table[, n := n]
      fit_table[, ngroups := ngroups]

      if(file.exists("output/model-comparisons/model-summary-stats.csv")) {
        current_table <- fread("output/model-comparisons/model-summary-stats.csv")

        fit_table[,Index := max(current_table$Index) + 1]
      } else {
        fit_table[,Index := 1]
      }


      setcolorder(fit_table, c("Index", "Model", "fun", "k", "n", "ngroups",
                               "R2", "bias", "MSE", "TrueTreatmentVariance")
      )

      fwrite(fit_table, "output/model-comparisons/model-summary-stats.csv", append = T)
    }

  }
}
