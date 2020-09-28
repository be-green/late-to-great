
make_fit_table <- function(actual_treatment, pred, outcomes, group_index,
                           est_tau, tau) {

  tsls <- data.table(actual_treatment, pred = pred$fit, outcomes, group_index) %>%
    split(by = "group_index") %>%
    lapply(function(x) data.table(Estimate = coef(lm(outcomes ~ pred, data = x))[2])) %>%
    rbindlist(idcol = T) %>% setnames(".id", "group_index") %>%
    .[,group_index := as.integer(group_index)] %>%
    merge(data.table(tau, group_index, rand_assign, natural_uptake) %>%
            .[natural_uptake == 0 & rand_assign == 1] %>%
            .[,.(tau = tau), by = group_index]) %>%
    .[!is.na(Estimate)] %>%
    .[,.(R2 = cor(tau, Estimate)^2, bias = mean(tau - Estimate),
         MSE = mean((tau - Estimate)^2),
         TrueTreatmentVariance = sd(tau)^2)]
  tsls[,Model := "2SLS"]


  draws <- transpose(as.data.table(est_tau$tau)) %>%
    .[,Ind := .I] %>% melt("Ind", variable.name = "Draw")
  draws <- draws[,.(UCI = quantile(value, 0.95),
                    UC50 = quantile(value, 0.75),
                    LC50 = quantile(value, 0.25),
                    LCI = quantile(value, 0.05),
                    mean = mean(value)), by = c("Ind")]

  plot_data <- data.table(actual_treatment, natural_uptake,
                          Ind = group_index,
                          rand_assign, TrueEffect = as.numeric(tau)) %>%
    .[natural_uptake == 0 & rand_assign == 1] %>%
    merge(draws, all.y = T)

  setorder(plot_data, mean)

  plot_data[,Ind := factor(Ind, levels = setorder(unique(plot_data[,.(Ind, mean)]), mean)$Ind)]

  plot_data[,Ind := as.integer(Ind)]

  gmrf_fit <- plot_data %>%
    .[!is.na(TrueEffect)] %>%
    .[,.(R2 = cor(TrueEffect, mean)^2, bias = mean(TrueEffect - mean),
         MSE = mean((TrueEffect - mean)^2),
         TrueTreatmentVariance = sd(TrueEffect)^2)]

  gmrf_fit[,Model := "GMRF Prior"]

  rbind(tsls, gmrf_fit)
}

