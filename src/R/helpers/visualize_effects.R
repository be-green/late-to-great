
plot_effects <- function(est_tau, tau, actual_treatment, natural_uptake,
                         rand_assign,
                         group_index) {
  draws <- transpose(as.data.table(est_tau$tau)) %>%
    .[,Ind := .I] %>% melt("Ind", variable.name = "Draw")
  draws <- draws[,.(UCI = quantile(value, 0.95),
                    UC50 = quantile(value, 0.75),
                    LC50 = quantile(value, 0.25),
                    LCI = quantile(value, 0.05),
                    mean = mean(value)), by = c("Ind")]

  plot_data <- data.table(actual_treatment, natural_uptake,
                          Ind = group_index,
                          rand_assign, tau = as.numeric(tau)) %>%
    .[natural_uptake == 0 & rand_assign == 1] %>%
    .[,.(TrueEffect = mean(tau), N = length(tau)), by = Ind] %>%
    merge(draws, all.y = T)

  setorder(plot_data, mean)

  plot_data[,Ind := 1:nrow(plot_data)]

  ggplot(plot_data, aes(x = Ind, y = TrueEffect, color = N)) +
    geom_pointrange(aes(y = mean, ymin = LCI, ymax = UCI),
                    color = "#009ADF") +
    geom_crossbar(aes(y = mean, ymin = LC50, ymax = UC50), width = 0.5,
                  color = "#009ADF") +
    geom_point() +
    xlab("Index") + ylab("Treatment Effect") + theme_minimal() +
    scale_color_distiller(palette = "Oranges")
}

plot_pred_v_raw <- function(est_tau, tau, actual_treatment, natural_uptake,
                         rand_assign,
                         group_index) {
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

  ggplot(plot_data, aes(x = Ind, y = TrueEffect, group = Ind)) +
    geom_pointrange(aes(y = mean, ymin = LCI, ymax = UCI),
                    color = "#009ADF") +
    geom_crossbar(aes(y = mean, ymin = LC50, ymax = UC50), width = 0.5,
                  color = "#009ADF") +
    stat_summary(fun.data = "mean_se") +
    xlab("Index") + ylab("Treatment Effect") + theme_minimal() +
    ggtitle("Posterior (blue) vs. True Effect (black)",
            subtitle = "Line is 95% interval, box is 50% interval")
}


plot_pred_mean_v_raw <- function(est_tau, tau, actual_treatment, natural_uptake,
                            rand_assign,
                            group_index) {
  draws <- transpose(as.data.table(est_tau$tau)) %>%
    .[,Ind := .I] %>% melt("Ind", variable.name = "Draw")
  draws <- draws[,.(UCI = quantile(value, 0.95),
                    UC50 = quantile(value, 0.75),
                    LC50 = quantile(value, 0.25),
                    LCI = quantile(value, 0.05),
                    mean = mean(value)), by = c("Ind")]

  plot_data <- data.table(actual_treatment, natural_uptake,
                          Ind = group_index,
                          rand_assign, tau = as.numeric(tau)) %>%
    .[natural_uptake == 0 & rand_assign == 1] %>%
    .[,.(TrueEffect = mean(tau), N = length(tau)), by = Ind] %>%
    merge(draws, all.y = T)

  setorder(plot_data, mean)

  plot_data[,Ind := 1:nrow(plot_data)]

  ggplot(plot_data, aes(x = mean, y = TrueEffect)) +
    geom_point() +
    xlab("Posterior Mean") + ylab("Treatment Effect") + theme_minimal() +
    scale_color_distiller(palette = "Oranges")
}
