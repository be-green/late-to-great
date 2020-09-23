# ivregress 2sls `var' (ohp_all_ever_`source'= treatment) `controls' ``var'_controls' ///
#   if 1==1 `condition' [pw = `weight'], cluster(household_id)
#

library(estimatr)
library(data.table)
library(rstan)
library(magrittr)

options(mc.cores = 2)

ohie_data <- fread("data/interim/OHIE/cleaned_ohie_data.csv")

ohie_data <- ohie_data[,.(person_id, rx_any_12m,
                          cost_tot_oop_12m, household_id,
                          weight_12m,
                          cost_any_oop_12m, ins_any_12m,
                          draw_survey_12m,
                          treatment, hhinc_cat_12m,
                          edu_12m, num19_12m, employ_12m,
                          race_white_12m, race_hisp_12m,
                          race_black_12m,
                          race_amerindian_12m,
                          race_asian_12m,
                          race_pacific_12m,
                          race_other_qn_12m)] %>%
  .[complete.cases(.)]

regfit <- iv_robust(log1p(cost_tot_oop_12m) ~ 0 + ins_any_12m +
            factor(hhinc_cat_12m) +
            factor(edu_12m) + employ_12m +
            race_white_12m +
            race_hisp_12m + race_black_12m +
            race_amerindian_12m + race_asian_12m +
            race_pacific_12m + race_other_qn_12m |
            treatment + factor(hhinc_cat_12m) +
            factor(edu_12m) + employ_12m + race_white_12m +
            race_hisp_12m + race_black_12m +
            race_amerindian_12m + race_asian_12m +
            race_pacific_12m + race_other_qn_12m,
        se_type = "stata",
        weights = weight_12m,
        clusters = household_id,
        data = ohie_data)

stanversion <- stan_model(file = "src/models/stan/2sls.stan")

# matrix of exogenous x regressors
# that are not instruments
x_exog <- model.matrix(
  log1p(cost_tot_oop_12m) ~ 0  +
    factor(hhinc_cat_12m) +
    factor(edu_12m) + employ_12m +
    race_white_12m +
    race_hisp_12m + race_black_12m +
    race_amerindian_12m + race_asian_12m +
    race_pacific_12m + race_other_qn_12m,
  data = ohie_data
)

first_stage <- brms::bf(
  log1p(cost_tot_oop_12m) ~ 0  +
      factor(hhinc_cat_12m) +
      factor(edu_12m) + employ_12m +
      race_white_12m +
      race_hisp_12m + race_black_12m +
      race_amerindian_12m + race_asian_12m +
      race_pacific_12m + race_other_qn_12m)

second_stage <- brms::bf(
  ins_any_12m ~ treatment  +
    factor(hhinc_cat_12m) +
    factor(edu_12m) + employ_12m +
    race_white_12m +
    race_hisp_12m + race_black_12m +
    race_amerindian_12m + race_asian_12m +
    race_pacific_12m + race_other_qn_12m)

iv <- brm(first_stage + second_stage, data = ohie_data)


pred_data <- posterior_linpred(first_stage, nsamples = 100) %>%
  as.data.table %>%
  .[,I := .I] %>%
  split(by = "I") %>%
  lapply(function(x) {
    x[,I := NULL]
    data.table(ohie_data, x_pred = transpose(x))
  }) %>%
  lapply(function(x) setnames(x, "x_pred.V1", "x_pred"))


second_stage <- brms::brm_multiple(
  log1p(cost_tot_oop_12m) ~ 0 + x_pred +
    factor(hhinc_cat_12m) +
    factor(edu_12m) + employ_12m +
    race_white_12m +
    race_hisp_12m + race_black_12m +
    race_amerindian_12m + race_asian_12m +
    race_pacific_12m + race_other_qn_12m,
  data = pred_data)


# right now this only works with a single variable
# to be instrumented for
stanfit2 <-
  sampling(stanversion,
             data = list(
               N = nrow(ohie_data),
               K_x = ncol(x_exog),
               K_z = 1,
               x_end = ohie_data$ins_any_12m,
               x_exog = scale(x_exog),
               z = as.matrix(ohie_data$treatment),
               y = as.vector(scale(log1p(ohie_data$cost_tot_oop_12m))
             ), iter = 2000
           )
  )
