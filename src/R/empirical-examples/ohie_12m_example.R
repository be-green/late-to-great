library(data.table)
library(rstan)
library(magrittr)
library(estimatr)

# compile stan model
options(mc.cores = 4)
model <- stan_model("src/stan/full-causal-model.stan")

# read in oregon health insurance experiment data
ohie_data <- fread("data/interim/OHIE/cleaned_ohie_data.csv")

# subset to variables of interest
ohie_data <- ohie_data[!is.na(smk_avg_mod_12m) & hhsize_12m < 100][
  is.na(employ_hrs_12m), employ_hrs_12m := 0
][,.(person_id, rx_any_12m,
                          cost_tot_oop_12m, household_id,
                          cost_any_oop_12m, ins_any_12m,
                          draw_survey_12m,
                          treatment,
                          employ_hrs_12m = make_bins(employ_hrs_12m, 6),
                          health_gen_12m, physical_act_12m,
                          birthyear = make_bins(birthyear_list, 6),
                          hhinc_cat_12m = hhinc_cat_12m,
                          hhsize_12m = as.vector(make_bins(log(hhsize_12m), 6)),
                          smk_avg_mod_12m = as.vector(make_bins(log1p(smk_avg_mod_12m), 6)),
                          edu_12m)] %>%
  .[complete.cases(.)]

setnames(ohie_data, colnames(ohie_data), stringr::str_replace_all(colnames(ohie_data), "\\.X",""))

# need our x_mu and se for first-stage results
first_stage <- predict(
  lm_robust(ohie_data$ins_any_12m ~ ohie_data$treatment, se_type = "stata"),
  se.fit = T)

X_groups <- ohie_data[,.(hhinc_cat_12m, hhsize_12m, smk_avg_mod_12m,
                         employ_hrs_12m, health_gen_12m, birthyear, physical_act_12m,
                         edu_12m)]
ncols <- ncol(X_groups)

array_dim <- sapply(X_groups, function(x) length(unique(x)))
groups <- array(1:prod(array_dim), dim = array_dim)


groups <- array(1:(nbins^ncols), dim = rep(nbins, ncols))

neighbor_edgeset <- rbindlist(furrr::future_map(1:prod(array_dim), get_id_neighbors,
                                     groups = groups))[!is.na(node2)]

X_nodes <- get_ids(X_groups, groups)
N <- nrow(X_groups)
N_edges <- nrow(neighbor_edgeset)

group_index <- get_ids(X_groups, groups)
N_nodes <- length(unique(c(neighbor_edgeset$node1, neighbor_edgeset$node2)))

list(
  N = N,
  N_edges = N_edges,
  N_nodes = N_nodes,
  node1 = neighbor_edgeset$node1,
  node2 = neighbor_edgeset$node2,
  group_index = group_index,
  y = Y,
  x_mu = pred$fit,
  x_sigma = pred$se.fit
)


