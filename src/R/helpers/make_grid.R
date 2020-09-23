#' Bin observable characteristics for our model
bin <- function(vec, groups) {
  findInterval(vec, quantile(vec, probs = seq(0, 1, by = 1/groups)),
               all.inside = T)
}

#' Make covariate data into a set of bins
#' @param X covariates to use in grid creation
#' @param bins number of bins to be used in grid
#' @importFrom purrr map_dfc
#' @importFrom data.table as.data.table
make_bins <- function(X, nbins) {
    data.table::as.data.table(
      purrr::map_dfc(as.data.frame(X), bin, groups = nbins)
    )
}

#' Create array representing group locations
#' @param nbins number of bins used for grid for each variable
#' @param ncols number of variables used to create grid
make_group_array <- function(nbins, ncols) {
  array(1:(nbins^ncols), dim = rep(nbins, ncols))
}

vec_to_index <- function(index_vector) {
  paste0(index_vector, collapse = ",")
}

get_id <- function(index_vector, groups) {
  eval(parse(text = paste0("groups[", vec_to_index(index_vector),"]")))
}

get_index <- function(id, groups) {
  as.vector(which(groups == id, arr.ind = T))
}

#' Given a set of group indexes, find place in vector ordering
get_ids <- function(X_groups, groups) {
  apply(X_groups,
        MARGIN = 1,
        get_id,
        groups = groups)
}

get_index_neighbors <- function(index_vector, groups) {

  mat <- matrix(ncol = length(index_vector))

  nbr_tbl <- rbindlist(lapply(1:length(index_vector), function(i) {
    if(index_vector[i] < nrow(groups)) {
      moved <- index_vector
      moved[i] <- index_vector[i] + 1
      data.table(
        t(moved)
      )
    }
  }))

  data.table(node1 = get_id(index_vector, groups), node2 = get_ids(nbr_tbl, groups))
}

get_id_neighbors <- function(node_id, groups) {
  get_index_neighbors(get_index(node_id, groups), groups)
}


#' Format data for stancode
#' @param X covariates for grid
#' @param nbins number of bins to be used in grid construction
#' @param D observed treatment uptake
#' @param Z instrument
#' @param Y outcomes
format_data <- function(X, nbins, D, p, Y){
  X_groups <- make_bins(X, nbins)
  ncols <- ncol(X_groups)

  groups <- make_group_array(nbins, ncols)

  neighbor_edgeset <- rbindlist(lapply(1:(nbins^ncols), get_id_neighbors,
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
}

