# LATE to Great

This repository is the code behind an early-stage project for direct estimation and extrapolation of varying treatment effects. The basic approach is one taken from survey statistics, where subjects are divided into groups and have responses estimated by survey "cell" in order to estimate high-order interaction effects. To predict the responses of a new population, the cells are simply re-weighted (in political science using voter files or census data). Without structure or other regularization, the number of dimensions will soon overtake the number of observations, so the parameters are allowed to inform each other through a latent modeling structure known as a "multilevel model," with separate likelihoods for both the cell-wise response and the cross-cell variation.

Instead of just estimating the base response levels, I also estimate causal parameters of interest in a fashion that scales properly with dimension. The variation of the average level of an outcome around a mean is often used to smooth correlated noise across cells, and sometimes is used to directly extrapolate or forecast. In this case I don't just estimate the average level (an intercept term), but I also estimate the causal parameter of interest.

The basic approach is to estimate the first stage without variation (though you could, in principle, alter this as in [this paper](https://arxiv.org/pdf/1709.01577.pdf)). Instead of running the typical second stage regression, I split the space of the covariates of interest (call them **X**) into an even grid or lattice. To simplify computation, I presume that all of the prior information we might have about the estimate for a specific cell is given by its neighbors on that lattice. Whether this assumption is reasonable is probably application specific, but the general framework can be expanded to include non-local effects. But for now, let's pretend that this is good enough. 

In order to avoid massive computational hurdles, I use an [Intrinsic Conditional Autoregressive prior](https://mc-stan.org/users/documentation/case-studies/icar_stan.html) for the variation of the treatment effects and intercepts alon this lattice. What this implies is that if we don't see any information about the estimate for a given cell, we would presume that the estimate is an average of its neighbors. This variation is parameterized _around_ an average treatment effect across all cells. A helpful way to think about this is that in one dimension, we would estimate an average, and then we would have a "random walk" fixed to mean 0 around that average as a prior on the variation around that mean. That random walk gets updated as it sees more data. In principle you could also have an "exact" autoregressive prior where you enforce stationarity and estimate coefficients for the magnitude of the variation in each direction of the lattice, but this induces a number of computational hurdles in practice, and it's unclear that we would expect treatment effects in this context to be stationary.

The model and simulations in this context work through a fairly standard causal inference problem: an experiment with a binary instrument and a binary treatment. This is to demonstrate the method on a well-understood problem that comes up a lot in a variety of fields.

Why is the prior necessary? Let's imagine we have an instrumental variable with a dataset of 10,000 individuals. If we even have 5 covariates of interest which we want to break up into 10 groups each, we have 10 times more cells than observations! Without the prior, the model is not identified. Even in the case where we only have 3 covariates, the average cell will only have 10 people in it. If we just estimated separate 2-stage-least-squares estimates for each cell, the variance would be almost infinite.

Why wouldn't we just induce sparsity via some kind of penalized selection? This might be a good approach, but one thing to keep in mind is that we might want to understand smooth variation, rather than just a set of sparse estimates. We might expect there to be small variation across cells on a fine grid. This process also fully reflects the uncertainty in the estimation, and is fully transparent in the method by which the covariate partition happens.


# Estimation Method

I estimate outcomes via 2SLS, the just-identified case. In this context, an experiment with full compliance would just be a special case of non-compliance, so it seems fine. I estimate the first stage (for speed) via linear regression with robust standard errors, and then take the predictive means and standard errors for each observation and plug them in to the stan file.

Because I am interested in the full distribution of the treatment effects, and the high-dimensional nature of the problem, I implement the model in [Stan](https://mc-stan.org/), a probabilistic programming language that facilitates full Bayesian inference. Stan allows me to fit the model with an adaptive variant of Hamiltonian Monte Carlo sampling, which is well-suited to high-dimensional spaces. The posterior draws allow me to compute all of the summary statistics we normally like, like expectations, standard errors, predictive intervals, and the like.

For a set of ordered observable covariates (e.g. numeric or ordinal), I split them into an even lattice across all dimensions, and estimate the model with a random-walk prior on all variables. The Stan file estimates the cell-wise treatment effect, which is estimated as the average treatment effect across cells plus the idiosyncratic cell-wise component. To extrapolate to a new population, or calculate a LATE for this population, you would simply weight the various cells by the population weights of interest.

In constructing and coding this stan file I referred extensively to [this excellent paper](https://www.sciencedirect.com/science/article/abs/pii/S1877584518301175). If you want a non-paywalled variant, there is a [very similar guide here](https://mc-stan.org/users/documentation/case-studies/icar_stan.html). While the Stan code is legible, it uses some tricks for sampling efficiency that aren't the most obvious (like having a non-centered parameterization for the variation around the mean, or pairwise differences for the prior).

# Requirements

You will need R, as well as the relevant packages:

1. data.table
2. magrittr
3. estimatr
4. rstan
5. ggplot2

You will also need Stan installed, so follow the relevant install guide, which you can find in the [RStan "Getting Started" guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)

You'll also need a C++ compiler, but that is covered in the RStan getting started wiki. If you are on windows, RTools is a good option.

# Code Structure

`src/stan` includes the models coded in Stan. The main one is called `full-causal-model.stan`, and the others are either special cases or related code written by others that I used to fix some of my syntax and/or make my implementation more efficient.

`src/R` contains two useful subdirectories, `simulations` and `helpers`. 

`R/helpers`:
The `helpers` directory has scripts that are used to construct a lattice given a set of covariates, maintaining the ability to map between the location in the X-covariate cell and a vector of ids for those locations. It also has functions for constructing a list of data in the format the stan code expects, and scripts for visualizing the effects.

`R/simulations`:
This is where the real work happens! The computational experiments are in `2sls-late-sim-multiple-dimensions.R`. 
