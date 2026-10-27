# xPGTauLeap
Data and codes for exchangeable particle Gibbs for time-discretised Markov jump processes


# Exchangeable particle Gibbs for time-discretised Markov jump processes

This repository contains code and data accompanying the paper
*Scalable Bayesian inference for time-discretised Markov jump processes
via exchangeable particle Gibbs*

The code implements exchangeable particle Gibbs (xPG) and particle Gibbs
for Bayesian inference in time-discretised Markov jump processes.

## Repository contents

- `[xpgTau]/`: Rcpp package implementing the algorithms.
- `[scripts-folder]/`: R scripts for running the experiments.
- `[figures-folder]/`: Code for producing the figures.

## Requirements

The code requires R and the following packages:

- Rcpp
- RcppArmadillo
- [Add any other required packages.]

Building the Rcpp package from source also requires a suitable C++
compiler.

## Installation

Install the required R packages:

```r
install.packages(c("Rcpp", "RcppArmadillo", "remotes"))
