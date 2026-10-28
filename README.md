# xPGTauLeap

## Introduction
This repository contains code and data accompanying the paper
*Scalable Bayesian inference for time-discretised Markov jump processes
via exchangeable particle Gibbs*

The code implements exchangeable particle Gibbs (xPG) and particle Gibbs
for Bayesian inference in time-discretised Markov jump processes.

## Repository contents

- `xpgTau`: Rcpp package implementing the algorithms.
- `data`: the true latent trajectory and the corresponding observations.
- `execution`: R scripts for running the experiments.
- `functions`: R functions implementing the algorithms.
- `simulation_results`: experiment results of running xPG and xPGBS for SIS, SIR, SEIR and autoregulatory models.
- `generate_figures.R`: Code for producing the figures.

## Requirements

The code requires R and the following packages:

- xpgTau
- Rcpp
- RcppArmadillo
- parallel

## Installation

Install the required R packages:

```r
install.packages(c("Rcpp", "RcppArmadillo", "parallel"))

Rcpp::compileAttributes("~/xpgTau")

devtools::load_all("~/xpgTau")



