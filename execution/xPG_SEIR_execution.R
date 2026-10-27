.libPaths("~/Rlibs")
library("parallel")
library(coda)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(extraDistr)

Rcpp::compileAttributes("~/exPG/xpgTau")

devtools::load_all("~/exPG/xpgTau")

source("SEIR_functions_cpp.R")

load("SEIR_high4_eta_0.6_tau_0.1_data.RData")

D <- 4
A_D <- matrix(0, nrow = 4*D, ncol=3*D) 
for(d in 1:D){
  for(j in 1:3){
    col_index <- (j-1)*D+d
    A_D[c(1, D+1)+col_index-1 , col_index] <- c(-1, 1)
  }
}

# across infection rates
SEIR.coupled.A <- matrix(c(1, 0.2, 0.1, 0.1, 0.2, 1, 0.15, 0.2, 0.1, 0.15, 1, 0.1, 0.1, 0.2, 0.1, 1), nrow = 4, ncol = 4)


x0 <- c(890, 840, 870, 880, 80, 100, 90, 90, 30, 60, 40, 30, 0, 0, 0, 0)

Ns <- rep(NA, D)
for(i in 1:D){
  Ns[i] <- sum(x0[c(i,D+i,2*D+i,3*D+i)])
}

# thetas <- SEIR.thetas.true
coupled.A <- SEIR.coupled.A
A <- A_D
Tend <- 20
L <- 20
tau <- 0.1

n <- round(Tend/tau)
tau.times <- (1:round(Tend/tau))*tau

index_S <- 1:D
index_E <- (D+1):(2*D)
index_I <- (2*D+1):(3*D)
index_R <- (3*D+1):(4*D)

obs.times <- 1:Tend
obs_index <- round(obs.times/tau)+1
time_steps <- c(1, obs_index) # including the initial state

loglike <- function(xs, y) {
  loglike.bin(xs, y, p = 0.9)
}

thetas.initial <- c(0.5, 1, 0.8)
# truth: 0.8, 0.6, 0.2
shape_prior <- c(1, 1, 20)
rate_prior <- c(1, 1, 100)


M_list <- c(1, 3, 5, 10)
# M_list <- c(130, 150, 180, 200, 250)
# M_list <- c(300)
delta_list <- c(0.01, 0.02, 0.05, (1:5)*0.1, 0.7, 1) # length=10


# Generate all combinations of M and delta
param_grid <- expand.grid(M = M_list, delta = delta_list)

B <- 10^4
num <- 10^6
thin_interval <- 10
num_thin <- round(num/thin_interval)

RNGkind("L'Ecuyer-CMRG")
set.seed(3)
# set.seed(12345)

results_list <- mclapply(1:nrow(param_grid), function(i) {
  tryCatch({
    M <- param_grid[["M"]][i]
    delta <- param_grid[["delta"]][i]
    
    # Optional: log progress
    message(sprintf("Running i = %d with M = %s, delta = %s", i, M, delta))
    
    xPG_SEIR_rcpp(num = num, 
                  delta = delta, 
                  M = M, 
                  A = A,
                  thetas.initial = thetas.initial, 
                  tau = tau, 
                  x0 = x0,
                  X.initial = X.ref, 
                  N.initial = N.ref, 
                  loglike = loglike, 
                  ys = ys, 
                  obs.times = obs.times, 
                  shape_prior = shape_prior, 
                  rate_prior = rate_prior, 
                  Ns = Ns, 
                  coupled.A = coupled.A,
                  thin_interval = thin_interval
    )
  }, error = function(e) {
    message(sprintf("Error in i = %d: %s", i, e$message))
    return(NULL)
  }
  )
}, mc.cores = parallel::detectCores() - 1)
# Optionally, store results in a structured format
names(results_list) <- paste0("M_", param_grid$M, "_delta_", param_grid$delta)

df <- do.call(rbind, lapply(names(results_list), function(name) {
  result <- results_list[[name]]
  if (is.null(result)) return(NULL)  # Skip if the result is NULL due to errors
  
  # Extract M and delta from the name string
  parts <- strsplit(name, "_")[[1]]
  M_val <- as.numeric(parts[2])  # Extract M value
  delta_val <- as.numeric(parts[4])  # Extract delta value
  
  # Compute acceptance rates
  acc_rates_ref <- colMeans(result$Z_ref)
  names(acc_rates_ref) <- paste0("alpha_ref", 1:L)
  Theta_means <- rowMeans(result$Thetas[ , B:(num_thin+1)])
  Theta_sd <- apply(result$Thetas[ , B:(num_thin+1)], 1, sd)
  names(Theta_means) <- paste0("theta_means", 1:3)
  names(Theta_sd) <- paste0("theta_sd", 1:3)
  Theta_ess <- apply(result$Thetas, 1, effectiveSize)
  names(Theta_ess) <- paste0("theta_ess", 1:3)
  na_mean_val <- sum(result$num_NA)
  
  # Compute ESS for X0 to X50 (i.e., Paths[, 1:51, ])
  ess_lists <- list()
  for (k in 1:(4*D)) {
    ess_k <- apply(result$Paths[k, 2:(L+1), ], 1, effectiveSize)
    names(ess_k) <- paste0("ESS_X", k, "_t", 1:L)
    ess_lists[[k]] <- ess_k
  }
  
  # Combine into a data frame row
  df <- data.frame(M = M_val,
                   delta = delta_val,
                   na_mean = na_mean_val,
                   as.list(acc_rates_ref),
                   as.list(Theta_means),
                   as.list(Theta_sd),
                   as.list(Theta_ess),
                   as.list(ess_lists[[1]]),  # X1
                   as.list(ess_lists[[2]]),  # X2
                   as.list(ess_lists[[3]]),  # X3
                   as.list(ess_lists[[4]]),  # X4
                   as.list(ess_lists[[5]]),  # X5
                   as.list(ess_lists[[6]]),  # X6
                   as.list(ess_lists[[7]]),  # X7
                   as.list(ess_lists[[8]]),  # X8
                   as.list(ess_lists[[9]]),  # X9
                   as.list(ess_lists[[10]]), # X10
                   as.list(ess_lists[[11]]), # X11
                   as.list(ess_lists[[12]]), # X12
                   as.list(ess_lists[[13]]),
                   as.list(ess_lists[[14]]),
                   as.list(ess_lists[[15]]),
                   as.list(ess_lists[[16]])
  )
}))

# Ensure df is a proper data frame
df <- as.data.frame(df)
rownames(df) <- NULL

save(df, file = "SEIR_xpg_final_tau_0.1_high4.RData")
