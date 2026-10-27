.libPaths("~/Rlibs")
library("parallel")
library(coda)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

Rcpp::compileAttributes("~/exPG/xpgTau")

devtools::load_all("~/exPG/xpgTau")

source("SIS_functions_cpp.R")

load("SIS_high3_tau_0.2_data.RData")
# data 
# lambda = 1.2 
# mu = 0.6

D <- 3
A_D <- cbind(
  diag(D),      # infection: I_d -> I_d + 1
  -diag(D)      # recovery:  I_d -> I_d - 1
)

lambda = 1.2 # 0.08*10 
mu.true = 0.6 # 0.06*10
Ns <- rep(1000, D)
prob.initial <- 0.05
I.0 <- Ns*prob.initial # initial susceptible number of the ref path
X.0 <- I.0
tau <- 0.2
Tend <- 20

C <- 10
loglike = make_loglike(w = C)


M_list <- c(1, 5, 10, 30, 50, 70, 100)
logdelta <- seq(-2.5, -0.1, by = 0.2)
delta_list_1 <- 10**logdelta # length=13
delta_list_2 <- c(1) # length=1

delta_list <- c(delta_list_1, delta_list_2) # length=14
delta_list <- round(delta_list, 5)

# Generate all combinations of M and delta
param_grid <- expand.grid(M = M_list, delta = delta_list)

thetas.initial = c(1, 1)
shape_prior = c(1, 1)
rate_prior = c(1, 1)
num = 10^5
thin_interval = 1
B <- 10^4
num_thin <- round(num/thin_interval)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

results_list <- mclapply(1:nrow(param_grid), function(i) {
  tryCatch({
    M <- param_grid[["M"]][i]
    delta <- param_grid[["delta"]][i]
    
    # Optional: log progress
    message(sprintf("Running i = %d with M = %s, delta = %s", i, M, delta))
    
    xPG_SIS(
      num = num, 
      delta = delta, 
      M = M, 
      A = A_D, 
      thetas.initial = thetas.initial, 
      tau = tau, 
      prob.initial = prob.initial, 
      X.initial = X.ref, 
      N.initial = N.ref, 
      loglike = loglike, 
      ys = ys, 
      obs.times = obs.times, 
      shape_prior = shape_prior, 
      rate_prior = rate_prior, 
      Ns = Ns, 
      thin_interval = thin_interval
    )
  }, error = function(e) {
    message(sprintf("Error in i = %d: %s", i, e$message))
    return(NULL)
  })
}, mc.cores = parallel::detectCores() - 1, mc.set.seed = TRUE)
# Optionally, store results in a structured format
names(results_list) <- paste0("M_", param_grid$M, "_delta_", round(param_grid$delta, 5))

df <- do.call(rbind, lapply(names(results_list), function(name) {
  result <- results_list[[name]]
  if (is.null(result)) return(NULL)
  
  # Extract M and delta
  parts <- strsplit(name, "_")[[1]]
  M_val <- as.numeric(parts[2])
  delta_val <- as.numeric(parts[4])
  
  # Acceptance rates
  acc_rates_ref <- colMeans(result$Z_ref)
  names(acc_rates_ref) <- paste0("alpha_ref", obs.times)
  
  # Theta summaries
  Theta_means <- rowMeans(result$Thetas[, B:(num_thin + 1)])
  Theta_sd <- apply(result$Thetas[, B:(num_thin + 1)], 1, sd)
  Theta_ess <- apply(result$Thetas, 1, effectiveSize)
  
  n_theta <- nrow(result$Thetas)
  names(Theta_means) <- paste0("theta_means", 1:n_theta)
  names(Theta_sd) <- paste0("theta_sd", 1:n_theta)
  names(Theta_ess) <- paste0("theta_ess", 1:n_theta)
  
  na_mean_val <- sum(result$num_NA)
  
  # ESS for all D latent processes
  ess_lists <- lapply(1:D, function(k) {
    ess_k <- apply(result$Paths[k, , ], 1, effectiveSize)
    names(ess_k) <- paste0("ESS_X", k, "_t", obs.times)
    ess_k
  })
  
  # Combine automatically for any D
  ess_all <- unlist(ess_lists, use.names = TRUE)
  
  data.frame(
    M = M_val,
    delta = delta_val,
    na_mean = na_mean_val,
    as.list(acc_rates_ref),
    as.list(Theta_means),
    as.list(Theta_sd),
    as.list(Theta_ess),
    as.list(ess_all),
    check.names = FALSE
  )
}))

# Ensure df is a proper data frame
df <- as.data.frame(df)
rownames(df) <- NULL

save(df, file = "SIS_high3_final_df.RData")


