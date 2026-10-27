.libPaths("~/Rlibs")
library("parallel")
library(coda)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
# library(devtools)

# detach("package:XPG", unload = TRUE, force = TRUE)
# devtools::load_all("/home/yangl35/STOR-601-env/xPG_MJP_cpp/XPG")
# devtools::load_all("~/exPG/XPG")

library(XPG)
setwd("/home/yangl35/STOR-601-env/XPG_MJPs")
source("SIR_functions_cpp.R")

u <- 3
v <- 2
A <- matrix(c(-1, 1, 0, 0, -1, 1), nrow = u, ncol = v)
lambda = 0.6
mu.true = 0.2
n <- 100*2
L <- 10*2
N <- 1000
prob.initial <- 0.1
I.0 <- N*prob.initial
S.0 <- N-I.0
R.0 <- 0
X.0 <- c(S.0, I.0, R.0)
tau <- 0.1
C <- 5
D <- 2

# generate the true state process and the initial path
set.seed(12345)
list.ref <- tau_leap_SIR_multi(A, lambda, mu.true, X.0, n, tau, N)
X.ref <- list.ref$X
times <- list.ref$times
N.1.ref <- list.ref$R.1 # number of R1 in each time interval
N.2.ref <- list.ref$R.2 # number of R2 in each time interval

ind <- (1:L)*round(n/L)+1 # the actual index of the states at obs times
ind <- c(1, ind)

X.obs <- X.ref[ , ind]

w <- C
probs <- dnorm(c(1:(2*w+1), (2*w):1), mean = 2*w+1, sd = w, log = FALSE)
probs <- probs/sum(probs)

set.seed(12345)
ys<- X.obs[2,] + sample((-2*w):(2*w), length(X.obs[1,]), prob = probs, replace = TRUE)
ys <- matrix(ys, nrow = 1)
obs.times <- (ind-1)*tau # ind is actually ts, we assume obs.times will be given in normal cases

# generate obs
X.obs <- X.ref[ ,ind] # the hidden states at the obs times
w <- C
probs <- dnorm(c(1:(2*w+1), (2*w):1), mean = 2*w+1, sd = w, log = FALSE)
probs <- probs/sum(probs)

# ys.S <- X.obs[1,] + sample((-2*w):(2*w), length(X.obs[1,]), prob = probs, replace = TRUE)
# ys.I <- X.obs[2,] + sample((-2*w):(2*w), length(X.obs[1,]), prob = probs, replace = TRUE)
# ys.R <- X.obs[3,] + sample((-2*w):(2*w), length(X.obs[1,]), prob = probs, replace = TRUE)
# # ys <- rbind(ys.S, ys.I, N-ys.S-ys.I)
# ys <- rbind(ys.S, ys.I, ys.R)
# rownames(ys) <- NULL

#png("/home/yangl35/h-drive/xPG_MJP/plot/SIR_truth_new11.png", width = 800, height = 600)
#ylim_range <- range(X.ref, ys)
#plot(times, X.ref[1, ], type = "l", col = "blue",
#     xlab = "time (days)", ylab = "number of people",  main = paste("lambda=",lambda, "mu=", mu.true), ylim = ylim_range)
#lines(times, X.ref[2, ], col = "red")
#lines(times, X.ref[3, ], col = "green")
#points(obs.times, ys[1,], col = "black", pch = 1)
#points(obs.times, ys[2,], col = "purple", pch = 1)
#points(obs.times, ys[3,], col = "green", pch = 1)
#dev.off()

M_list <- c(1, 5, 10, 30, 50, 70, 100)
delta_list <- c(0.01, 0.05, (1:10)*0.1) # length=12

# Generate all combinations of M and delta
param_grid <- expand.grid(M = M_list, delta = delta_list)

lambda.initial = 1
mu.initial = 1
shape_prior_lambda = 1
rate_prior_lambda = 1
shape_prior_mu = 1
rate_prior_mu = 1

X.ref <- array(X.ref, dim = c(1, 3, 201))
N.1.ref <- matrix(N.1.ref, nrow = 1)
N.2.ref <- matrix(N.2.ref, nrow = 1)

results_list <- mclapply(1:nrow(param_grid), function(i) {
  tryCatch({
    M <- param_grid[["M"]][i]
    delta <- param_grid[["delta"]][i]
    
    # Optional: log progress
    message(sprintf("Running i = %d with M = %s, delta = %s", i, M, delta))
    
    xPG_multi(
      num = 10,  # Replace with actual value if needed
      delta = delta, 
      M = M, 
      A = A,
      lambda.initial = lambda.initial, 
      mu.initial = mu.initial, 
      tau = tau, 
      X.initial = X.ref, 
      N.1.initial = N.1.ref, 
      N.2.initial = N.2.ref, 
      loglike = make_loglike(w = C), 
      ys = ys, 
      obs.times = obs.times, 
      N = N, 
      prob.initial = prob.initial, 
      C = C,
      shape_prior_lambda = shape_prior_lambda, 
      rate_prior_lambda = rate_prior_lambda, 
      shape_prior_mu = shape_prior_mu, 
      rate_prior_mu = rate_prior_mu
    )
  }, error = function(e) {
    message(sprintf("Error in i = %d: %s", i, e$message))
    return(NULL)
  })
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
  
  # Compute acc_rates and ESS of X at the first obs time point 
  acc_rates_ref <- mean(result$U)
  acc_rates_val_X0 <- colMeans(result$Z)[1]
  acc_rates_val_X1 <- colMeans(result$Z)[2]
  ESS_lambda <- effectiveSize(result$Lambda)
  ESS_mu <- effectiveSize(result$Mu)
  
  # Compute ESS for X0 to X10 (i.e., Paths[, 1:11, ])
  ess_list_S <- sapply(1:(L+1), function(i) effectiveSize(result$Paths.S[, round((i-1)*10+1), ]))
  names(ess_list_S) <- paste0("ESS_X", 0:L, "_S")
  ess_list_I <- sapply(1:(L+1), function(i) effectiveSize(result$Paths.I[, round((i-1)*10+1), ]))
  names(ess_list_I) <- paste0("ESS_X", 0:L, "_I")
  ess_list_R <- sapply(1:(L+1), function(i) effectiveSize(result$Paths.R[, round((i-1)*10+1), ]))
  names(ess_list_R) <- paste0("ESS_X", 0:L, "_R")
  
  # Combine into a data frame row
  data.frame(M = M_val,
             delta = delta_val,
             ESS_lambda = ESS_lambda,
             ESS_mu = ESS_mu,
             acc_rates_ref = acc_rates_ref,
             acc_rates_val_X0 = acc_rates_val_X0,
             acc_rates_val_X1 = acc_rates_val_X1,
             as.list(ess_list_S),
             as.list(ess_list_I),
             as.list(ess_list_R))
}))

# Ensure df is a proper data frame
df <- as.data.frame(df)
rownames(df) <- NULL

M_100_list <- results_list[grepl("^M_100", names(results_list))]

delta_list <- c(0.01, 0.05, (1:10)*0.1)
x <- 0:L
all_data <- data.frame()

for(i in seq_along(delta_list)) {
  Acc <- apply(M_100_list[[i]]$Z, 2, mean)
  temp_df <- data.frame(
    x = x,
    Acc = Acc,
    delta = as.factor(delta_list[i])  # Use factor for legend
  )
  all_data <- bind_rows(all_data, temp_df)
}

# save(results_list, file = "grid_10obs_SIR_extra1_10w.RData")
save(all_data, file = "SIR_accDecay_new12.RData")
save(df, file = "df_SIR_new12.RData")




