.libPaths("~/Rlibs")
library("parallel")
library(coda)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

devtools::load_all("~/exPG/xpgTau")

library(xpgTau)

source("MJPsystems.r")
source("MJPforwards.r")
source("MJPutils.r")
source("MJPLanya.R")

load("autoreg_data.RData")

vals="Lanya"

if (vals=="GW05") { ## GW05 Biometrics
  thetas<-c(0.1,0.7,0.35,0.2,0.1,0.9,0.3,0.1) 
  x0<-c(5,8,8,8,2)
  Tend=50
}
if (vals=="WilkBk") { ## Darren's book
  thetas<-c(1,10,.01,10,1,1,.1,.1) ## Wilkinson Book
  x0<-c(0,0,0,0,10)
  Tend=200
}
if (vals=="SG22") { ## Sherlock and Golightly 2022, JCGS
  x0=c(5,5,5,1,1) 
  thetas<-c(0.1,0.7,0.35,0.2,0.1,0.9,0.3,0.1)
  Tend=50
}

if (vals=="Lanya") {
  thetas<-c(0.1,0.7,0.35,0.2,0.1,0.9,0.3,0.1) *2
  x0<-c(5,8,8,8,2) ## NB, so DNA = 8+2-DNA.P2
  # Tend=50
  Tend=20
}

tau <- 0.2
tau.times <- (1:round(Tend/tau))*tau
n <- length(tau.times)

A <- Autoreg.S
loglike <- loglike_autoreg
obs.times <- ts
time_steps <- round((obs.times/tau) + 1)

u <- nrow(A)
v <- ncol(A)

######## run xPG ############
X.initial <- X.ref
N.initial <- N.ref

thetas.true <- thetas # truth

# initial values of thetas
thetas.initial <- c(0.5, 1, 0.3, 0.8, 0.5, 2, 0.3, 0.1)

# Give theta1 and theta5 tight priors
# theta1.true = 0.2, theta5.true = 0.2
shape_prior <- c(100, 1, 1, 1, 100, 1, 1, 1)
rate_prior <- c(500, 1, 1, 1, 500, 1, 1, 1)


# M_list <- c(200, 500)
M_list <- c(1000)
logdelta <- seq(-2.4, -0, by = 0.2) # length = 13
delta_list <- signif(10**logdelta, 3)

B <- 5000
num <- 10^5

# Generate all combinations of M and delta
param_grid <- expand.grid(M = M_list, delta = delta_list)

results_list <- mclapply(1:nrow(param_grid), function(i) {
  M <- param_grid$M[i]
  delta <- param_grid$delta[i]
  
  xPG_rcpp(num = num, 
           delta = delta, 
           M = M, 
           A = A, 
           thetas.initial = thetas.initial, 
           tau = tau, 
           x0 = x0, 
           X.initial = X.initial, 
           N.initial = N.initial, 
           loglike = loglike, 
           ys = ys, 
           obs.times = obs.times, 
           shape_prior = shape_prior, 
           rate_prior = rate_prior)
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
  acc_rates_val <- colMeans(result$Z_val)
  acc_rates_ref <- colMeans(result$Z_ref)
  names(acc_rates_val) <- paste0("alpha_val", 0:L)
  names(acc_rates_ref) <- paste0("alpha_ref", 0:L)
  Theta_means <- rowMeans(result$Thetas[ , B:(num+1)])
  Theta_sd <- apply(result$Thetas[ , B:(num+1)], 1, sd)
  names(Theta_means) <- paste0("theta_means", 1:v)
  names(Theta_sd) <- paste0("theta_sd", 1:v)
  na_mean_val <- mean(result$num_NA)
  
  # Compute ESS for X0 to X50 (i.e., Paths[, 1:51, ])
  ess_list_x3 <- apply(result$Paths[3, time_steps, ], 1, effectiveSize)
  names(ess_list_x3) <- paste0("ESS_X3_t", 0:L)
  ess_list_x4 <- apply(result$Paths[4, time_steps, ], 1, effectiveSize)
  names(ess_list_x4) <- paste0("ESS_X4_t", 0:L)
  ess_list_mean <- (ess_list_x3 + ess_list_x4)/2
  names(ess_list_mean) <- paste0("ESS_mean_t", 0:L)
  
  # Combine into a data frame row
  data.frame(M = M_val,
             delta = delta_val,
             na_mean = na_mean_val,
             as.list(acc_rates_val),
             as.list(acc_rates_ref),
             as.list(Theta_means),
             as.list(Theta_sd),
             as.list(ess_list_x4),
             as.list(ess_list_mean)
  )
}))

# Ensure df is a proper data frame
df <- as.data.frame(df)
rownames(df) <- NULL

save(df, file = "autoreg_df.RData")




