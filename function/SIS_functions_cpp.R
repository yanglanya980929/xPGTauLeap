SIS.indep.Rates <-function(r,xs,thetas, Ns) {
  # xs: length D*u
  # couple.A: D*D
  # N: length D
  # rat: length D
  D <- round(length(xs)/3)
  index_S <- 1:D
  index_I <- (D+1):(2*D)
  R1 = thetas[1]*xs[index_S]*xs[index_I] # (M+1)*D
  R2 = thetas[2]*xs[index_I]
  
  for (d in 1:D) {
    if (r == d) {
      rat <- R1[d] / Ns[d]
    } else if (r == D + d) {
      rat <- R2[d]
    } 
  }
  return (rat)
}

tau_leap_SIS_multi <- function(A, lambda, mu, xs.initial, Tend, tau, Ns) {
  u <- nrow(A)
  v <- ncol(A)
  D <- u
  index_inf <- 1:D
  index_rec <- (D + 1):(2 * D)
  n <- round(Tend/tau)
  
  X <- matrix(NA_real_, nrow = D, ncol = n + 1)
  X[, 1] <- xs.initial
  
  Rs <- matrix(NA_real_, nrow = 2 * D, ncol = n)
  
  flag_upper <- 0
  flag_lower <- 0
  for (k in 1:n) {
    I <- X[, k]
    S <- Ns - I
    # Frozen tau-leap rates
    rate_inf <- lambda * S * I / Ns
    rate_rec <- mu * I
    
    rs <- c(
      rpois(D, lambda = rate_inf * tau),
      rpois(D, lambda = rate_rec * tau)
    )
    
    r_inf <- rs[index_inf]
    r_rec <- rs[index_rec]
    
    # Need:
    # 0 <= I + r_inf - r_rec <= Ns
    # Upper boundary: I_new <= Ns
    bad_upper <- r_inf > S + r_rec
    if (any(bad_upper)) {
      r_inf[bad_upper] <- S[bad_upper] + r_rec[bad_upper]
      flag_upper <- flag_upper + sum(bad_upper)
    }
    # Lower boundary: I_new >= 0
    bad_lower <- r_rec > I + r_inf
    if (any(bad_lower)) {
      r_rec[bad_lower] <- I[bad_lower] + r_inf[bad_lower]
      flag_lower <- flag_lower + sum(bad_lower)
    }
    rs[index_inf] <- r_inf
    rs[index_rec] <- r_rec
    
    # Store the corrected reaction counts
    Rs[, k] <- rs
    # State update
    X[, k + 1] <- X[, k] + drop(A %*% rs)
  }
  times <- tau * (0:n)
  return(list(
    X = X,
    times = times,
    Rs = Rs,
    flag_upper = flag_upper,
    flag_lower = flag_lower
  ))
}

systematic_sampling <- function(weights, size) {
  # arguments: weights must be normalised weights
  # size: the number of samples
  M <- length(weights)
  cumulative_sum <- cumsum(weights)
  step <- 1 / size
  start <- runif(1, min = 0, max = step)
  positions <- seq(start, 1, by = step)
  new_samples <- integer(size)
  i <- 1
  for (j in 1:(size)){
    while (positions[j] > cumulative_sum[i]){
      i <- i + 1
    }
    new_samples[j] <- i
  }
  return(new_samples)
}

loglike <- function(x, y, w){
  probs <- dnorm(c(1:(2*w+1), (2*w):1), mean = 2*w+1, sd = w, log = FALSE)
  probs <- probs/sum(probs)
  if(abs(x-y)<=2*w){
    L <- probs[round(2*w+1-abs(x-y))]
  }
  else{
    L <- 0
  }
  return(log(L))
}

make_loglike <- function(w) {
  probs <- dnorm(c(1:(round(2*w)+1), round(2*w):1), mean = 2*w+1, sd = w, log = FALSE)
  probs <- probs/sum(probs)
  function(xt, y) {
    d <- abs(xt - y)
    if (all(d <= 2*w)) {
      return(sum(log(probs[round(2*w + 1 - d)])))
    } else {
      return(-Inf)
    }
  }
}


trace_back <- function(obs.times, tau, M, X.ref, con_SIS) {
  time_steps <- round(obs.times / tau) + 1
  L <- length(obs.times)
  D <- nrow(X.ref)
  n <- ncol(X.ref)
  stopifnot(
    time_steps[1] == 1,
    time_steps[L] == n,
    dim(con_SIS$X.proposals)[1] == M + 1,
    dim(con_SIS$X.proposals)[2] == D
  )
  # One reconstructed path for each of the M+1 final particles
  Paths <- array(NA_real_, dim = c(M + 1, D, n))
  for (k in 1:(M + 1)) {
    index <- k
    # Last observation interval
    k.2 <- n
    k.1 <- time_steps[L - 1] + 1
    Paths[k, , k.1:k.2] <- matrix(
      con_SIS$X.proposals[index, , k.1:k.2],
      nrow = D, ncol = k.2 - k.1 + 1
    )
    # Earlier observation intervals
    if (L > 2) {
      for (j in seq_len(L - 2)) {
        index <- con_SIS$ancestors[index, L - j]
        k.2 <- k.1 - 1
        k.1 <- time_steps[L - j - 1] + 1
        Paths[k, , k.1:k.2] <- matrix(
          con_SIS$X.proposals[index, , k.1:k.2],
          nrow = D, ncol = k.2 - k.1 + 1
        )
      }
    }
    # Initial state at time 0
    index <- con_SIS$ancestors[index, 1]
    Paths[k, , 1] <- con_SIS$X.proposals[index, , 1]
  }
  return(Paths)
}



transition_prob <- function(i, j, N, lambda, mu, tau) {
  lambda_plus <- lambda * (N-i) * i * tau/N  # Mean of susceptible Poisson
  lambda_minus <- mu * i * tau    # Mean of recovery Poisson
  # Sum over k (Poisson convolution)
  prob <- 0
  if(j==N){
    for(new in N:(N+100)){
      j <- new
      for (k in 0:200) {  # Truncate sum for computational efficiency
        if (k + j - i >= 0) {
          prob <- prob + dpois(k + j - i, lambda_plus) * dpois(k, lambda_minus)
        }
      }
    }
  }
  if(j==0){
    for(new in -100:0){
      j <- new
      for (k in 0:200) {  # Truncate sum for computational efficiency
        if (k + j - i >= 0) {
          prob <- prob + dpois(k + j - i, lambda_plus) * dpois(k, lambda_minus)
        }
      }
    }
  }
  else{
    for (k in 0:200) {  # Truncate sum for computational efficiency
      if (k + j - i >= 0) {
        prob <- prob + dpois(k, lambda_minus) * dpois(k + j - i, lambda_plus)
      }
    }
  }
  return(prob)
}

# multi/single SIS xPG;
conditional_SimSIS <- function(A, tau, M, X.ref, N.ref, delta, loglike, ys, 
                               obs.times, prob.initial, thetas, Ns){
  n <- ncol(N.ref)
  if (is.null(dim(ys))) {
    ys <- matrix(ys, nrow = 1)
  }
  L <- ncol(ys)
  u <- nrow(A) # D
  v <- ncol(A) # 2*D
  
  D <- u
  
  N.proposals <- array(NA, dim = c(M+1, v, n))
  N.proposals[1, , ] <- N.ref

  I.0 <- X.ref[ , 1] # 10 dim 
  
  # generate M exchangeable initial states
  N.0 <- floor(Ns*(1-delta))
  I.0.proposals <- matrix(NA, nrow = M, ncol = u)
  for(i in 1:D){
    I.core <- rhyper(1, I.0[i], Ns[i]-I.0[i], N.0[i])
    I.0.proposals[ , i] <- I.core + rbinom(M, Ns[i]-N.0[i], prob = prob.initial)
  }
  
  X.proposals <- array(NA, dim = c(M+1, u, n+1)) # (M+1)*D*201-dim
  X.proposals[1, , ] <- X.ref # 201
  X.proposals[2:(M+1), , 1] <- matrix(I.0.proposals, nrow = M, ncol = u)
  
  time_steps <- round((obs.times/tau) + 1) # includes 1
  
  stopifnot(
    time_steps[1] == 1,      # first observation at t = 0
    all(diff(time_steps) > 0),
    time_steps[L] == n + 1   # if final observation = end of trajectory
  )
  
  ancestors <- matrix(NA, nrow = M+1, ncol = L-1) # (M+1)*(L-1)
  ancestors[1,] <- 1
  
  X.ancestors <- matrix(NA, nrow = M+1, u) # we still need the ancestor of the reference to get mu^{(0)}

  # do resampling at time zero
  logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, , 1], ys[, 1]), numeric(1))
  logW.max <- max(logW)
  wstars <- exp(logW - logW.max) # length M
  y.loglikely <- log(mean(wstars)) + logW.max # maximum marginal loglikelhood
  w.nor <- wstars/sum(wstars) 
  indexes <- sample(1:(M+1), M, prob = w.nor, replace = TRUE) # use multinomial sampling 
  ancestors[2:(M+1),1] <- indexes
  X.ancestors[1,] <- X.proposals[1, , 1]
  X.ancestors[2:(M+1),] <- X.proposals[indexes, , 1]
  for(l in 1:(L-1)){
    # propagate 
    N.ref_local <- N.ref[ , time_steps[l]:(time_steps[l+1]-1), drop = FALSE] 
    X.ref_local <- X.ref[ , (time_steps[l]+1):time_steps[l+1], drop = FALSE] 
    propagation <- propagate_rcpp_SIS(l, X.ancestors, X.ref_local, N.ref_local, delta, M, time_steps, A, Ns, tau, thetas)
    
    X.proposals[ , , (time_steps[l]+1):time_steps[l+1]] <- propagation$X.proposals_local # (M+1)*u*4
    N.proposals[ , , time_steps[l]:(time_steps[l+1]-1)] <- propagation$N.proposals_local # (M+1)*v*4
    # resample; update X.ancestors
    if(l < (L-1)){
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, , time_steps[l+1]], ys[ , l+1]), numeric(1))
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      w.nor <- wstars/sum(wstars)
      indexes <- sample(1:(M+1), M, prob = w.nor, replace = TRUE) # use multinomial sampling 
      ancestors[2:(M+1), l+1] <- indexes
      X.ancestors[1,] <- X.proposals[1, , time_steps[l+1]]
      X.ancestors[2:(M+1),] <- X.proposals[indexes, , time_steps[l+1]]
    }else{
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, , time_steps[L]], ys[ , L]), numeric(1))
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      w.nor <- wstars/sum(wstars)
    }
  }
  na_num <- sum(apply(X.proposals[ , , time_steps], 1, function(row) any(is.na(row))))
  return(list(weights=w.nor, X.proposals=X.proposals, N.proposals=N.proposals, ancestors=ancestors, na_num=na_num))
}

update_theta_SIS <- function(vv, N.path, path, shape_prior, rate_prior, tau, Ns) {
  n <- ncol(N.path)
  D <- length(Ns)
  
  stopifnot(
    nrow(path) == D,
    ncol(path) == n + 1,
    nrow(N.path) == 2 * D,
    vv %in% 1:2
  )
  # ------------------------------------------------------------
  # States at the beginning of each tau interval
  # ------------------------------------------------------------
  I <- path[, 1:n, drop = FALSE]
  S <- matrix(Ns, nrow = D, ncol = n) - I
  
  # ------------------------------------------------------------
  # Integrated rate multipliers
  #
  # Infection: lambda_d = beta * S_d * I_d / N_d
  #
  # Recovery: lambda_d = gamma * I_d
  # ------------------------------------------------------------
  rate_beta <- sum((S * I) /matrix(Ns, nrow = D, ncol = n))
  rate_gamma <- sum(I)
  
  index_inf <- 1:D
  index_rec <- (D + 1):(2 * D)
  
  count_beta <- sum(N.path[index_inf, , drop = FALSE])
  count_gamma <- sum(N.path[index_rec, , drop = FALSE])
  if (vv == 1) {
    theta_new <- rgamma(
      1,
      shape = shape_prior[1] + count_beta,
      rate  = rate_prior[1] + tau * rate_beta
    )
  } else {
    theta_new <- rgamma(
      1,
      shape = shape_prior[2] + count_gamma,
      rate  = rate_prior[2] + tau * rate_gamma
    )
  }
  return(theta_new)
}

xPG_SIS <- function(num, delta, M, A, thetas.initial, tau, prob.initial, X.initial, N.initial, loglike, 
                    ys, obs.times, shape_prior, rate_prior, Ns, thin_interval){
  num_thin <- floor(num / thin_interval)
  
  thetas <- thetas.initial
  num_NA <- rep(NA, num)
  n <- ncol(N.initial) # 200
  if (is.null(dim(ys))) {
    ys <- matrix(ys, nrow = 1)
  }
  L <- ncol(ys)
  u <- nrow(A)
  v <- ncol(A)
  D <- u  
  Thetas <- matrix(NA, nrow = length(thetas.initial), ncol = num_thin + 1)
  Thetas[ , 1] <- thetas.initial 
  
  time_steps <- c(round((obs.times/tau) + 1)) #indexes of ys in X
  Paths <- array(NA, dim=c(u, L, num_thin+1))
  Paths[ , ,1] <- X.initial[ , time_steps]
  # Paths <- array(NA, dim=c(u, n+1, num_thin+1))
  # Paths[ , ,1] <- X.initial
  
  N.path <- N.initial
  path <- X.initial
  Z_ref <- matrix(NA, nrow = num_thin, ncol = L) # for alpha_ref
  for (i in 1:num){
    # simulate a path based on lambda
    conPF <- conditional_SimSIS(A, tau, M, path, N.path, delta, loglike, ys, 
                                obs.times, prob.initial, thetas, Ns)
    weights <- conPF$weights
    ancestors <- conPF$ancestors
    X.proposals <- conPF$X.proposals
    N.proposals <- conPF$N.proposals
    num_NA[i] <- conPF$na_num
    k <- sample(1:(M+1), 1, prob = weights, replace = TRUE) # multinomial sampling
    
    ## obtain the new path back from k
    path <- matrix(NA, nrow = u, ncol = n+1)
    N.path <- matrix(NA, nrow = v, ncol = n)
    k.2 <- n+1
    k.1 <- time_steps[length(time_steps)-1]+1
    Index <- rep(NA, L)
    index <- k
    Index[L] <- index
    path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u, ncol = k.2-k.1+1)
    N.path[ , (k.1-1):(k.2-1)] <- matrix(N.proposals[index, , (k.1-1):(k.2-1)], nrow = v, ncol = k.2-k.1+1)
    if (L > 2) {
      for (j in seq_len(L - 2)) {
        index <- ancestors[index, L - j]
        Index[L - j] <- index
        k.2 <- k.1 - 1
        k.1 <- time_steps[L - j - 1] + 1
        path[, k.1:k.2] <- matrix(X.proposals[index, , k.1:k.2], nrow = u, ncol = k.2 - k.1 + 1)
        N.path[, (k.1 - 1):(k.2 - 1)] <- matrix(N.proposals[index, , (k.1 - 1):(k.2 - 1)],
                                                nrow = 2 * u, ncol = k.2 - k.1 + 1)
      }
    }
    index <- ancestors[index,1]
    Index[1] <- index
    path[ , 1] <- matrix(X.proposals[index, ,1], nrow = u, ncol = 1)
    # update all the thetas 
    for(vv in 1:2){
      theta_new <- update_theta_SIS(vv, N.path, path, shape_prior, rate_prior, tau, Ns)
      thetas[vv] <- theta_new
    }
    if(i %% thin_interval == 0){
      i_thin <- floor(i/thin_interval)
      Paths[, , i_thin+1] <- path[ , time_steps]
      z_ref <- ifelse(Index != 1, 1, 0) # for computing alpha_ref
      Z_ref[i_thin, ] <- z_ref
      Thetas[ , i_thin+1] <- thetas
    }
  }
  return(list(Paths=Paths, Z_ref=Z_ref, Thetas=Thetas, num_NA=num_NA, N.path=N.path, path=path))
}

