safe_ess <- function(x) {
  if (length(unique(x)) <= 1 || sd(x, na.rm = TRUE) == 0) {
    return(0)
  } else {
    return(effectiveSize(x))
  }
}

mean_check <- function(x1, x2){
  return(abs(mean(x1)/mean(x2) -1))
}
sd_check <- function(x1, x2){
  return(abs(sd(x1)/sd(x2) -1))
}

compute_CI <- function(x1, x2){
  mu_hat <- mean(x1)-mean(x2)  # 0.04195217
  sigma2 <- var(x1)/effectiveSize(x1)+var(x2)/effectiveSize(x2)
  sigma  <- sqrt(sigma2)   # theoretical variance
  CI <- c(mu_hat - 1.96 * sigma,
          mu_hat + 1.96 * sigma)
  return(CI)
}

SEIR.Rates.single <-function(r,xs,thetas, N) {
  # xs: length u
  R1 = thetas[1]*xs[1]*xs[3]
  R2 = thetas[2]*xs[2]
  R3 = thetas[3]*xs[3]
  if (r == 1) {
    rat <- R1/N
  } else if (r == 2) {
    rat <- R2
  } else if (r == 3) {
    rat <- R3
  }
  return (rat)
}

SEIR.Rates.coupled <-function(r,xs,thetas, coupled.A, Ns) {
  # xs: length D*u
  # couple.A: D*D
  # N: length D
  # rat: length D
  D <- round(length(xs)/4)
  index_S <- 1:D
  index_E <- (D+1):(2*D)
  index_I <- (2*D+1):(3*D)
  R1 = thetas[1]*xs[index_S]*(coupled.A %*% xs[index_I]) # (M+1)*D
  R2 = thetas[2]*xs[index_E]
  R3 = thetas[3]*xs[index_I]
  for (d in 1:D) {
    if (r == d) {
      rat <- R1[d] / Ns[d]
    } else if (r == D + d) {
      rat <- R2[d]
    } else if (r == 2 * D + d) {
      rat <- R3[d]
    }
  }
  return (rat)
}

SEIR.Rates.coupled.vec <- function(r, xs, thetas, coupled.A, Ns) {
  # xs: matrix (M+1)*(D*u)
  # couple.A: D*D
  # Ns: length D
  # rat: length (M+1)
  # r is 1,2,...,3*D
  stopifnot(ncol(xs) %% 4 == 0)
  D <- ncol(xs) / 4
  index_S <- 1:D
  index_E <- (D + 1):(2 * D)
  index_I <- (2 * D + 1):(3 * D)
  S <- xs[, index_S, drop = FALSE]
  E <- xs[, index_E, drop = FALSE]
  I <- xs[, index_I, drop = FALSE]
  
  infection_pressure <- I %*% t(coupled.A)
  
  R1 <- thetas[1] * S * infection_pressure
  R2 <- thetas[2] * E
  R3 <- thetas[3] * I
  if (r <= D) {
    d <- r
    rat <- R1[, d] / Ns[d]
  } else if (r <= 2 * D) {
    d <- r - D
    rat <- R2[, d]
  } else {
    d <- r - 2 * D
    rat <- R3[, d]
  }
  return(rat)
}

tau_leap_single <- function(A, thetas, xs.initial, n, tau, N){
  X <- matrix(NA, nrow = 4, ncol = n+1)
  X[ , 1] <- xs.initial
  Rs <- matrix(NA, nrow = 3, ncol = n)
  for(k in 1:n){
    Mu <- sapply(1:3, function(i) SEIR.Rates.single(i, X[ , k], thetas, N))*tau 
    rs <- rpois(3, lambda = as.vector(Mu))
    Rs[ , k] <- rs
    if(any(rs[1] > X[1, k])){
      rs[1] <- X[1, k]
    }
    if(any(rs[2] > X[2, k] + rs[1])){
      rs[2] <- X[2, k] + rs[1]
    }
    if(any(rs[3] > X[3, k] + rs[2])){
      rs[3] <- X[3, k] + rs[2]
    }
    X[ , k+1] <- X[ , k]+rs[1]*A[,1]+rs[2]*A[,2]+rs[3]*A[,3]
  }
  times <- (1:n)*tau
  times <- c(0, times)
  return(list(X = X, times = times, Rs = Rs))
}

tau_leap_couple <- function(A, thetas, xs.initial, n, tau, Ns, coupled.A){
  # initial_states: D*u
  u <- nrow(A)
  v <- ncol(A) # number of the types of reactions
  D <- round(u/4)
  index_S <- 1:D
  index_E <- (D+1):(2*D)
  index_I <- (2*D+1):(3*D)
  X <- matrix(NA, nrow = u, ncol = n+1)
  X[ , 1] <- xs.initial
  Rs <- matrix(NA, nrow = v, ncol = n)
  for(k in 1:n){
    Mu <- sapply(1:v, function(i) SEIR.Rates.coupled(i, X[ , k], thetas, coupled.A, Ns))*tau # D*v
    rs <- rpois(v, lambda = as.vector(Mu))
    Rs[ , k] <- rs
    if(any(rs[index_S] > X[index_S, k])){
      rs[index_S] <- X[index_S, k]
    }
    if(any(rs[index_E] > X[index_E, k] + rs[index_S])){
      rs[index_E] <- X[index_E, k] + rs[index_S]
    }
    if(any(rs[index_I] > X[index_I, k] + rs[index_E])){
      rs[index_I] <- X[index_I, k] + rs[index_E]
    }
    X[ , k+1] <- X[ , k]+rs[index_S]%*%t(A[,index_S])+rs[index_E]%*%t(A[,index_E])+rs[index_I]%*%t(A[,index_I])
  }
  times <- (1:n)*tau
  times <- c(0, times)
  return(list(X = X, times = times, Rs = Rs))
}

make_loglike <- function(w) {
  probs <- dnorm(c(1:(round(2*w) + 1), round(2*w):1),
                 mean = 2*w + 1, sd = w, log = FALSE)
  probs <- probs / sum(probs)
  
  function(xt, y) {
    d <- abs(xt - y)
    
    # check for NA values in inputs
    if (anyNA(xt) || anyNA(y)) {
      return(-Inf)
    } else if (all(d <= 2*w)) {
      return(sum(log(probs[round(2*w + 1 - d)])))
    } else {
      return(-Inf)
    }
  }
}

loglike.bin <- function(xs, y, p){
  # xs: length D
  # y: length D
  if(anyNA(xs)||anyNA(y)){
    return(-Inf)
  } else{
    l <- dbinom(y, size = xs, prob = p, log = TRUE)
    return(sum(l))
  }
}

loglike.SEIR.NBin <- function(xs, y, dispersion_para){
  # xs: length 2
  # y: length 2
  if(anyNA(xs)||anyNA(y)){
    return(-Inf)
  } else{
    l <- dnbinom(y, size = dispersion_para, mu = xs, log = TRUE)
    return(sum(l))
  }
}

update_theta_SEIR <- function(vv, N.path, path, shape_prior, rate_prior, tau, coupled.A, Ns){
  n <- ncol(N.path)
  D <- round(nrow(path)/4)
  S_index <- 1:D
  E_index <- (D+1):(2*D)
  I_index <- (2*D+1):(3*D)
  rate1 <- sum(rowSums(coupled.A %*% path[I_index, 1:n]*path[S_index, 1:n])/Ns) # D*n
  rate2 <- sum(path[E_index, 1:n])
  rate3 <- sum(path[I_index, 1:n])
  if(vv == 1){
    theta_new <- rgamma(1, shape = shape_prior[vv] + sum(N.path[S_index, ]), 
                        rate = rate_prior[vv] + tau*rate1)
  }
  if(vv == 2){
    theta_new <- rgamma(1, shape = shape_prior[vv] + sum(N.path[E_index, ]), 
                        rate = rate_prior[vv] + tau*rate2)
  }
  if(vv == 3){
    theta_new <- rgamma(1, shape = shape_prior[vv] + sum(N.path[I_index, ]), 
                        rate = rate_prior[vv] + tau*rate3)
  }
  return(theta_new)
}


# xPG
conditional_SEIR_rcpp <- function(A, tau, M, X.ref, N.ref, delta, loglike, ys, 
                                           obs.times, x0, thetas, Ns, coupled.A){
  n <- ncol(N.ref) # 200
  if (is.null(dim(ys))) {
    ys <- matrix(ys, nrow = 1)
  }
  L <- ncol(ys)
  u <- nrow(A)
  v <- ncol(A)
  D <- round(u/4)
  
  N.proposals <- array(NA, dim = c(M+1, v, n))
  N.proposals[1, , ] <- N.ref
  
  X.proposals <- array(NA, dim = c(M+1, u, n+1)) # (M+1)*4*201-dim
  X.proposals[1, , ] <- X.ref # 4*201
  X.proposals[2:(M+1), , 1] <- matrix(x0, nrow = M, ncol = u, byrow = TRUE) # M*1-dim
  
  index_S <- 1:D
  index_E <- (D+1):(2*D)
  index_I <- (2*D+1):(3*D)
  index_R <- (3*D+1):(4*D)
  
  ancestors <- matrix(NA, nrow = M+1, ncol = L-1) # (M+1)*(L-1)
  ancestors[1,] <- 1
  
  X.ancestors <- matrix(NA, nrow = M+1, u) # we still need the ancestor of the reference to get mu^{(0)}
  X.ancestors[1, ] <- X.proposals[1 , , 1]
  X.ancestors[2:(M+1) , ] <- X.proposals[2:(M+1) , , 1]
  
  time_steps <- c(1, round((obs.times/tau) + 1))
  for(l in 1:L){
    # propagate 
    N.ref_local <- N.ref[ , time_steps[l]:(time_steps[l+1]-1), drop = FALSE] 
    X.ref_local <- X.ref[ , (time_steps[l]+1):time_steps[l+1], drop = FALSE]
    propagation <- propagate_rcpp_SEIR(l, X.ancestors, X.ref_local, N.ref_local, delta, M, time_steps, A, Ns, tau, thetas, coupled.A)
    
    X.proposals[ , , (time_steps[l]+1):time_steps[l+1]] <- propagation$X.proposals_local # (M+1)*u*4
    N.proposals[ , , time_steps[l]:(time_steps[l+1]-1)] <- propagation$N.proposals_local # (M+1)*v*4
    # resample; update X.ancestors
    if(l < L){
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, index_I, time_steps[l+1]], ys[ , l]), numeric(1))
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      w.nor <- wstars/sum(wstars)
      indexes <- sample(1:(M+1), M, prob = w.nor, replace = TRUE) # use multinomial sampling 
      ancestors[2:(M+1),l] <- indexes
      X.ancestors[1,] <- X.proposals[1, , time_steps[l+1]]
      X.ancestors[2:(M+1),] <- X.proposals[indexes, , time_steps[l+1]]
    }else{
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, index_I, ncol(X.ref)], ys[ , L]), numeric(1))
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      w.nor <- wstars/sum(wstars)
    }
  }
  na_num <- sum(apply(X.proposals[ , , time_steps], 1, function(row) any(is.na(row))))
  return(list(weights=w.nor, X.proposals=X.proposals, N.proposals=N.proposals, ancestors=ancestors, na_num=na_num))
}

xPG_SEIR_rcpp <- function(num, delta, M, A, thetas.initial, tau, x0, X.initial, N.initial, loglike, 
                                   ys, obs.times, shape_prior, rate_prior, Ns, coupled.A, thin_interval){
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
  D <- round(u/4)  
  Thetas <- matrix(NA, nrow = length(thetas.initial), ncol = num_thin + 1)
  Thetas[ , 1] <- thetas.initial 
  
  time_steps <- c(1, round((obs.times/tau) + 1)) #indexes of ys in X
  Paths <- array(NA, dim=c(u, L+1, num_thin+1))
  Paths[ , ,1] <- X.initial[ , time_steps]
  # Paths <- array(NA, dim=c(u, n+1, num_thin+1))
  # Paths[ , ,1] <- X.initial
  
  N.path <- N.initial
  path <- X.initial
  
  Z_ref <- matrix(NA, nrow = num_thin, ncol = L) # for alpha_ref
  for (i in 1:num){
    # simulate a path based on lambda
    conPF <- conditional_SEIR_rcpp(A, tau, M, path, N.path, delta, loglike, ys, 
                                            obs.times, x0, thetas, Ns, coupled.A)
    weights <- conPF$weights
    ancestors <- conPF$ancestors
    X.proposals <- conPF$X.proposals
    N.proposals <- conPF$N.proposals
    num_NA[i] <- conPF$na_num
    # k <- systematic_sampling(weights, 1)
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
    for(j in 1:(L-2)){
      index <- ancestors[index,L-j]
      Index[L-j] <- index
      k.2 <- k.1-1
      k.1 <- time_steps[length(time_steps)-j-1]+1
      path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u, ncol = k.2-k.1+1)
      N.path[ , (k.1-1):(k.2-1)] <- matrix(N.proposals[index, , (k.1-1):(k.2-1)], nrow=v, ncol=k.2-k.1+1)
    }
    index <- ancestors[index,1]
    Index[1] <- index
    k.2 <- k.1-1
    k.1 <- 1
    path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u, ncol = k.2-k.1+1)
    N.path[ , 1:(k.2-1)] <- matrix(N.proposals[index, , 1:(k.2-1)], nrow=v, ncol=k.2-1)
    # update all the thetas 
    for(vv in 1:3){
      theta_new <- update_theta_SEIR(vv, N.path, path, shape_prior, rate_prior, tau, coupled.A, Ns)
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

# xPGBS
transition_logprob_SEIR <- function(x.proposals, x.ref, thetas, Ns, M, coupled.A, tau){ 
  # dim(x.proposals) = c(M+1, u); u = D*4
  # length(x.ref) = u
  u <- length(x.ref)
  D <- round(u/4)
  log_trans_prob <- rep(NA, M+1)
  index_S <- 1:D
  index_E <- (D+1):(2*D)
  index_I <- (2*D+1):(3*D)
  index_R <- (3*D+1):(4*D)
  idx_na <- which(apply(x.proposals, 1, function(row) any(is.na(row))))
  
  difference <- matrix(rep(x.ref, each = M + 1), nrow = M + 1) - x.proposals
  
  Reactions <- matrix(NA, nrow = M+1, ncol = round(3*D))
  
  Reactions[ , index_S] <- -difference[ , index_S, drop = FALSE] # R1  # (M+1, D) 
  Reactions[ , index_I] <- difference[ , index_R, drop = FALSE] # R3
  Reactions[ , index_E] <- Reactions[ , index_I, drop = FALSE] +  difference[ , index_I, drop = FALSE] # R2
  # R2 is also R1 - difference[ , index_E, drop = FALSE]
  
  Rates <- matrix(NA, nrow = M+1, ncol = round(3*D))
  for(i in 1:(3*D)){
    Rates[ ,i] <- SEIR.Rates.coupled.vec(i, x.proposals, thetas, coupled.A, Ns)*tau 
  }
  # log Poisson probabilities, vectorised
  log_prob <- dpois(Reactions, Rates, log = TRUE) # (M+1)*(3*D)
  # sum over dimension d
  log_trans_prob <- rowSums(log_prob)
  log_trans_prob[idx_na] <- -Inf
  log_trans_prob[is.na(log_trans_prob)] <- -Inf
  return(log_trans_prob)
}


xPG_SEIR_BS <- function(num, delta, M, A, thetas.initial, tau, x0, X.initial, N.initial, loglike, 
                                ys, obs.times, shape_prior, rate_prior, Ns, coupled.A, thin_interval){
  num_thin <- floor(num/thin_interval)
  
  thetas <- thetas.initial
  num_NA <- rep(NA, num)
  n <- ncol(N.initial) # 200
  
  if (is.null(dim(ys))) {
    ys <- matrix(ys, nrow = 1)
  }
  L <- ncol(ys)
  
  u <- nrow(A)
  v <- ncol(A)
  Thetas <- matrix(NA, nrow = length(thetas.initial), ncol = num_thin + 1)
  Thetas[ , 1] <- thetas.initial 
  
  time_steps <- c(1, round((obs.times/tau) + 1)) #indexes of ys in X
  Paths <- array(NA, dim=c(u, L+1, num_thin+1))
  Paths[ , ,1] <- X.initial[ , time_steps]
  # Paths <- array(NA, dim=c(u, n+1, num_thin+1))
  # Paths[ , ,1] <- X.initial
  
  N.path <- N.initial
  path <- X.initial
  
  Z_ref <- matrix(NA, nrow = num_thin, ncol = L) # for alpha_ref
  D <- round(u/4)
  index_S <- 1:D
  index_E <- (D+1):(2*D)
  index_I <- (2*D+1):(3*D)
  index_R <- (3*D+1):(4*D)
  # set.seed(12345)
  # i <- 1
  for (i in 1:num){
    # simulate a path based on lambda
    conPF <- conditional_SEIR_rcpp(A, tau, M, path, N.path, delta, loglike, ys,
                                            obs.times, x0, thetas, Ns, coupled.A)
    weights <- conPF$weights
    ancestors <- conPF$ancestors
    X.proposals <- conPF$X.proposals
    N.proposals <- conPF$N.proposals
    num_NA[i] <- conPF$na_num
    k <- sample(1:(M+1), 1, prob = weights, replace = TRUE) # multinomial sampling
    
    ## obtain the new path back from k
    # backward sampling
    path <- matrix(NA, nrow = u, ncol = n+1)
    k.2 <- n+1
    k.1 <- time_steps[length(time_steps)-1]+1
    Index <- rep(NA, L)
    index <- k
    Index[L] <- index
    path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u, ncol = k.2-k.1+1)
    for(j in 1:(L-2)){
      l <- L-j  # resample at the l-th obs
      # BS
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, index_I, time_steps[l+1]], ys[ , l]), numeric(1))
      log.tran.prob <- transition_logprob_SEIR(X.proposals[ , , time_steps[l+1]], path[ , time_steps[l+1]+1],thetas, Ns, M, coupled.A, tau) # length = (M+1) 
      
      log.w.tilde <- logW + log.tran.prob
      log.w.tilde.max <- max(log.w.tilde) 
      wstars.tilde <- exp(log.w.tilde - log.w.tilde.max)
      w.nor.tilde <-  wstars.tilde/sum(wstars.tilde)
      
      index <- sample(1:(M+1), 1, prob = w.nor.tilde, replace = TRUE)
      
      Index[l] <- index
      k.2 <- k.1-1
      k.1 <- time_steps[length(time_steps)-j-1]+1
      path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u, ncol = k.2-k.1+1)
    }
    # First observation
    logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, index_I, time_steps[2]], ys[ , 1]), numeric(1))
    log.tran.prob <- transition_logprob_SEIR(X.proposals[ , , time_steps[2]], path[ , time_steps[2]+1],thetas, Ns, M, coupled.A, tau) # length = (M+1) 
    
    log.w.tilde <- logW + log.tran.prob
    log.w.tilde.max <- max(log.w.tilde) 
    wstars.tilde <- exp(log.w.tilde - log.w.tilde.max)
    w.nor.tilde <-  wstars.tilde/sum(wstars.tilde)
    
    index <- sample(1:(M+1), 1, prob = w.nor.tilde, replace = TRUE)
    
    Index[1] <- index
    k.2 <- k.1-1
    k.1 <- 1
    path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u, ncol = k.2-k.1+1)
    
    N.path <- matrix(NA, nrow = v, ncol = n)
    diff <- path[ , 2:(n+1)] - path[ , 1:n] # u*n; u=8
    N.path[index_S, ] <- -diff[index_S, ]
    N.path[index_I, ] <- diff[index_R, ]
    N.path[index_E, ] <- N.path[index_I, ] + diff[index_I, ]
    
    # update all the thetas 
    for(vv in 1:3){
      theta_new <- update_theta_SEIR(vv, N.path, path, shape_prior, rate_prior, tau, coupled.A, Ns)
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


