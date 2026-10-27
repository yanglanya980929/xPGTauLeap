tau_leap_SIR <- function(A, lambda, mu, initial_values, n, tau, N){
  u <- nrow(A)
  X <- matrix(NA, nrow = u, ncol = n+1)
  X[,1] <- initial_values
  k = 1
  R.1 <- rep(NA, n)
  R.2 <- rep(NA, n)
  for(k in 1:n){
    mu.1 <- X[1,k]*X[2,k]*lambda*tau/N
    mu.2 <- X[2,k]*mu*tau
    r.1 <- rpois(1, lambda = mu.1) # number of reaction1
    r.2 <- rpois(1, lambda = mu.2) # number of reaction2
    if(r.1 > X[1,k]){
      r.1 <- X[1,k]
    }
    if(r.2 > X[2,k]){
      r.1 <- X[2,k]
    }
    R.1[k] <- r.1
    R.2[k] <- r.2
    X[,k+1] <- X[,k]+A[,1]*r.1+A[,2]*r.2
  }
  times <- (1:n)*tau
  times <- c(0, times)
  # X <- X[,1:n]
  return(list(X = X, times = times, R.1 = R.1, R.2 = R.2))
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

tau_leap_SIR_multi <- function(A, lambda, mu, initial_values, n, tau, N, M){
  u <- nrow(A)
  X <- array(NA, dim = c(M, u, n+1))
  # X <- matrix(NA, nrow = u, ncol = n+1)
  X[ , , 1] <- t(replicate(M, initial_values))
  k = 1
  R.1 <- matrix(NA, nrow = M, ncol = n)
  R.2 <- matrix(NA, nrow = M, ncol = n)
  for(k in 1:n){
    mu.1 <- X[ , 1, k]*X[ , 2, k]*lambda*tau/N
    mu.2 <- X[ , 2, k]*mu*tau
    r.1 <- rpois(M, lambda = mu.1) # number of reaction1
    r.2 <- rpois(M, lambda = mu.2) # number of reaction2
    if(any(r.1 > X[ , 1, k])){
      ind <- which(r.1 > X[ , 1, k])
      r.1[ind] <- X[ind, 1, k]
    }
    if(any(r.2 > X[ , 2, k])){
      ind <- which(r.2 > X[ , 2, k])
      r.2[ind] <- X[ind, 2, k]
    }
    R.1[ , k] <- r.1
    R.2[ , k] <- r.2
    X[ , , k+1] <- X[ , , k] + r.1%*%t(A[,1]) + r.2%*%t(A[,2])
  }
  times <- (1:n)*tau
  times <- c(0, times)
  return(list(X = X, times = times, R.1 = R.1, R.2 = R.2))
}

loglike.SIR.bin <- function(xs, y, p){
  # xs: length D
  # y: length D
  if(anyNA(xs)||anyNA(y)){
    return(-Inf)
  } else{
    l <- dbinom(y, size = xs, prob = p, log = TRUE)
    return(sum(l))
  }
}

SIR.Rates.coupled <-function(r,xs,thetas, coupled.A, Ns) {
  # xs: length D*u
  # couple.A: D*D
  # N: length D
  # rat: length D
  D <- round(length(xs)/3)
  index_S <- 1:D
  index_I <- (D+1):(2*D)
  R1 = thetas[1]*xs[index_S]*(coupled.A %*% xs[index_I]) # (M+1)*D
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

tau_leap_couple_SIR <- function(A, thetas, xs.initial, n, tau, Ns, coupled.A){
  # initial_states: D*u
  u <- nrow(A)
  v <- ncol(A) # number of the types of reactions
  D <- round(u/3)
  index_S <- 1:D
  index_I <- (D+1):(2*D)
  X <- matrix(NA, nrow = u, ncol = n+1)
  X[ , 1] <- xs.initial
  Rs <- matrix(NA, nrow = v, ncol = n)
  for(k in 1:n){
    Mu <- sapply(1:v, function(i) SIR.Rates.coupled(i, X[ , k], thetas, coupled.A, Ns))*tau # D*v
    rs <- rpois(v, lambda = as.vector(Mu))
    Rs[ , k] <- rs
    if(any(rs[index_S] > X[index_S, k])){
      rs[index_S] <- X[index_S, k]
    }
    if(any(rs[index_I] > X[index_I, k] + rs[index_S])){
      rs[index_I] <- X[index_I, k] + rs[index_S]
    }
    X[ , k+1] <- X[ , k]+rs[index_S]%*%t(A[,index_S])+rs[index_I]%*%t(A[,index_I])
  }
  times <- (1:n)*tau
  times <- c(0, times)
  return(list(X = X, times = times, Rs = Rs))
}

# multi/single xPG;
conditional_SimSIR <- function(A, tau, M, X.ref, N.ref, delta, loglike, ys, 
                                         obs.times, x0, thetas, Ns, coupled.A){
  n <- ncol(N.ref)
  if (is.null(dim(ys))) {
    ys <- matrix(ys, nrow = 1)
  }
  L <- ncol(ys)
  u <- nrow(A)
  v <- ncol(A) # 2*D
  
  stopifnot(u %% 3 == 0)
  D <- round(u/3)
  
  N.proposals <- array(NA, dim = c(M+1, v, n))
  N.proposals[1, , ] <- N.ref
  
  X.proposals <- array(NA, dim = c(M+1, u, n+1)) # (M+1)*3*201-dim
  X.proposals[1, , ] <- X.ref # 3*201
  X.proposals[2:(M+1), , 1] <- matrix(x0, nrow = M, ncol = u, byrow = TRUE) # M*1-dim
  
  index_S <- 1:D
  index_I <- (D+1):(2*D)
  index_R <- (2*D+1):(3*D)
  
  ancestors <- matrix(NA, nrow = M+1, ncol = L-1) # (M+1)*(L-1)
  ancestors[1,] <- 1

  X.ancestors <- matrix(NA, nrow = M+1, u) # we still need the ancestor of the reference to get mu^{(0)}
  X.ancestors[1, ] <- X.proposals[1 , , 1]
  X.ancestors[2:(M+1) , ] <- X.proposals[2:(M+1) , , 1]
  
  time_steps <- c(1, round((obs.times/tau) + 1)) # obs.times does not include 0
  for(l in 1:L){
    # propagate 
    N.ref_local <- N.ref[ , time_steps[l]:(time_steps[l+1]-1), drop = FALSE] 
    X.ref_local <- X.ref[ , (time_steps[l]+1):time_steps[l+1], drop = FALSE] 
    propagation <- propagate_rcpp_SIR(l, X.ancestors, X.ref_local, N.ref_local, delta, M, time_steps, A, Ns, tau, thetas, coupled.A)
    
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
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, index_I, time_steps[L+1]], ys[ , L]), numeric(1))
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      w.nor <- wstars/sum(wstars)
    }
  }
  na_num <- sum(apply(X.proposals[ , , time_steps], 1, function(row) any(is.na(row))))
  return(list(weights=w.nor, X.proposals=X.proposals, N.proposals=N.proposals, ancestors=ancestors, na_num=na_num))
}

update_theta_SIR <- function(vv, N.path, path, shape_prior, rate_prior, tau, coupled.A, Ns){
  n <- ncol(N.path)
  D <- round(nrow(path) / 3)
  
  index_S <- 1:D
  index_I <- (D + 1):(2 * D)
  
  # states at the beginning of each tau interval
  S <- path[index_S, 1:n, drop = FALSE]
  I <- path[index_I, 1:n, drop = FALSE]
  
  # coupled infection pressure
  # dimension: D x n
  infection_pressure <- coupled.A %*% I
  
  # beta multiplier:
  # sum_{k,d} S_{d,k} *
  #              sum_j A_{d,j} I_{j,k} / N_d
  rate_beta <- sum(rowSums(S * infection_pressure) / Ns)
  
  # gamma multiplier:
  # sum_{k,d} I_{d,k}
  rate_gamma <- sum(I)
  
  # Reaction ordering:
  # 1:D             infection
  # (D+1):(2D)      recovery
  count_beta <-  sum(N.path[index_S, , drop = FALSE])
  count_gamma <- sum(N.path[index_I, , drop = FALSE])
  
  if (vv == 1) {
    theta_new <- rgamma(1, shape = shape_prior[1] + count_beta, 
                        rate  = rate_prior[1] + tau * rate_beta)
  } else if (vv == 2) {
    theta_new <- rgamma(1, shape = shape_prior[2] + count_gamma,
                        rate  = rate_prior[2] + tau * rate_gamma)
  }
  return(theta_new)
}

xPG_SIR <- function(num, delta, M, A, thetas.initial, tau, x0, X.initial, N.initial, loglike, 
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
  D <- floor(u/3)  
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
    conPF <- conditional_SimSIR(A, tau, M, path, N.path, delta, loglike, ys, 
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
    for(vv in 1:2){
      theta_new <- update_theta_SIR(vv, N.path, path, shape_prior, rate_prior, tau, coupled.A, Ns)
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


# SIR: BS + xPG
SIR.Rates.coupled.vec <- function(r, xs, thetas, coupled.A, Ns) {
  
  # xs: (M+1) x (3D)
  # state ordering: S_1,...,S_D, I_1,...,I_D, R_1,...,R_D
  # coupled.A: D x D
  # Ns: length D
  # r = 1,...,2D
  #
  # reactions:
  # 1:D         infection
  # (D+1):(2D)  recovery
  stopifnot(ncol(xs) %% 3 == 0)
  D <- floor(ncol(xs) / 3)
  stopifnot(
    length(thetas) == 2,
    length(Ns) == D,
    all(dim(coupled.A) == c(D, D)),
    r >= 1,
    r <= 2 * D
  )
  index_S <- 1:D
  index_I <- (D + 1):(2 * D)
  S <- xs[, index_S, drop = FALSE]
  I <- xs[, index_I, drop = FALSE]
  
  # (M+1) x D
  infection_pressure <- I %*% t(coupled.A)
  
  # infection rates
  R1 <- thetas[1] * S * infection_pressure
  # recovery rates
  R2 <- thetas[2] * I
  if (r <= D) {
    d <- r
    rat <- R1[, d] / Ns[d]
  } else {
    d <- r - D
    rat <- R2[, d]
  }
  return(rat)
}

transition_logprob_SIR <- function(x.proposals, x.ref, thetas, Ns, M, coupled.A, tau){ 
  # dim(x.proposals) = c(M+1, u); u = D*3
  # length(x.ref) = u
  u <- length(x.ref)
  D <- floor(u/3)
  log_trans_prob <- rep(NA, M+1)
  index_S <- 1:D
  index_I <- (D+1):(2*D)
  index_R <- (2*D+1):(3*D)
  idx_na <- which(apply(x.proposals, 1, function(row) any(is.na(row))))
  
  difference <- matrix(rep(x.ref, each = M + 1), nrow = M + 1) - x.proposals
  
  Reactions <- matrix(NA, nrow = M+1, ncol = round(2*D))
  
  Reactions[ , index_S] <- -difference[ , index_S, drop = FALSE] # R1  # (M+1, D) 
  Reactions[ , index_I] <- difference[ , index_R, drop = FALSE] # R2
  
  Rates <- matrix(NA, nrow = M+1, ncol = round(2*D))
  for(i in 1:(2*D)){
    Rates[ ,i] <- SIR.Rates.coupled.vec(i, x.proposals, thetas, coupled.A, Ns)*tau 
  }
  # log Poisson probabilities, vectorised
  log_prob <- dpois(Reactions, Rates, log = TRUE) # (M+1)*(2*D)
  # sum over dimension d
  log_trans_prob <- rowSums(log_prob)
  log_trans_prob[idx_na] <- -Inf
  log_trans_prob[is.na(log_trans_prob)] <- -Inf
  return(log_trans_prob)
}

xPG_SIR_BS <- function(num, delta, M, A, thetas.initial, tau, x0, X.initial, N.initial, loglike, 
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
  D <- round(u/3)  
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
  index_S <- 1:D
  index_I <- (D+1):(2*D)
  index_R <- (2*D+1):(3*D)
  for (i in 1:num){
    # simulate a path based on lambda
    conPF <- conditional_SimSIR(A, tau, M, path, N.path, delta, loglike, ys, 
                                obs.times, x0, thetas, Ns, coupled.A)
    weights <- conPF$weights
    ancestors <- conPF$ancestors
    X.proposals <- conPF$X.proposals
    N.proposals <- conPF$N.proposals
    num_NA[i] <- conPF$na_num
    # k <- systematic_sampling(weights, 1)
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
      log.tran.prob <- transition_logprob_SIR(X.proposals[ , , time_steps[l+1]], path[ , time_steps[l+1]+1],thetas, Ns, M, coupled.A, tau) # length = (M+1) 
      
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
    log.tran.prob <- transition_logprob_SIR(X.proposals[ , , time_steps[2]], path[ , time_steps[2]+1], thetas, Ns, M, coupled.A, tau) # length = (M+1) 
    
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
    
    # update all the thetas 
    for(vv in 1:2){
      theta_new <- update_theta_SIR(vv, N.path, path, shape_prior, rate_prior, tau, coupled.A, Ns)
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


