get.type.reactions <- function(path, matrix.S){
  num_reacts <- nrow(path$xs)
  reacts <- path$xs[2:num_reacts, ]-path$xs[1:(num_reacts-1), ] # 2300*5
  
  reaction_types <- t(matrix.S)
  # Create storage for one-hot encoding
  type.reacts <- matrix(0, nrow = nrow(reacts), ncol = nrow(reaction_types))
  # Match each row of reacts to the reaction type
  for (i in 1:nrow(reacts)) {
    # find which reaction type matches this row
    match_idx <- which(apply(reaction_types, 1, function(r) all(r == reacts[i, ])))
    if (length(match_idx) == 1) {
      type.reacts[i, match_idx] <- 1
    } else {
      warning(sprintf("Row %d did not match uniquely", i))
    }
  }
  return(type.reacts)
}

get.number.reacts.tau <-function(tau.times, path, type.reacts) {
  react.times <- path$ts[2: length(path$ts)] # exclude time zero
  ntau = length(tau.times) # 200
  v=ncol(type.reacts) # 8
  nreacts=matrix(nrow=v, ncol=ntau)
  npatht=nrow(path$xs) 
  ind <- 0
  for (i in 1:ntau) {
    ind.start <- ind + 1
    # cat("ind.start is ", ind.start)
    while ((ind < npatht) && (react.times[ind+1] < tau.times[i])) {
      ind <- ind + 1
    }
    # cat("ind is ", ind)
    nreacts[ , i] <- colSums(type.reacts[ind.start:ind, , drop = FALSE])
  }
  return(nreacts)
}

Autoreg.Rates.vec<-function(r,xs,thetas,G) {
  # rat=-1
  # dim(xs)=(M+1)*(u-1)
  if (r==1) {  ## DNA+P2->DNA.P2
    rat=thetas[1]*(G-xs[ , 4])*xs[ , 3]
  }
  else if (r==2) { ## DNA.P2->DNA+P2
    rat=thetas[2]*xs[ , 4]
  }
  else if (r==3) { ## DNA->DNA+RNA
    rat=thetas[3]*(G-xs[ , 4])
  }
  else if (r==4) { ## RNA->RNA+P
    rat=thetas[4]*xs[ , 1]
  }
  else if (r==5) { ## 2P->P2
    rat=thetas[5]*xs[ , 2]*(xs[ , 2]-1)/2
  }
  else if (r==6) {
    rat=thetas[6]*xs[ , 3] ## P2->2P
  }
  else if (r==7) { ## RNA->0
    rat=thetas[7]*xs[ , 1]
  }
  else if (r==8) { ## P->0
    rat=thetas[8]*xs[ , 2]
  }
  # length(rat)=M+1 
  return (rat)
}

update_proposals <- function(indices, direction, sorted_indices, Mu, tilde_N_proposals, delta) {
  if (direction == "forward") {
    for (i in indices) {
      current <- sorted_indices[i]
      next_i <- sorted_indices[i + 1]
      tilde_N_proposals[next_i] <- tilde_N_proposals[current] + rpois(1, (Mu[next_i] - Mu[current]) * (1 - delta))
    }
  } else {
    for (i in indices) {
      current <- sorted_indices[i]
      prev_i <- sorted_indices[i - 1]
      if (Mu[current] > 0) {
        tilde_N_proposals[prev_i] <- rbinom(1, tilde_N_proposals[current], Mu[prev_i] / Mu[current])
      } else {
        tilde_N_proposals[prev_i] <- 0
      }
    }
  }
  tilde_N_proposals
}

# Not using it anymore
Autoreg_boundary <- function(N.proposals_local, X.ancestors, X.ref_local, k, u, M, G){
  N.trunc <- N.proposals_local[ , , k] # (M+1)*8
  X.potential <- matrix(NA, nrow =M+1, ncol = u-1)
  X.potential[1, ] <- X.ref_local[ ,k]
  X.potential[2:(M+1), ] <- round(X.ancestors[2:(M+1), ] + N.trunc[2:(M+1), ] %*% t(A[1:4,]))
  if(any(N.trunc[2:(M+1), 7] > X.ancestors[2:(M+1), 1])){ # R7>X1
    idx <- which(N.trunc[, 7] > X.ancestors[, 1]) 
    N.trunc[idx, 7] <- X.ancestors[idx, 1]
  }
  if(any(X.potential[2:(M+1), ]<0)||any(X.potential[2:(M+1), 4]>G)){# R1>X5+R2
    gamma1 <- X.ancestors[, 2]/(2*N.trunc[, 5]+N.trunc[, 8]-N.trunc[, 4])
    gamma2 <- X.ancestors[, 3]/(N.trunc[, 1]-N.trunc[, 2]-N.trunc[, 5]+N.trunc[, 6])
    gamma3 <- X.ancestors[, 4]/(N.trunc[, 2]-N.trunc[, 1])
    gamma4 <- (G-X.ancestors[, 4])/(N.trunc[, 1]-N.trunc[, 2])
    gamma_mat <- cbind(gamma1, gamma2, gamma3, gamma4)
    
    idx_row <- which( rowSums(X.potential < 0) > 0 | X.potential[,4] > G )
    idx_min_pos <- apply(gamma_mat[idx_row, ,drop = FALSE], 1, function(x) {
      pos_idx <- which(x > 0)
      if (length(pos_idx) == 0) NA else pos_idx[which.min(x[pos_idx])]
    })
    gamma_min_pos <- gamma_mat[cbind(idx_row, idx_min_pos)]
    N.trunc[idx_row, c(1:6, 8)] <- floor(gamma_min_pos*N.trunc[idx_row, c(1:6, 8)])
  }
  X.proposals_local_k_trunc <- round(X.ancestors[2:(M+1), ] + N.trunc[2:(M+1), ] %*% t(A[1:4,]))
  if(any(X.proposals_local_k_trunc<0)){
    X.proposals_local_k_trunc[X.proposals_local_k_trunc < 0] <- 0
  }
  return(X.proposals_local_k_trunc)
}

propagate <- function(l, X.ancestors, X.ref_local, 
                      N.ref_local, delta, M, time_steps, A, G, tau, thetas){
  u <- nrow(A)
  v <- ncol(A)
  N.proposals_local <- array(NA, dim = c(M+1, v, time_steps[l+1]-time_steps[l]))
  X.proposals_local <- array(NA, dim = c(M+1, u-1, time_steps[l+1]-time_steps[l]))
  
  N.proposals_local[1, , ] <- N.ref_local # 8*4
  X.proposals_local[1, , ] <- X.ref_local # 4*4
  
  for(k in (1:(time_steps[l+1]-time_steps[l]))){
    Mu <- sapply(1:v, function(i) {
      Autoreg.Rates.vec(i, X.ancestors, thetas, G)*tau
    }) # (M+1)*v
    Mu[is.na(Mu)] <- 0 # replace NA with zeros
    sorted_indices <- apply(Mu, 2, order) # (M+1)*v
    
    tilde.N.ref <- rbinom(v, N.ref_local[ , k], 1-delta) # length=8
    tilde.N.proposals <- matrix(NA, nrow = M+1, ncol = v) # (M+1)*v
    
    c <- apply(sorted_indices, 2, function(col) which(col == 1)) # length = v
    tilde.N.proposals[1, ] <- tilde.N.ref
    for(i in 1:v){
      if (c[i] == 1) {
        # Forward pass
        tilde.N.proposals[,i] <- update_proposals(c[i]:M, "forward", sorted_indices[,i], Mu[,i], tilde.N.proposals[,i], delta)
      } else if (c[i] == (M + 1)) {
        # Backward pass
        tilde.N.proposals[,i] <- update_proposals(c[i]:2, "backward", sorted_indices[,i], Mu[,i], tilde.N.proposals[,i], delta)
      } else if (c[i] < (M + 1) && c[i] > 1) {
        # Backward pass for part of the range
        tilde.N.proposals[,i] <- update_proposals(c[i]:2, "backward", sorted_indices[,i], Mu[,i], tilde.N.proposals[,i], delta)
        # Forward pass for part of the range
        tilde.N.proposals[,i] <- update_proposals(c[i]:M, "forward", sorted_indices[,i], Mu[,i], tilde.N.proposals[,i], delta)
      }
    }
    N.proposals_local[2:(M+1), , k] <- tilde.N.proposals[2:(M+1),] + matrix(rpois(M*v, lambda = delta * Mu[2:(M+1),]), nrow = M, ncol = v)
    # if(any(is.na(N.proposals_local[2:(M+1), , k]))){
    #   print("Mu is ")
    #   print(Mu)
    #   print("tilde.N.proposals is")
    #   print(tilde.N.proposals)
    #   print("X.ancestors is ")
    #   print(X.ancestors)
    # }
    N.trunc <- N.proposals_local[ , , k] # (M+1)*8
    
    X.potential <- matrix(NA, nrow =M+1, ncol = u) # (M+1)*u
    X.potential[1, 1:(u-1)] <- X.ref_local[ ,k]
    X.potential[2:(M+1), 1:(u-1)] <- round(X.ancestors[2:(M+1), ] + N.trunc[2:(M+1), ] %*% t(A[1:4,])) # may has NA
    X.potential[ , u] <- G - X.potential[ , u-1]
    if(any(X.potential[2:(M+1), ]<0)){
      neg_rows <- which(apply(X.potential, 1, function(row) any(row < 0, na.rm = TRUE)))
      pos_rows <- setdiff(1:(M + 1), neg_rows)
      X.proposals_local[neg_rows, , k] <- NA
      N.proposals_local[neg_rows, , k] <- NA
    }else{
      pos_rows <- 2:(M+1)
    }
    X.proposals_local[pos_rows, , k] <- X.potential[pos_rows, 1:(u-1)]
    # update X.ancestors
    X.ancestors[ , ] <- X.proposals_local[ , , k]
    # print(X.ancestors)
    # print(k)
  }
  return(list(X.proposals_local=X.proposals_local, N.proposals_local=N.proposals_local))
}

loglike_autoreg <- function(xt, y){
  # length(xt) = 3 
  # length(y) = 3  
  diff <- sum(abs(xt - y))
  tol <- 1e-8
  likelihood <- as.integer(diff < tol)
  return(log(likelihood))
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

conditional_Autoreg <- function(A, tau, M, X.ref, N.ref, delta, loglike, ys, 
                                obs.times, x0, thetas){
  n <- ncol(N.ref) # 200
  L <- ncol(ys)
  u <- nrow(A)
  v <- ncol(A)
  G <- x0[4] + x0[5]
  
  N.proposals <- array(NA, dim = c(M+1, v, n))
  N.proposals[1, , ] <- N.ref
  
  X.proposals <- array(NA, dim = c(M+1, u-1, n+1)) # (M+1)*4*201-dim
  X.proposals[1, , ] <- X.ref # 4*201
  X.proposals[2:(M+1), , 1] <- matrix(x0[1:(u-1)], nrow = M, ncol = u-1, byrow = TRUE) # M*1-dim
  
  ancestors <- matrix(NA, nrow = M+1, ncol = L-1) # (M+1)*(L-1)
  ancestors[1,] <- 1
  ancestors[2:(M+1), 1] <- 2:(M+1)
  
  X.ancestors <- matrix(NA, nrow = M+1, u-1) # we still need the ancestor of the reference to get mu^{(0)}
  X.ancestors[ , ] <- X.proposals[ , , 1]
  
  obs.times <- ts
  time_steps <- round((obs.times/tau) + 1)
  for(l in 2:L){
    # propagate 
    N.ref_local <- N.ref[ , time_steps[l-1]:(time_steps[l]-1)] # 8*4
    X.ref_local <- X.ref[ , (time_steps[l-1]+1):time_steps[l]] # 4*4
    propagation <- propagate(l-1, X.ancestors, X.ref_local, N.ref_local, delta, M, time_steps, A, G, tau, thetas)
    # X.ancestors[ , ] <- X.proposals_local[ , , k]
    X.proposals[ , , (time_steps[l-1]+1):time_steps[l]] <- propagation$X.proposals_local # (M+1)*u*4
    N.proposals[ , , time_steps[l-1]:(time_steps[l]-1)] <- propagation$N.proposals_local # (M+1)*v*4
    # resample; update X.ancestors
    if(l < L){
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, 1:3, time_steps[l]], ys[ , l]), numeric(1))
      # logW <- loglike(X.proposals[,1,time_steps[l+1]], ys[l], N, prob.obs)
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      y.loglikely <- log(mean(wstars)) + logW.max # maximum marginal loglikelhood
      w.nor <- wstars/sum(wstars)
      indexes <- systematic_sampling(w.nor, M) # sample M from 1:(M+1)
      ancestors[2:(M+1),l] <- indexes
      X.ancestors[1,] <- X.proposals[1, , time_steps[l]]
      X.ancestors[2:(M+1),] <- X.proposals[indexes, , time_steps[l]]
    }else{
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, 1:3, ncol(X.ref)], ys[1:3, L]), numeric(1))
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      y.loglikely <- log(mean(wstars)) + logW.max # maximum marginal loglikelhood
      w.nor <- wstars/sum(wstars)
    }
  }
  na_num <- sum(apply(X.proposals[ , , n+1], 1, function(row) any(is.na(row))))
  return(list(weights=w.nor, X.proposals=X.proposals, N.proposals=N.proposals, ancestors=ancestors, na_num=na_num))
}


xPG <- function(num, delta, M, A, thetas.initial, tau, x0, X.initial, N.initial, loglike, 
                ys, obs.times, shape_prior_theta1, rate_prior_theta1, shape_prior_theta3, rate_prior_theta3){
  G <- x0[4] + x0[5]
  
  Theta1 <- rep(NA, num+1)
  Theta1[1] <- thetas.initial[1] 
  
  Theta3 <- rep(NA, num+1)
  Theta3[1] <- thetas.initial[3]
  
  thetas <- thetas.initial
  num_NA <- rep(NA, num)
  n <- ncol(N.ref) # 200
  L <- ncol(ys)
  u <- nrow(A)
  v <- ncol(A)
  Paths <- array(NA, dim=c(u-1, n+1, num+1))
  Paths[ , ,1] <- X.initial
  # Paths <- array(NA, dim=c(1, L, num))
  N.path <- N.initial
  path <- X.initial
  time_steps <- round((obs.times/tau) + 1) #indexes of ys in X
  Z_val <- matrix(NA, nrow = num, ncol = L) # for alpha_val of X4
  Z_ref <- matrix(NA, nrow = num, ncol = L) # for alpha_ref
  set.seed(12345)
  i <- 1
  for (i in 1:num){
    # simulate a path based on lambda
    path.ref <- path # store the reference path for computing acc rate
    conPF <- conditional_Autoreg(A, tau, M, path, N.path, delta, loglike=loglike_autoreg, ys,
                                 obs.times, x0, thetas)
    weights <- conPF$weights
    ancestors <- conPF$ancestors
    X.proposals <- conPF$X.proposals
    N.proposals <- conPF$N.proposals
    num_NA[i] <- conPF$na_num
    k <- systematic_sampling(weights, 1)
    # k <- sample(1:(M+1), 1, prob = weights)
    ## obtain the new path back from k
    path <- matrix(NA, nrow = u-1, ncol = n+1)
    N.path <- matrix(NA, nrow = v, ncol = n)
    k.2 <- n+1
    k.1 <- time_steps[length(time_steps)-1]+1
    Index <- rep(NA, L)
    index <- k
    Index[L] <- index
    path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u-1, ncol = k.2-k.1+1)
    N.path[ , (k.1-1):(k.2-1)] <- matrix(N.proposals[index, , (k.1-1):(k.2-1)], nrow = v, ncol = k.2-k.1+1)
    # N.1.path[(k.1-1):(k.2-1)] <- con_SIS$N.1.proposals[index, (k.1-1):(k.2-1)]
    # N.2.path[(k.1-1):(k.2-1)] <- con_SIS$N.2.proposals[index, (k.1-1):(k.2-1)]
    for(j in 1:(L-2)){
      index <- ancestors[index,L-j]
      Index[L-j] <- index
      k.2 <- k.1-1
      k.1 <- time_steps[length(time_steps)-j-1]+1
      path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u-1, ncol = k.2-k.1+1)
      N.path[ , (k.1-1):(k.2-1)] <- matrix(N.proposals[index, , (k.1-1):(k.2-1)], nrow=v, ncol=k.2-k.1+1)
    }
    # index <- ancestors[index,1]
    Index[1] <- index
    path[ , 1] <- matrix(X.proposals[index, ,1], nrow = u-1, ncol = 1)
    Paths[, ,i+1] <- path
    # Paths[, ,i] <- path[1, time_steps]
    z_val <- colSums(ifelse(path[ , time_steps] != path.ref[ , time_steps], 1, 0)) # for alpha_val
    Z_val[i, ] <- z_val
    z_ref <- ifelse(Index != 1, 1, 0) # for computing alpha_ref
    Z_ref[i, ] <- z_ref
    # update mu given N.1.path, N.2.path and path
    # suppose the prior distribution of mu is Gamma(2, 4)
    theta1 <- rgamma(1, shape = shape_prior_theta1 + sum(N.path[1, ]), 
                     rate = rate_prior_theta1 + tau*sum((G-path[4, 1:n])*path[3, 1:n]))
    
    theta3 <- rgamma(1, shape = shape_prior_theta3 + sum(N.path[3, ]), 
                     rate = rate_prior_theta3 + tau*sum(G-path[4, 1:n]))
    thetas[1] <- theta1
    thetas[3] <- theta3
    Theta1[i+1] <- theta1
    Theta3[i+1] <- theta3
  }
  return(list(Paths=Paths, Z_val=Z_val, Z_ref=Z_ref, Theta1=Theta1, Theta3=Theta3, num_NA=num_NA))
}

plot.path.lanya <- function(tau.times, obs.times, X.proposals, ys, i, M, textsc, ...){
  y_vals <- range(X.proposals[, i, ])
  if (i <= nrow(ys)) {
    y_vals <- range(y_vals, ys[i, ])
  }
  plot(c(0,tau.times), X.proposals[1, i, ], type = "l", ylim = y_vals,
       cex.lab=textsc,cex.main=textsc,cex.axis=textsc,col="red", ...)
  if (i <= nrow(ys)) {
    points(obs.times, ys[i, ])
  }
  for(m in 2:(M+1)){
    lines(c(0,tau.times), X.proposals[m, i, ], type = "l", lty = 2, col="grey35")
  }
}

### Rcpp version of conditional particle filter and xPG ####

update_theta <- function(vv, N.path, path, shape_prior, rate_prior, tau, G){
  n <- ncol(N.path)
  if(vv == 1){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[1, ]), 
                        rate = rate_prior + tau*sum((G-path[4, 1:n])*path[3, 1:n]))
  }
  if(vv == 2){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[2, ]), 
                        rate = rate_prior + tau*sum(path[4, 1:n]))
  }
  if(vv == 3){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[3, ]), 
                        rate = rate_prior + tau*sum((G-path[4, 1:n])))
  }
  if(vv == 4){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[4, ]), 
                        rate = rate_prior + tau*sum(path[1, 1:n]))
  }
  if(vv == 5){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[5, ]), 
                        rate = rate_prior + tau*sum(path[2, 1:n]*(path[2, 1:n]-1)/2))
  }
  if(vv == 6){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[6, ]), 
                        rate = rate_prior + tau*sum(path[3, 1:n]))
  }
  if(vv == 7){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[7, ]), 
                        rate = rate_prior + tau*sum(path[1, 1:n]))
  }
  if(vv == 8){
    theta_new <- rgamma(1, shape = shape_prior + sum(N.path[8, ]), 
                        rate = rate_prior + tau*sum(path[2, 1:n]))
  }
  return(theta_new)
}

conditional_Autoreg_rcpp <- function(A, tau, M, X.ref, N.ref, delta, loglike, ys, 
                                     obs.times, x0, thetas){
  n <- ncol(N.ref) # 200
  L <- ncol(ys)
  u <- nrow(A)
  v <- ncol(A)
  G <- x0[4] + x0[5] # suppose we do not know x0[4] and x0[5] but we know G
  
  N.proposals <- array(NA, dim = c(M+1, v, n))
  N.proposals[1, , ] <- N.ref
  
  X.proposals <- array(NA, dim = c(M+1, u-1, n+1)) # (M+1)*4*201-dim
  X.proposals[1, , ] <- X.ref # 4*201
  # X.proposals[2:(M+1), , 1] <- matrix(x0[1:(u-1)], nrow = M, ncol = u-1, byrow = TRUE) # M*1-dim
  X.proposals[2:(M+1), 1:(u-2), 1] <- matrix(x0[1:(u-2)], nrow = M, ncol = u-2, byrow = TRUE) # x1, x2 and x3 are known
  x4.initial <- sample(0:G, M, replace = TRUE) # uniformly sample the value of x4 at time zero
  X.proposals[2:(M+1), 4, 1] <- x4.initial
  
  ancestors <- matrix(NA, nrow = M+1, ncol = L-1) # (M+1)*(L-1)
  ancestors[1,] <- 1
  ancestors[2:(M+1), 1] <- 2:(M+1)
  
  X.ancestors <- matrix(NA, nrow = M+1, u-1) # we still need the ancestor of the reference to get mu^{(0)}
  X.ancestors[ , ] <- X.proposals[ , , 1]
  
  obs.times <- ts
  time_steps <- round((obs.times/tau) + 1)
  for(l in 2:L){
    # propagate 
    N.ref_local <- N.ref[ , time_steps[l-1]:(time_steps[l]-1)] # 8*4
    X.ref_local <- X.ref[ , (time_steps[l-1]+1):time_steps[l]] # 4*4
    propagation <- propagate_rcpp(l-1, X.ancestors, X.ref_local, N.ref_local, delta, M, time_steps, A, G, tau, thetas, rate_func = Autoreg_Rates_vec_cpp)
    # X.ancestors[ , ] <- X.proposals_local[ , , k]
    X.proposals[ , , (time_steps[l-1]+1):time_steps[l]] <- propagation$X.proposals_local # (M+1)*u*4
    N.proposals[ , , time_steps[l-1]:(time_steps[l]-1)] <- propagation$N.proposals_local # (M+1)*v*4
    # resample; update X.ancestors
    if(l < L){
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, 1:2, time_steps[l]], ys[ , l]), numeric(1))
      # logW <- loglike(X.proposals[,1,time_steps[l+1]], ys[l], N, prob.obs)
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      y.loglikely <- log(mean(wstars)) + logW.max # maximum marginal loglikelhood
      w.nor <- wstars/sum(wstars)
      # indexes <- systematic_sampling(w.nor, M) # sample M from 1:(M+1)
      indexes <- sample(1:(M+1), M, prob = w.nor, replace = TRUE) # use multinomial sampling 
      ancestors[2:(M+1),l] <- indexes
      X.ancestors[1,] <- X.proposals[1, , time_steps[l]]
      X.ancestors[2:(M+1),] <- X.proposals[indexes, , time_steps[l]]
    }else{
      logW <- vapply(1:(M + 1), function(i) loglike(X.proposals[i, 1:2, ncol(X.ref)], ys[ , L]), numeric(1))
      logW[is.na(logW)] <- -Inf
      logW.max <- max(logW)
      wstars <- exp(logW - logW.max) # length M
      y.loglikely <- log(mean(wstars)) + logW.max # maximum marginal loglikelhood
      w.nor <- wstars/sum(wstars)
    }
  }
  na_num <- sum(apply(X.proposals[ , , n+1], 1, function(row) any(is.na(row))))
  return(list(weights=w.nor, X.proposals=X.proposals, N.proposals=N.proposals, ancestors=ancestors, na_num=na_num))
}


xPG_rcpp <- function(num, delta, M, A, thetas.initial, tau, x0, X.initial, N.initial, loglike, 
                     ys, obs.times, shape_prior, rate_prior){
  G <- x0[4] + x0[5]
  
  thetas <- thetas.initial
  num_NA <- rep(NA, num)
  n <- ncol(N.ref) # 200
  L <- ncol(ys)
  u <- nrow(A)
  v <- ncol(A)
  Thetas <- matrix(NA, nrow = v, ncol = num + 1)
  Thetas[ , 1] <- thetas.initial 
  Paths <- array(NA, dim=c(u-1, n+1, num+1))
  Paths[ , ,1] <- X.initial
  # Paths <- array(NA, dim=c(1, L, num))
  N.path <- N.initial
  path <- X.initial
  time_steps <- round((obs.times/tau) + 1) #indexes of ys in X
  Z_val <- matrix(NA, nrow = num, ncol = L) # for alpha_val of X4
  Z_ref <- matrix(NA, nrow = num, ncol = L) # for alpha_ref
  set.seed(12345)
  i <- 1
  for (i in 1:num){
    # simulate a path based on lambda
    path.ref <- path # store the reference path for computing acc rate
    conPF <- conditional_Autoreg_rcpp(A, tau, M, path, N.path, delta, loglike=loglike_autoreg, ys,
                                      obs.times, x0, thetas)
    weights <- conPF$weights
    ancestors <- conPF$ancestors
    X.proposals <- conPF$X.proposals
    N.proposals <- conPF$N.proposals
    num_NA[i] <- conPF$na_num
    # k <- systematic_sampling(weights, 1)
    k <- sample(1:(M+1), 1, prob = weights, replace = TRUE) # multinomial sampling
    ## obtain the new path back from k
    path <- matrix(NA, nrow = u-1, ncol = n+1)
    N.path <- matrix(NA, nrow = v, ncol = n)
    k.2 <- n+1
    k.1 <- time_steps[length(time_steps)-1]+1
    Index <- rep(NA, L)
    index <- k
    Index[L] <- index
    path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u-1, ncol = k.2-k.1+1)
    N.path[ , (k.1-1):(k.2-1)] <- matrix(N.proposals[index, , (k.1-1):(k.2-1)], nrow = v, ncol = k.2-k.1+1)
    # N.1.path[(k.1-1):(k.2-1)] <- con_SIS$N.1.proposals[index, (k.1-1):(k.2-1)]
    # N.2.path[(k.1-1):(k.2-1)] <- con_SIS$N.2.proposals[index, (k.1-1):(k.2-1)]
    for(j in 1:(L-2)){
      index <- ancestors[index,L-j]
      Index[L-j] <- index
      k.2 <- k.1-1
      k.1 <- time_steps[length(time_steps)-j-1]+1
      path[ , k.1:k.2] <- matrix(X.proposals[index, ,k.1:k.2], nrow = u-1, ncol = k.2-k.1+1)
      N.path[ , (k.1-1):(k.2-1)] <- matrix(N.proposals[index, , (k.1-1):(k.2-1)], nrow=v, ncol=k.2-k.1+1)
    }
    # index <- ancestors[index,1]
    Index[1] <- index
    path[ , 1] <- matrix(X.proposals[index, ,1], nrow = u-1, ncol = 1)
    Paths[, ,i+1] <- path
    # Paths[, ,i] <- path[1, time_steps]
    z_val <- colSums(ifelse(path[ , time_steps] != path.ref[ , time_steps], 1, 0)) # for alpha_val
    Z_val[i, ] <- z_val
    z_ref <- ifelse(Index != 1, 1, 0) # for computing alpha_ref
    Z_ref[i, ] <- z_ref
    # update all the thetas 
    for(vv in 1:v){
      theta_new <- update_theta(vv, N.path, path, shape_prior[vv], rate_prior[vv], tau, G)
      Thetas[vv, i+1] <- theta_new
    }
    thetas <-Thetas[ , i+1]
  }
  return(list(Paths=Paths, Z_val=Z_val, Z_ref=Z_ref, Thetas=Thetas, num_NA=num_NA))
}


