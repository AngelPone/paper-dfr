library(Matrix)
library(dplyr)
library(quadprog)


construct_Q <- function(r, q, pi_hat) {
  # r: cardinality of reconciled probability
  # q: cardinality of probability
  # pi_hat: q*T matrix
  Q <- pi_hat %*% t(pi_hat)
  bdiag(replicate(r, Q, simplify = FALSE))
}


cons_realDummy <- function(series, coherent_domain) {
  #' transform hierarchical time series into dummy matrix used in Brier Score.
  #' @param series the observed hierarchical time series m * T
  #' @param coherent_domain the domain of 
  r = dim(coherent_domain)[1]
  m = dim(series)[2]
  dummy_mat <- matrix(0, r, dim(series)[1])
  for (i in 1:dim(series)[1]) {
    dummy_mat[which(apply(coherent_domain[, 1:m], 1, function(x) {
      all(x == series[i, ])
    })), i] = 1
  }
  dummy_mat
}
cal_distanceMatrix <- function(incoherent_domain, coherent_domain) {
  r = dim(coherent_domain)[1]
  q = dim(incoherent_domain)[1]
  distance = matrix(NA, nrow = r, ncol = q)
  for (j in 1:r) {
    for (k in 1:q) {
      distance[j, k] = sum(abs(coherent_domain[j, ] - incoherent_domain[k, ]))
    }
  }
  distance
}

#QR with solve.qr

opt_fun <- function(pi_hat, real, distance, lambda) {
  # pi_hat: q * T
  # real: r * T
  # distance: r * q
  r = dim(real)[1]
  q = dim(pi_hat)[1]
  time_window = dim(pi_hat)[2]
  Dmat <- 2 * construct_Q(r, q, pi_hat) / time_window
  Z <- matrix(real, nrow = 1)
  D <- NULL
  for (i in 1:time_window) {
    D <-
      rbind(D, bdiag(replicate(r, matrix(pi_hat[, i], nrow = 1), simplify = FALSE)))
  }
  dvec <-
    -t(lambda * matrix(t(distance), nrow = 1) - 2 * Z %*% D / time_window)
  # check the result in the optimization function equals to the result of matrix form.
  # total_opt = -t(dvec) %*% vec_A + 1/2 * t(vec_A) %*% Dmat%*% vec_A + t(vec_Z) %*% vec_Z/time_window
  A1 <- diag(q)
  for (i in 1:(r - 1)) {
    A1 <- cbind(A1, diag(q))
  }
  b1 <- rep(1, q)
  
  A2 <- diag(r * q)
  b2 <- rep(0, r * q)
  
  A3 <- -diag(r * q)
  b3 <- rep(-1, r * q)
  Amat <- t(rbind(A1, A2, A3))
  bvec <- c(b1, b2, b3)
  
  solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = q)
  
  t(matrix(solution$solution, 12, 4))
}


simulate_series <- function(param_a,
                            n,
                            intercept = 0.0,
                            innovation_mu = c(0, 0),
                            innovation_Sigma = matrix(c(0.1, 0.05, 0.05, 0.1), 2, 2)) {
  # simulate time series using State Space model
  error_term <- MASS::mvrnorm(n + 1, innovation_mu, innovation_Sigma)
  # latent state
  for (i in (length(param_a) + 1):dim(error_term)[1]) {
    error_term[i, ] = intercept +  + error_term[i, ]
    for (j in 1:length(param_a)) {
       error_term[i, ] = error_term[i, ] + param_a[j] * error_term[i - j, ]
    }
    
  }
  # pi
  pi = 1 / (1 + 1 / exp(error_term))
  ifelse(pi > 0.5, 1, 0)[2:(n + 1),]
}

glarma_base_forecast <- function(series, total = FALSE) {
  if (total) {
    series <- cbind(series, 2 - series)
  } else{
    series <- cbind(series, 1 - series)
  }
  model <-
    glarma::glarma(
      series,
      X = matrix(1, nrow = dim(series)[1]),
      phiLags = c(1),
      phiInit = c(0.7),
      type = 'Bin'
    )
  fore <- glarma::forecast(model, h = 1, newdata = matrix(1))
  
  if (total) {
    pi <- 1 / (1 + 1 / exp(fore$W))
    return(c((1 - pi) ^ 2, 2 * pi * (1 - pi), pi ^ 2))
  } else{
    return(c(1 - fore$mu, fore$mu))
  }
}


brier_score <- function(probf, real) {
  time_window <- dim(probf)[2]
  total = 0
  for (i in 1:time_window) {
    total = total + sum((probf[, i] - real[, i]) ^ 2)
  }
  total / time_window
}
# equality constraint

# >=0

# >= -1



opt_simplified <- function(distance, pi_hat, real){
  r = dim(distance)[1]
  q = dim(distance)[2]
  time_window = dim(pi_hat)[2]
  
  getVIndexs <- function(distance){
    Avec <- matrix(0, r, q)
    indicators <- list()
    for (i in 1:q){
      if (min(distance[,i]) > 0){
        Avec[which(distance[,i]==min(distance[,i])), i] = 1
      }
    }
    for (i in 1:r){
      indicators[[i]] <- which(Avec[i, ] == 1)
    }
    indicators
  }
  construct_Q <- function(pi_hat, indicators){
    Qs <- lapply(indicators, function(indicator){
      pi_hat[indicator] %*% t(pi_hat[indicator])
    })
    bdiag(Qs)
  }
  construct_D <- function(pi_hat, indicators){
    D_mat <- NULL
    for (time in 1:time_window){
      Ds <- lapply(indicators, function(indicator){
        matrix(pi_hat[indicator, time], nrow=1)
      })
      D_mat <- rbind(D_mat, bdiag(Ds))
    }
    D_mat
  }
  construct_E <- function(distance, indicators){
    n <- sum(sapply(indicators, function(x){length(x)}))
    E <- NULL
    for (j in 1:dim(distance)[2]){
      if (min(distance[, j]) == 0) next;
      E_jRow <- matrix(0, ncol=n)
      currentVar = 0
      for (i in 1:length(indicators)) {
        E_jRow[which(indicators[[i]] == j) + currentVar] = 1
        currentVar = currentVar + length(indicators[[i]])
      }
      E <- rbind(E, E_jRow)
    }
    E
  }
  construct_res <- function(distance, solution, indicators){
    res <- replace(distance, distance==0, -1)
    res <- replace(res, distance>0, 0)
    currentVar = 1
    for (i in 1:length(indicators)){
      res[i, indicators[[i]]] = solution[currentVar: (currentVar + length(indicators[[i]]) - 1)]
      currentVar = currentVar + length(indicators[[i]])
    }
    res
  }
  indicators <- getVIndexs(distance)
  Q <- construct_Q(pi_hat, indicators)/time_window
  real <- as.vector(real)
  D <- construct_D(pi_hat, indicators)
  E <- construct_E(distance, indicators)
  n_eq <- dim(E)[1]
  n_var <- dim(E)[2]
  A <- rbind(E, 
             diag(rep(1, n_var)),
             -diag(rep(1, n_var)))
  
  Dmat <- 2 * Q
  dvec <- t(2 / time_window * t(real) %*% D)
  Amat <- t(A)
  bvec <- c(rep(1, n_eq),
            rep(0, n_var),
            rep(-1, n_var))
  solution <- solve.QP(Dmat, dvec, Amat, bvec, meq = n_eq)
  construct_res(distance, solution, indicators)
}
opt_simplified(distance, pi_hat, real)




## Binomail AR(1) fit

tpbin <- function(k,l,p,rho, n){
  beta <- p*(1-rho)
  alpha <- beta+rho
  
  tp <- 0
  for(j in c(max(0,k+l-n):min(k,l))){
    tp <- tp + dbinom(j,l,alpha)*dbinom(k-j,n-l,beta)
  }
  tp
}


#Log-likelihood of binomial AR(1) model:
llbar1 <- function(par,data, n){
  #par is vector (p,rho)
  T <- length(data)
  value <- -log(dbinom(data[1], n, par[1])) #full likelihood, otherwise use 0 here
  
  for(t in c(2:T)) {
    value <- value-log(tpbin(data[t], data[t-1], par[1], par[2], n))
  }
  value
}

binARforecast <- function(data, h, total=FALSE){
  T <- length(data)
  maxval = 1
  n = 1
  if (total){ 
    maxval = 2 
    n = 2
  }
  
  estml <- suppressWarnings(optim(c(sum(data>0)/T, 0.65), llbar1, method="L-BFGS-B", lower=c(0.0001,0.0001), upper=c(0.9999,0.9999), control=list(ndeps=c(1e-4,1e-4)), data=data, n=maxval, hessian=TRUE))
  pestml <- estml$par[[1]]
  rhoestml <- estml$par[[2]]
  
  forecasts <- array(0, c(h, maxval+1))
  for(i in c(1:h)){
    for(k in c(0:maxval)){
      forecasts[i,k+1] <- tpbin(k, data[T],pestml,rhoestml^h, maxval)
    }
  }
  forecasts <- t(apply(forecasts, 1, function(x){x/sum(x)}))
  forecasts[1,]
}


