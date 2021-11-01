library(Matrix)
library(dplyr)
library(quadprog)


construct_Q <- function(r, q, pi_hat){
  # r: cardinality of reconciled probability
  # q: cardinality of probability
  # pi_hat: q*T matrix
  Q <- pi_hat %*% t(pi_hat)
  bdiag(replicate(r, Q, simplify = FALSE))
}

# test
# Q <- construct_Q(r, q, pi_hat)
# A = rnorm(r*q) %>% matrix(r, q)

## result of matrix form
# Avec <- as.vector(t(A))
# pi_z <- matrix(0, r*q, t_window)
# for (i in 1:t_window){
#  pi_z[,i] <- rep(pi_hat[,i], r) * rep(Z[,i], each=q)
# }

# D <- lambda*as.vector(t(d)) - 2*apply(pi_z, 1, sum)

# matrix_result <- t(Avec) %*% as.matrix(Q) %*% Avec + t(D) %*% Avec + t(as.vector(t(Z))) %*% as.vector(t(Z))

# raw result
# S = 0
#  for (t in 1:100){
#   for (k in 1:r) {
#     S = S + (sum(A[k,]*pi_hat[,t]) - Z[k,t])^2
#   }
# }
# S = S + lambda*sum(d*A)
# S



#QR with solve.qr

opt_fun <- function(pi_hat, real, distance, lambda){
  # pi_hat: q * T
  # real: r * T
  # distance: r * q
  r = dim(real)[1]
  q = dim(pi_hat)[1]
  time_window = dim(pi_hat)[2]
  Dmat <- 2*construct_Q(r, q, pi_hat)
  pi_z <- matrix(0, r*q, t_window)
  for (i in 1:time_window){
    pi_z[,i] <- rep(pi_hat[,i], r) * rep(Z[,i], each=q)
  }
  dvec <-  2*apply(pi_z, 1, sum) - lambda*as.vector(t(distance))
  
  A1 <- diag(q)
  for (i in 1:(r-1)){
    A1 <- cbind(A1, diag(q))
  }
  b1 <- rep(1, q)
  
  A2 <- diag(r*q)
  b2 <- rep(0, r*q)
  
  A3 <- -diag(r*q)
  b3 <- rep(-1, r*q)
  Amat <- t(rbind(A1, A2, A3))
  bvec <- c(b1, b2, b3)
  
  solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = r)
  
  t(matrix(solution$solution, 12, 4))
}



simulate_series <- function(param_a, n, intercept = 0.0,
                            innovation_mu=c(0, 0), 
                            innovation_Sigma=matrix(c(0.1, 0.05, 0.05, 0.1), 2, 2)){
  # simulate time series using State Space model
  error_term <- MASS::mvrnorm(n+1, innovation_mu, innovation_Sigma)
  # latent state
  for (i in 2:dim(error_term)[1]){
    error_term[i,] = intercept + param_a * error_term[i-1,] + error_term[i,]
  }
  # pi
  pi = 1 / (1 + 1/exp(error_term))
  ifelse(pi>0.5, 1, 0)[2:(n+1), ]
}

glarma_base_forecast <- function(series, total=FALSE){
  if (total){
    series <- cbind(series, 2-series)
  }else{
    series <- cbind(series, 1-series)
  }
  model <- glarma::glarma(series, X=matrix(1, nrow = dim(series)[1]), 
                          phiLags = c(1), phiInit = c(0.5), type = 'Bin')
  fore <- glarma::forecast(model, h=1, newdata = matrix(1))
  
  if (total){
    pi <- 1/(1+1/exp(fore$W))
    return(c((1-pi)^2, 2*pi*(1-pi), pi^2))
  }else{
    return(c(1-fore$mu, fore$mu))
  }
}


# equality constraint

# >=0

# >= -1




# computation consideration
# getVIndexs <- function(distance){
#   Avec <- rep(0, r*q)
#   for (i in 1:q){
#     if (min(distance[,i]) == 0){
#       Avec[(which.min(distance[,i])-1)*q+i] = 1
#     } else {
#       Avec[(which(d[,i]==min(d[,i]))-1)*q+i] = -1
#     }
#   }
#   Avec
# }
# 
# indicator <- getVIndexs(d)
# # select elements to be optimized and their corresponding Dmat, Amat, bvec, dvec
# variable_index <- which(indicator==-1)
# Dmat_opt <- 2*Q[variable_index, variable_index]
# dvec_opt <- dvec[variable_index]
# A1_opt <- A1[, variable_index]
# b1 <- rep(1, length(which(!apply(A1_opt, 1, function(x){all(x==0)}))))
# A1_opt <- A1_opt[which(!apply(A1_opt, 1, function(x){all(x==0)})),]
# A2_opt <- diag(length(variable_index))
# b2 <- rep(0, length(variable_index))
# A3_opt <- -diag(length(variable_index))
# b3 <- rep(-1, length(variable_index))
# Amat_opt <- rbind(A1_opt, A2_opt, A3_opt)
# bvec_opt <- c(b1,b2,b3)
# 
# # optimize
# solution_opt <- solve.QP(Dmat_opt, dvec_opt, t(Amat_opt), bvec_opt, 
#                          meq = dim(A1_opt)[1])
# 
# # fill result and back to matrix
# indicator[variable_index] = solution_opt$solution
# A_solved <- t(matrix(indicator, q, r))
# 
# # check result
# apply(A_solved %*% pi_hat, 2, sum)


