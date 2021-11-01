library(dplyr)
source('qp.R')

# Domain Setting
incoherent_domain = expand.grid(0:1, 0:1, 0:2)
coherent_domain = expand.grid(0:1, 0:1)
coherent_domain$Var3 = coherent_domain$Var1 + coherent_domain$Var2

r = nrow(coherent_domain)
q = nrow(incoherent_domain)

distance = matrix(0, nrow = r, ncol = q)
for (j in 1:r){
  for (k in 1:q){
    distance[j, k] = sum(abs(coherent_domain[j,] - incoherent_domain[k,]))
  }
}

history_length = 300
window_length = 60
test_length = 30

n = history_length+window_length+test_length

foo <- simulate_series(param_a = runif(1, 0.5, 0.9), 
                       n = n, 
                       intercept = -0.25,
                       innovation_mu = c(0, 0),
                       innovation_Sigma = matrix(c(0.1, 0.05, 0.05, 0.2), 2, 2))

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

loop_FUN <- function(i){
  f1 <- glarma_base_forecast(foo[i:(i+history_length-1), 1])
  f2 <- glarma_base_forecast(foo[i:(i+history_length-1), 2])
  f3 <- glarma_base_forecast(apply(foo[i:(i+history_length-1),], 1, sum), total = TRUE)
  list(pred=apply(expand.grid(f1, f2, f3), 1, prod),
       real=(foo[i+history_length,]),
       f1=f1,
       f2=f2,
       f3=f3)
}

loop_result <- function(res, l = window_length){
  pred <- sapply(res, function(x){x$pred})
  real <- sapply(res, function(x){x$real}) %>% t()
  # dummy real observations matrix
  real_mat <- matrix(0, l, 4)
  for (i in 1:l){
    real_mat[i, which((coherent_domain$Var1 == real[i,1]) & 
                     (coherent_domain$Var2 == real[i,2]))] <- 1
  }
  list(pred=pred, real=t(real_mat))
}

train_res <- foreach(i=1:(window_length)) %dopar% 
  loop_FUN(i)
train_res <- loop_result(train_res)

# train model
lambda = 1

Ahat <- opt_fun(train_res$pred, train_res$real, distance, lambda = lambda)
Ahat[Ahat < 1e-6] = 0


# test
test_base <- foreach(i=(window_length+1):(window_length+test_length)) %dopar% 
  loop_FUN(i)
test_base <- loop_result(test_base, l=test_length)

base <- apply(test_base$pred, 2, which.max)
real <- apply(test_base$real, 2, function(x){which(x==1)})
reconciled <- apply(Ahat %*% test$pred, 2, which.max)

sum(reconciled == real)/length(real)
sum(base == real)/length(real)

