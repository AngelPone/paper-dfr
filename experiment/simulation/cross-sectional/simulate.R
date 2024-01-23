library(dplyr)
library(DiscreteRecon)
set.seed(42)
setwd("experiment/simulation")

simulate_series <- function(param_a,
                            n,
                            intercept = 0.0,
                            innovation_mu = c(0, 0),
                            innovation_Sigma = matrix(c(0.1, 0.05, 0.05, 0.1), 2, 2)) {
  # simulate time series using State Space model
  error_term <- MASS::mvrnorm(n + 1, innovation_mu, innovation_Sigma)
  # latent state
  for (i in 2:NROW(error_term)) {
    error_term[i, ] = error_term[i, ] + param_a * error_term[i - 1, ]
  }
  # pi
  pi = 1 / (1 + 1 / exp(error_term))
  ifelse(pi > 0.5, 1, 0)[2:(n + 1),]
}

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




history_length = 150
window_length = 300
test_length = 30

n = history_length+window_length+test_length


library(foreach)
num.cores <- 8

if (exists("num.cores")) {
  cl <- parallel::makeCluster(num.cores)
} else {
  cl <- parallel::makeCluster(1)
}

doParallel::registerDoParallel(cl)


cal_singleBase <- function(hist){
  f1 <- binARforecast(hist[, 1], 1)
  f2 <- binARforecast(hist[, 2], 1)
  f3 <- binARforecast(rowSums(hist), 1, total = TRUE)
  list(f1, f2, f3)
}

library(dplyr)

cal_Base <- function(series, window_length, history_length){
  trainBase <- list()
  t <- NROW(series)
  for (i in window_length : 1){
    trainBase[[(window_length - i + 1)]] <- 
      cal_singleBase(series[(t-i-history_length+1):(t - i),])
  }
  base <- list()
  base[[1]] <- t(sapply(trainBase, function(x){x[[3]]}))
  colnames(base[[1]]) <- c("0", "1", "2")
  base[[2]] <- t(sapply(trainBase, function(x){x[[1]]}))
  colnames(base[[2]]) <- c("0", "1")
  base[[3]] <- t(sapply(trainBase, function(x){x[[2]]}))
  colnames(base[[3]]) <- c("0", "1")
  base
}

s_mat <- rbind(c(1, 1), diag(2))
domain <- rbind(c(0, 0), c(1, 1))

res <- foreach(i=1:1000, .packages = c("DiscreteRecon"), 
               .export = c("simulate_series", "llbar1", "binARforecast", "tpbin")) %dopar% {
  output <- list()
  a1 = runif(1, 0.4, 0.5)
  a2 = runif(1, 0.3, 0.5)
  series <- simulate_series(c(a1, a2), n)
  bf <- try(cal_Base(series, window_length + test_length, history_length), silent = TRUE)
  if (is.element('try-error', class(output$train)) | is.element('try-error', class(output$test))){
    print(class(output$train))
    output <- list()
    return(output)
  }
  bf_train <- lapply(bf, function(x){ x[1:window_length,] })
  bf_test <- lapply(bf, function(x){ x[(window_length+1):(window_length+test_length),] })
  list(series=series, 
       fcasts = list(base_train=bf_train, base_test=bf_test))
}

saveRDS(res, 'results/cs-res.rds')

res <- foreach(dt=iterators::iter(res), .packages = c("DiscreteRecon")) %dopar% {
  obs <- dt$series[n - ((window_length + test_length):1) + 1,]
  
  obs_train <- obs[1:window_length,]
  obs_test <- obs[(window_length+1):(window_length + test_length),]
  
  bf_train <- dt$fcasts$base_train
  bf_test <- dt$fcasts$base_test
  
  ht <- dhier(s_mat, domain)
  mdl <- dfr(ht, "dfr", bf_train = bf_train, obs_train = obs_train)
  recf <- reconcile(mdl, bf_test)
  
  td <- dfr(ht, "td", obs_train = obs_train)
  td <- reconcile(td, bf_test)
  
  bu <- dfr(ht, "bu")
  bu <- reconcile(bu, bf_test)
  
  emp <- dfr(ht, method = "emp", obs_train = obs_train)
  emp <- reconcile(emp, obs_test, h = test_length)
  
  bs_vec <- function(x){ unname(c(x$series, sum(x$hierarchy))) }
  output <- list()
  
  output$metric <- list(dfr=bs_vec(brier_score(recf, obs_test, ht)),
                        bu=bs_vec(brier_score(bu, obs_test, ht)),
                        td=bs_vec(brier_score(td, obs_test, ht)),
                        base=bs_vec(brier_score(bf_test, obs_test, ht)),
                        emp=bs_vec(brier_score(emp, obs_test, ht)))
  output$fcasts <- list(dfr=recf, bu=bu, td=td, emp=emp, 
                        base_train=bf_train,
                        base_test=bf_test)
  output$series <- dt$series
  output
}

saveRDS(res, 'results/cs-res.rds')

