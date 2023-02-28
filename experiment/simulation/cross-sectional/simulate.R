library(dplyr)
library(DiscreteRecon)
set.seed(42)


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
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)


cal_singleBase <- function(hist){
  # generate base forecast
  f1 <- binARforecast(hist[, 1], 1)
  f2 <- binARforecast(hist[, 2], 1)
  f3 <- binARforecast(rowSums(hist), 1, total = TRUE)
  list(f1, f2, f3)
}

library(dplyr)

cal_Base <- function(series, history_length, window_length){
  trainBase <- list()
  for (i in 1:window_length){
    trainBase[[i]] <- cal_singleBase(series[i:(i+history_length),])
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

# brier_score of reconciled forecast

res <- foreach(i=1:1000, .packages = c("DiscreteRecon"), 
               .export = c("simulate_series", "llbar1", "binARforecast", "tpbin")) %dopar% {
  output <- list()
  a1 = runif(1, 0.4, 0.5)
  a2 = runif(1, 0.3, 0.5)
  output$a <- c(a1, a2)
  series <- simulate_series(output$a, n)
  output$train <- try(cal_Base(series, history_length, window_length), silent = TRUE)
  output$test <- try(cal_Base(series[(window_length+1): n,], history_length, test_length), silent = TRUE)
  if (is.element('try-error', class(output$train)) | is.element('try-error', class(output$test))){
    print(class(output$train))
    output <- list()
    return(output)
  }
  
  train_dhts <- dhts(series[(history_length+1): (history_length+window_length),],
                     s_mat = s_mat,
                     domain)
  test_dhts <- dhts(series[(history_length+window_length+1): (history_length+test_length+window_length),],
                    s_mat = s_mat,
                    domain)
  dfr <- reconcile_train(output$train, train_dhts, step_wise = FALSE)
  dfr <- reconcile(dfr, output$test, meta=train_dhts$meta)
  
  td <- topdown.train(train_dhts)
  td <- reconcile(td, output$test)
  
  bu <- marginal2Joint(output$test, test_dhts$meta, method='bu')
  
  bs_vec <- function(x){ unname(c(x$series, sum(x$hierarchy))) }
  
  output$metric <- list(dfr=bs_vec(brier_score(dfr, test_dhts)),
                       bu=bs_vec(brier_score(bu, test_dhts)),
                       td=bs_vec(brier_score(td, test_dhts)),
                       base=bs_vec(brier_score(output$test, test_dhts)))
  output
}

# save results
saveRDS(res, 'res.rds')

# filter out errors
res <- Filter(function(x){!is.null(x$metric)}, res)


# summarise
accs_list <- list()
for (m in c("base", "bu", "td", "dfr")){
  accs_list[[m]] <- sapply(res, function(x){
    x$metric[[m]]
  })
}
accs_sum <- data.frame(lapply(accs_list, rowMeans))
row.names(accs_sum) <- c("y_0", "y_1", "y_2", "Y")
accs_sum

library(tsutils)



plot_data <- list()
for (i in 1:4){
  plot_data[[i]] <- sapply(accs_list, function(x){
    x[i,]
  })
}

pdf(file="mcb.pdf", width = 12, height = 4,
    pointsize = 16)
par(mfrow=c(1,3))
nemenyi(plot_data[[1]], plottype = "vmcb", 
        labels = c("Base", "Bottom-up", "Top-down", "DFR"), 
        main = "MCB Test for total series")
nemenyi(rbind(plot_data[[2]], plot_data[[3]]), plottype = "vmcb", 
        labels = c("Base", "Bottom-up", "Top-down", "DFR"), 
        main = "MCB Test for bottom series")
nemenyi(plot_data[[4]], plottype = "vmcb", 
        labels = c("Base", "Bottom-up", "Top-down", "DFR"), 
        main = "MCB Test for hierarchy")
dev.off()

