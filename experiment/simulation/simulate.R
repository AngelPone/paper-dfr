library(dplyr)
source('qp.R')

# Domain Setting
incoherent_domain = expand.grid(0:1, 0:1, 0:2)
incoherent_string = do.call(paste0, incoherent_domain)
coherent_domain = expand.grid(0:1, 0:1)
coherent_domain$Var3 = coherent_domain$Var1 + coherent_domain$Var2
coherent_string = do.call(paste0, coherent_domain)

r = nrow(coherent_domain)
q = nrow(incoherent_domain)

distance <- cal_distanceMatrix(incoherent_domain, coherent_domain)
colnames(distance) <- incoherent_string
rownames(distance) <- coherent_string

history_length = 100
window_length = 150
test_length = 30

n = history_length+window_length+test_length

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

cal_singleBase <- function(hist){
  # generate base forecast
  f1 <- glarma_base_forecast(hist[, 1])
  f2 <- glarma_base_forecast(hist[, 2])
  f3 <- glarma_base_forecast(rowSums(hist), total = TRUE)
  list(base=apply(expand.grid(f1, f2, f3), 1, prod),
       bu=apply(expand.grid(f1, f2), 1, prod)
       #prob=list(f1, f2, f3)
       )
}

cal_singleBase2 <- function(hist){
  # generate base forecast
  f1 <- binARforecast(hist[, 1], 1)
  f2 <- binARforecast(hist[, 2], 1)
  f3 <- binARforecast(rowSums(hist), 1, total = TRUE)
  list(base=apply(expand.grid(f1, f2, f3), 1, prod),
       bu=apply(expand.grid(f1, f2), 1, prod)
       #prob=list(f1, f2, f3)
  )
}


cal_Base <- function(series, history_length, window_length){
  library(dplyr)
  cal_singleBase <- cal_singleBase2
  trainBase <- foreach(i=1:(window_length), .export=c('binARforecast', 'llbar1', 'tpbin')) %dopar% 
    cal_singleBase(series[i:(i+history_length),])
  base <- sapply(trainBase, function(x){x$base})
  bu <- sapply(trainBase, function(x){x$bu})
  real <- cons_realDummy(series[(history_length+1):(window_length+history_length),], coherent_domain)
  list(base=base, real=real, bu=bu)
}

# brier_score of reconciled forecast

res <- list()
for (i in (length(res)+1):(length(res)+20)){
  res[[i]] <- list()
  a1 = runif(1, 0.4, 0.5)
  a2 = runif(1, 0.3, 0.5)
  res[[i]]$a <- c(a1, a2)
  series <- simulate_series(res[[i]]$a, n)
  res[[i]]$train <- try(cal_Base(series, history_length, window_length), silent = TRUE)
  res[[i]]$test <- try(cal_Base(series[(window_length+1): n,], history_length, test_length), silent = TRUE)
  if (is.element('try-error', class(res[[i]]$train)) | is.element('try-error', class(res[[i]]$test))){
    print(i)
    res[[i]] <- list()
    next
  }
  metric <- tibble(bs=numeric(), method=character(), lambda=numeric())
  res[[i]]$Ahat <- list()
  for (lambda in c(0, 0.001, 0.005, 0.01)){
    Amat <- 
      opt_fun(res[[i]]$train$base, res[[i]]$train$real, distance, lambda = lambda)
    metric <- metric %>% 
      add_row(lambda=lambda, method='bu', bs=brier_score(res[[i]]$test$bu, res[[i]]$test$real)) %>%
      add_row(lambda=lambda, method='reconciliation', bs=brier_score(Amat %*% res[[i]]$test$base, res[[i]]$test$real))
    Amat[Amat<1e-10] = 0
    res[[i]]$Ahat[[as.character(lambda)]] <- Amat
  }
  res[[i]]$series <- series
  res[[i]]$metric <- metric
}
res <- Filter(function(x){!is.null(x$metric)}, res)

accs <- res %>% lapply(function(x){x$metric}) %>% 
  do.call(rbind.data.frame, .)

accs %>% group_by(method, lambda) %>% 
  summarise(bs = mean(bs))
