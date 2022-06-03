# file including config of simulation experiments in temporal setting.
library(dplyr)
library(foreach)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
# summing matrix
s_mat <- rbind(rep(1, 7), diag(rep(1, 7)))

# domain of hierarchy
domain <- rbind(rep(0, 7), rep(1, 7))

## utility functions
logit <- function(x) exp(x)/(1+exp(x))
temporalAgg <- function(x){
  mult <- floor(length(x)/7)*7
  x <- x[(length(x) - mult + 1) : length(x)]
  sapply(split(ifelse(x>=1, 1, 0), ceiling((1:mult)/7)), sum)
}


reduction2lst <- function(series){
  res <- list()
  for (i in 1:(length(series) - 6 - lb_window)){
    res[[i]] <- list(
      x = series[i: (i+lb_window-1)],
      y = series[(i+lb_window) : (i + lb_window + 6)]
    )
  }
  res
}
reduction2df <- function(series){
  x <- NULL
  for (i in 1:(length(series) - 6)){
    x <- rbind(x, series[i:(i+6)])
  }
  x <- data.frame(x)
  colnames(x)[7] <- "y"
  x
}

train_test_prepare <- function(ts_lst){
  tmpf1 <- function(x){
    cs = reduction2df(x$x)
    model <- glm(y~., data=cs)
    pred <- c()
    new_x <- cs[dim(cs)[1], 2:7]
    colnames(new_x) <- paste0("X", 1:6)
    for (i in 1:7){
      pred <- c(pred, predict(model, new_x))
      new_x[1, 1:5] <- new_x[1, 2:6]
      new_x[1, 6] <- ifelse(pred[length(pred)] > 0.5, 1, 0)
    }
    names(pred) = NULL
    list(x=pred, y=x$y)
  }
  logit_res <- foreach(x=iterators::iter(ts_lst), .export = c("reduction2df"), .errorhandling = "pass") %dopar% tryCatch(tmpf1(x))
  tmpf <- function(x){
    cs = temporalAgg(x$x)
    model <- tryCatch(tscount::tsglm(cs, list(past_obs=1, past_mean=1:3), link=c("identity"), distr = c("poisson")))
    if ('error' %in% class(model)){
      return(model)
    }
    pred <- dpois(0:6, predict(model)$pred)
    pred <- c(pred, 1-sum(pred))
    list(x=pred, y=sum(x$y))
  }
  total_res <- foreach(x=iterators::iter(ts_lst), .export = c("temporalAgg"), .errorhandling = "pass") %dopar% tmpf(x)
  
  wrong_idx <- which(sapply(total_res, function(x){'error' %in% class(x)}))
  if (length(wrong_idx) > 0){
    total_res <- total_res[-wrong_idx]
    logit_res <- logit_res[-wrong_idx]
  }
  
  x_train = list()
  x_train[[1]] = t(sapply(total_res[1:train_size], function(x){x$x}))
  
  x_test = list()
  x_test[[1]] = t(sapply(total_res[(train_size+1):length(total_res)], function(x){x$x}))
  colnames(x_train[[1]]) <- 0:7
  
  for (i in 2:8){
    x_train[[i]] <- t(sapply(logit_res[1:train_size], function(x){
      c(1-x$x[i-1], x$x[i-1])
    }))
    colnames(x_train[[i]]) <- c(0, 1)
    x_test[[i]] <- t(sapply(logit_res[(train_size+1):length(logit_res)], function(x){
      c(1-x$x[i-1], x$x[i-1])
    }))
    colnames(x_test[[i]]) <- c(0, 1)
  }
  
  y_train = t(sapply(logit_res[1:train_size], function(x){x$y}))
  y_test = t(sapply(logit_res[(train_size+1):length(logit_res)], function(x){x$y}))
  list(x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test, error_idx=wrong_idx)
}

