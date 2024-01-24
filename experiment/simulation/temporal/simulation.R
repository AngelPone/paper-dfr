library(gratis)
library(feasts)
library(tsibble)
setwd("experiment/simulation")
source("temporal/config.R")

generate_bts <- function(x){
  d <- rep(0, 7)
  d[order(rbeta(7, 1:7, 4), decreasing = TRUE)[1:x]] <- 1
  d
}

simulate_temporal <- function(time_window, times){
  total <- arima_model(4, p=3, d=0, D=0) %>% 
    generate(time_window + 100, times) %>%
    filter(index >= make_yearquarter(26, 1)) %>%
    group_by_key() %>%
    mutate(value = (value - min(value))/(max(value) - min(value)) * 2.5 +  2) %>%
    ungroup() %>%
    mutate(value = rpois(nrow(.), value)) %>%
    mutate(obs=ifelse(value > 7, 7, value))
  
  tmp <- function(x){
    matrix(sapply(as.list(x), generate_bts, simplify = TRUE), nrow = 1)[1,]
  }
  bts <- group_by_key(total) %>%
    group_map(~ tmp(.x$obs))
}

train_test_prepare <- function(ts_lst){
  tmpf1 <- function(x){
    cs = reduction2df(x$x)
    model <- glm(y~., data=cs)
    pred <- c()
    new_x <- cs[dim(cs)[1], 2:13]
    colnames(new_x)[1:6] <- paste0("X", 1:6)
    week_idx <- diag(7)[,1:6]
    for (i in 1:7){
      new_x[1, 7:12] <- week_idx[i,]
      pred <- c(pred, predict(model, new_x))
      new_x[1, 1:5] <- new_x[1, 2:6]
      new_x[1, 6] <- ifelse(pred[length(pred)] > 0.5, 1, 0)
    }
    names(pred) = NULL
    list(x=pred, y=x$y)
  }
  logit_res <- foreach(x=iterators::iter(ts_lst), .export = c("reduction2df"), .errorhandling = "pass") %do% tryCatch(tmpf1(x))
  tmpf <- function(x){
    cs = temporalAgg(x$x)
    model <- tryCatch(tscount::tsglm(cs, list(past_obs=1:3, past_mean=1:3), link=c("identity"), distr = c("poisson")))
    if ('error' %in% class(model)){
      return(model)
    }
    pred <- dpois(0:6, predict(model)$pred)
    pred <- c(pred, 1-sum(pred))
    list(x=pred, y=sum(x$y))
  }
  total_res <- foreach(x=iterators::iter(ts_lst), .export = c("temporalAgg"), .errorhandling = "pass") %do% tmpf(x)
  
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
  colnames(x_test[[1]]) <- 0:7
  
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
  week_idx <- lapply(1:(length(series)/7), function(x){diag(7)[,1:6]}) %>%
    do.call(rbind, .)
  x <- cbind(x, week_idx[7:length(series),])
  x <- data.frame(x)
  
  colnames(x)[7] <- "y"
  x
}


# simulate time series
series <- simulate_temporal(time_window, nsim)

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

# generate base forecast and prepare train set and test set
basef_temporal <- foreach(s=iterators::iter(series), 
                          .packages = c("dplyr", "foreach")) %dopar%
  train_test_prepare(reduction2lst(s))

saveRDS(basef_temporal, paste0("results/basef.rds"))
