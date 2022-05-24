# Simulate Daily Time Series
library(dplyr)
library(foreach)
setwd('simulation/temporal')
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

# config
lb_window = 300
train_size = 500
test_size = 21
h = 7
simulation_size = 1000

time_window = train_size + lb_window + test_size + 6

## utility functions
logit <- function(x) exp(x)/(1+exp(x))
temporalAgg <- function(x){
  mult <- floor(length(x)/7)*7
  x <- x[(length(x) - mult + 1) : length(x)]
  sapply(split(ifelse(x>=1, 1, 0), ceiling((1:mult)/7)), sum)
}

# random_params <- function(){
#   a <- abs(rnorm(1))
#   a <- c(a, a * 0.7^(1:5))
#   a/sum(a) * 0.9
# }
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



# simulate function to simulate daily time series
simfunc <- function(time_window, n){
  a = tsintermittent::simID(n, time_window, idi = 2, level = 1.5)
  a[a>1] = 1
  a
}

# function to construct train and test model
train_test_prepare <- function(ts_lst){
  logit_res <- lapply(ts_lst, function(x){
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
  })
  tmpf <- function(x){
    cs = temporalAgg(x$x)
    model <- tscount::tsglm(cs, list(past_obs=1:3, past_mean=1:3), link=c("identity"), distr = c("poisson"))
    pred <- dpois(0:6, predict(model)$pred)
    pred <- c(pred, 1-sum(pred))
    list(x=pred, y=sum(x$y))
  }
  total_res <- foreach(x=iterators::iter(ts_lst), .export = c("temporalAgg")) %dopar% tmpf(x)
  
  x_train = list()
  x_train[[1]] = t(sapply(total_res[1:500], function(x){x$x}))
  
  y_train = list()
  x_test = list()
  x_test[[1]] = t(sapply(total_res[501:length(total_res)], function(x){x$x}))
  colnames(x_train[[1]]) <- 0:7
  y_test = list()
  for (i in 2:8){
    x_train[[i]] <- t(sapply(logit_res[1:500], function(x){
      c(1-x$x[i-1], x$x[i-1])
    }))
    colnames(x_train[[i]]) <- c(0, 1)
    x_test[[i]] <- t(sapply(logit_res[501:length(logit_res)], function(x){
      c(1-x$x[i-1], x$x[i-1])
    }))
    colnames(x_test[[i]]) <- c(0, 1)
  }
  y_train = t(sapply(logit_res[1:500], function(x){x$y}))
  y_test = t(sapply(logit_res[501:length(logit_res)], function(x){x$y}))
  list(x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test)
}


# simulate series
series <- simfunc(time_window, simulation_size)

tt <- list()
for (i in 1:100){
  print(i)
  tt[[i]] = train_test_prepare(reduction2lst(series[,i]))
}

saveRDS(tt, 'sim_temporal_base.rds')


