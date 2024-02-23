library(tscount)
library(foreach)
setwd("experiment/dc_crime/")
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

TEST_START <- c(2022, 1)
TEST_END <- c(2022, 20)
TRAIN_END <- c(2021, 53)
AGGREGATE_LENGTH <- 4
WINDOW_SIZE <- 53 * 1
VALID_SIZE <- 4

temporalAgg <- function(x){
  x <- x[(length(x) - length(x) %/% AGGREGATE_LENGTH * AGGREGATE_LENGTH + 1):length(x)]
  idx <- rep(1:(length(x) %/% AGGREGATE_LENGTH), each = AGGREGATE_LENGTH)
  unname(sapply(split(unclass(x), idx), sum))
}

rolling_train <- function(tss){
  stopifnot(is.ts(tss))
  tss <- window(tss, end = TEST_END)
  if (length(tss) >= WINDOW_SIZE + VALID_SIZE){
    rolling_size <- length(tss) - WINDOW_SIZE - VALID_SIZE + 1
    output <- lapply(1:rolling_size, function(x){
      s1 <- time(tss)[1]
      s2 <- time(tss)[x + WINDOW_SIZE - 1]
      s3 <- time(tss)[x + WINDOW_SIZE]
      s4 <- time(tss)[x + WINDOW_SIZE + VALID_SIZE - 1]
      list(train=window(tss, start = s1, end=s2) , test=window(tss, start=s3, end=s4))
    })
  }
}

data <- readRDS('data/data.rds')


basef_series <- function(x){
  x <- rolling_train(x)
  fcasts <- list()
  test_lst <- list()
  for (i in seq_along(x)) {
    train <- x[[i]]$train
    test <- x[[i]]$test
    mdl <- tsglm(train, model = list(past_obs = 1:3, past_mean = 1:4), link = "identity",
          distr = "poisson")
    mean_bottom <- predict(mdl, n.ahead = 4)$pred
    
    # total
    
    x_total <- temporalAgg(train)
    mdl <- tsglm(x_total, model = list(past_obs = 1:2, past_mean = 1:2), link = "identity",
                 distr = "poisson")
    fcasts[[i]] <- c(predict(mdl, n.ahead = 1)$pred, mean_bottom)
    test_lst[[i]] <- test
  }
  
  list(fcasts = do.call(rbind, fcasts), y = do.call(rbind, test_lst))
}

basef <- foreach(dt = iterators::iter(data), .packages = "tscount", .errorhandling = "pass") %dopar% {
  basef_series(dt)
}

saveRDS(basef, "data/basef.rds")



