library(tscount)
library(rslurm)

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

data_rolling <- lapply(data, rolling_train)
# for (i in seq_along(data_rolling)){
#   data_rolling[[i]] <- lapply(data_rolling[[i]], function(x){
#     x$id <- i
#     x
#   })
# }


# filter all zero series
data_rolling <- lapply(data_rolling, function(x){Filter(function(y){any(y$train>0)}, x)})
data_rolling <- do.call(c, data_rolling)

# get train and test
train_rolling <- Filter(function(x){max(time(x$train)) < 2021.999}, data_rolling)
test_rolling <- Filter(function(x){max(time(x$train)) >= 2021.999}, data_rolling)


Filter(function(x){}, data_rolling)
# saveRDS(list(train=train_rolling, test=test_rolling), "data/rolling_data.rds")


cal_basef <- function(obj){
  train <- obj$train
  mdl <- tsglm(train, model = list(past_obs = 1:3, past_mean = 1:4), link = "identity",
          distr = "poisson")
  obj$mean <- predict(mdl, n.ahead = 4)$pred
  obj
}

cal_basef_total <- function(obj){
  x <- temporalAgg(obj$train)
  mdl <- tsglm(x, model = list(past_obs = 1:2, past_mean = 1:2), link = "identity",
               distr = "poisson")
  obj$mean <- c(predict(mdl, n.ahead = 1)$pred, obj$mean)
  obj
}



slurm_map(train_rolling,
          function(x){cal_basef_total(cal_basef(x))},
          global_objects = ls(),
          jobname = "dfr_base_train",
          nodes = 2,
          cpus_per_node = 32,
          submit = TRUE)
slurm_map(test_rolling, 
          function(x){cal_basef_total(cal_basef(x))},
          jobname = "dfr_base_test",
          global_objects = ls(),
          nodes = 1,
          cpus_per_node = 32,
          submit = TRUE)

