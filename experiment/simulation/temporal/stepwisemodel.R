# setwd("simulation/temporal")
# rm(list = ls())

#devtools::load_all('../../../DiscreteRecon/')
library(DiscreteRecon)
source('./config.R')

# load pre-trained base forecasts
basef <- readRDS('basef_temporal_2.rds')


#' function to train single series
#' @param x x_train, x_test, y_train, y_test of single series
train_func <- function(x){
  y_train <- dhts(x$y_train, s_mat, domain)
  x_train <- x$x_train
  step_model <- reconcile.train(x_train, y_train)
}


# models <- foreach(i = 1:length(basef), .packages = "DiscreteRecon") %dopar% {
#   x <- basef[[i]]
#   x_train <- x$x_train
#   y_train <- dhts(x$y_train, s_mat, domain)
#   reconcile_train(x_train, y_train, optimized = TRUE)
# }
models <- list()
for(i in 1:length(basef)){
  x <- basef[[i]]
  x_train <- x$x_train
  y_train <- dhts(x$y_train, s_mat, domain)
  models[[i]] <- reconcile_train(x_train, y_train, optimized = TRUE)
}

saveRDS(models, "stepwisemodels2.rds")


