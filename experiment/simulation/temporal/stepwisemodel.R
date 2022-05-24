setwd("simulation/temporal")
rm(list = ls())

devtools::load_all('../../../DiscreteRecon/')
source('./config.R')

# load pre-trained base forecasts
basef <- readRDS('simu_temporal_base.rds')


#' function to train single series
#' @param x x_train, x_test, y_train, y_test of single series
train_func <- function(x){
  y_train <- dhts(x$y_train, s_mat, domain)
  x_train <- x$x_train
  step_model <- reconcile.train(x_train, y_train)
}

models <- list()
for (i in 1:length(basef)){
  x <- basef[[i]]
  x_train <- x$x_train
  y_train <- dhts(x$y_train, s_mat, domain)
  models[[i]] <- reconcile.train(x_train, y_train, optimized = TRUE)
}

saveRDS(models, "stepwisemodels.rds")


# 
# concat_ <- function(indx){
#   x_train <- list()
#   for (i in 1:8){
#     x_train[[i]] <- do.call(rbind, lapply(basef[indx], function(x){x$x_train[[i]]}))
#   }
#   y_train <- do.call(rbind, lapply(basef[indx], function(x){x$y_train}))
#   list(x_train, y_train)
# }

