source('./config.R')
library(DiscreteRecon)

basef <- readRDS('basef.rds')
# basef[[i]]: x_train, y_train, x_test, y_test
# x_train: list{8, 526 * 8}
# y_train: 526 * 7
# x_test: list{8, 21 * 8}
# y_test: 21 * 7


## Step1: train the reconciliation model #########
start = Sys.time()
cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)

step_models <- foreach::foreach(bf = iterators::iter(basef), .packages = "DiscreteRecon", .errorhandling = "pass") %dopar% {
  x_train <- bf$x_train
  y_train <- dhts(bf$y_train, smatrix, domain)
  reconcile_train(x_train, y_train, optimized = FALSE, lambda = 1,
                  step_wise=TRUE)
}
end = Sys.time()
sprintf("step_models 用时 %2d 小时", (end - start)/3600)

saveRDS(step_models, "step_models.rds")



## Step2: 调整预测值 #####
l <- lapply(step_models, function(x){length(x)})
idx <- which(l == 6)
step_models <- step_models[idx]
basef <- basef[idx]

start = Sys.time()
res <- foreach(i=1:length(step_models), .packages = "DiscreteRecon", .errorhandling = "pass") %dopar% {
  modeli <- step_models[[i]]
  x_test <- basef[[i]]$x_test
  y_test <- dhts(basef[[i]]$y_test, smatrix, domain)
  dist <- reconcile(modeli, x_test)
  list(dist=dist, y=y_test)
}
end = Sys.time()
sprintf("step_models 用时 %2d 小时", (end - start)/3600)
saveRDS(res, 'SFP.rds')


## Step3: Top-down reconciliation####
tdres <- foreach(i=1:length(step_models), .packages = "DiscreteRecon", .errorhandling = "pass") %dopar%{
  modeli <- topdown.train(dhts(basef[[i]]$y_train, smatrix, domain))
  x_test <- basef[[i]]$x_test
  y_test <- dhts(basef[[i]]$y_test, smatrix, domain)
  dist <- reconcile(modeli, x_test)
  list(dist=dist, y=y_test)
}

saveRDS(tdres, 'td.rds')
