setwd("experiment/simulation")
source('temporal/config.R')
library(DiscreteRecon)

basef <- readRDS('results/basef.rds')
# basef[[i]]: x_train, y_train, x_test, y_test
# x_train: list{8, 526 * 8}
# y_train: 526 * 7
# x_test: list{8, 21 * 8}
# y_test: 21 * 7


cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

ht <- dhier(smatrix, domain)
step_models <- foreach::foreach(bf = iterators::iter(basef), .packages = "DiscreteRecon", .errorhandling = "pass") %dopar% {
  bf_train <- bf$x_train
  obs_train <- bf$y_train
  mdl <- dfr(ht, method = "sdfr", obs_train = obs_train, bf_train = bf_train)
  mdl_bu <- dfr(ht, "bu")
  mdl_td <- dfr(ht, "td", obs_train)
  mdl_emp <- dfr(ht, "emp", obs_train)
  
  
  bf_test <- bf$x_test
  sdfr = reconcile(mdl, bf_test)
  bu = reconcile(mdl_bu, bf_test)
  td = reconcile(mdl_td, bf_test)
  emp = reconcile(mdl_emp, bf_test)
  list(
    fcasts = list(sdfr=sdfr, bu=bu, td=td, emp=emp, base=bf_test),
    y = bf$y_test
  )
}

saveRDS("results/temporal_recf.rds")
