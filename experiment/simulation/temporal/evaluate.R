devtools::load_all('../../../DiscreteRecon/')
source('./config.R')
#library(DiscreteRecon)
library(foreach)

models <- readRDS("stepwisemodels2.rds")
basef <- readRDS('basef_temporal_2.rds')
# basef[[i]]: x_train, y_train, x_test, y_test
# x_train: list{8, 500 * domain}
# y_train: 500 * 7
# x_test: list{8, 21 * domain}
# y_test: 21 * 7


# 
nsim = length(basef)

bs <- matrix(0, nsim, 2)
res <- list()
res <- foreach(i=1:nsim, .export = c("models", "basef"), .packages = "DiscreteRecon") %dopar% {
  modeli <- models[[i]]
  x_test <- basef[[i]]$x_test
  y_test <- dhts(basef[[i]]$y_test, s_mat, domain)
  dist <- reconcile(modeli, x_test)
  list(dist=dist, y=y_test)
}


# brier score with random reconciliation matrix
random_res <- function(x){
  A = matrix(rnorm(2^7 * 2^7 * 8), nrow = 2^7)
  A = apply(A, 2, function(x){abs(x)/sum(abs(x))})
  t(A %*% t(x))
}
random_bs <- t(sapply(as.list(1:nsim), function(i){
  x <- marginal2Joint(basef[[i]]$x_test)
  brier_score(random_res(x), res[[i]]$y)
}))



train_bs <- read.csv('run2.log')
train_bs %>% group_by(unknowns) %>% summarise(across(.fns = mean))
c(colMeans(test_bs), mean(random_bs))

# function to calculate brier score for each level
bs <- function(res, basef){
  domain <- res$y$domain$coherent_domain
  y <- res$y$bts
  tmp_bs <- function(dist, dummies){
    sum((dist - dummies)^2)/dim(dist)[1]
  }
  val_dum <- function(x, d){
    z <- matrix(0, length(x), d)
    for (i in 1:length(x)){
      z[i, x[i] + 1] <- 1
    }
    z
  }
  bss <- tmp_bs(Joint2Marginal(res$dist, domain, 1), val_dum(rowSums(y), 8))
  for(i in 2:dim(domain)[2]){
    bss <- c(bss, tmp_bs(Joint2Marginal(res$dist, domain, i), val_dum(y[,i-1], 2)))
  }
  
  output <- list(rec = bss)
  # base
  bss <- tmp_bs(basef$x_test[[1]], val_dum(rowSums(y), 8))
  for (i in 2:dim(domain)[2]){
    bss <- c(bss, tmp_bs(basef$x_test[[i]], val_dum(y[,i-1], 2)))
  }
  output$base <- bss
  
  # bu
  bss[1] <- tmp_bs(marginal2Sum(basef$x_test[2:8], res$y$domain$domain_bts), val_dum(rowSums(y), 8))
  output$bu <- bss
  output
}
# function to calculate point forecast of each series
point <- function(x, basef){
  point <- point_forecast(x$dist, x$y$domain$coherent)
  base_point_lower <- sapply(basef$x_test[2:8], function(x){x[,2]})
  base_point_upper <- apply(basef$x_test[[1]], 1, function(x){sum(x*(0:7))})
  bu_point <- rowSums(base_point_lower)
  list(bu=cbind(bu_point, base_point_lower), base=cbind(base_point_upper, base_point_lower), rec=point)
}


# res: res[[i]]
# basef: basef[[i]]
# function to compute point metrics (i.e., RMSE and MAE)
point_metrics <- function(i){
  pointf <- point(res[[i]], basef[[i]])
  y = res[[i]]$y$bts
  total <- rowSums(y)
  rmse <- function(x, y){sqrt(mean((x-y)^2))}
  mae <- function(x, y){mean(abs(x-y))}
  rmses <-foreach(method = c('bu', 'base', 'rec'), .combine=rbind) %do% {
    c(rmse(pointf[[method]][,2:8], y), rmse(pointf[[method]][,1], y))
  }
  maes <- foreach(method = c('bu', 'base', 'rec'), .combine=rbind) %do% {
    c(mae(pointf[[method]][,2:8], y), mae(pointf[[method]][,1], y))
  }
  d <- rbind(rmses, maes) %>% data.frame()
  colnames(d) <- c('level1', 'level2')
  d$method <- rep(c('bu', 'base', 'rec'), 2)
  d$metric <- rep(c("rmse", "mae"), each=3)
  d$index <- i
  d
}

# point metric
metric_point <- foreach(i=1:100, .combine = rbind) %do% point_metrics(i)
metric_point %>% group_by(metric, method) %>% summarise(across(starts_with('level'),.fns=mean))

# brier score
metric_bs_level <- foreach(i=1:100) %do% bs(res[[i]], basef[[i]])
metric_bs_level_bu <- do.call('rbind', lapply(metric_bs_level, function(x){x$bu}))
metric_bs_level_base <- do.call('rbind', lapply(metric_bs_level, function(x){x$base}))
metric_bs_level_rec <- do.call('rbind', lapply(metric_bs_level, function(x){x$rec}))
metric_bs_hierarchy <- t(sapply(as.list(1:nsim), function(i){
  x <- res[[i]]
  basef <- basef[[i]]$x_test
  c(brier_score(x$dist, x$y), brier_score(marginal2Joint(basef[2:8]), x$y))
}))

bss <- matrix(0, 3, 3)
rownames(bss) <- c("Base", "BU", "SDFR")
colnames(bss) <- c("Bottom", "Total", "Hierarchy")
bss[3, 3] <- colMeans(metric_bs_hierarchy)[1]
bss[2, 3] <- colMeans(metric_bs_hierarchy)[2]
bss[1, 3] <- NA
bss[3, 2:1] <- c(colMeans(metric_bs_level_rec)[1], mean(colMeans(metric_bs_level_rec)[2:8]))
bss[2, 2:1] <- c(colMeans(metric_bs_level_bu)[1], mean(colMeans(metric_bs_level_bu)[2:8]))
bss[1, 2:1] <- c(colMeans(metric_bs_level_base)[1], mean(colMeans(metric_bs_level_base)[2:8]))


