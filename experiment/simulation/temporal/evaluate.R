
library(DiscreteRecon)
res <- readRDS("SFP.rds")
tdres <- readRDS('td.rds')
basef <- readRDS('basef.rds')
step_models <- readRDS('step_models.rds')

idx <- which(sapply(step_models, length) == 6)
basef <- basef[idx]

output_mae <- matrix(0, 8, 4)
lst_mae <- list()
mae <- function(x, y){mean(abs(x-y))}

output_rmse <- matrix(0, 8, 4)
lst_rmse <- list()
rmse <- function(x, y){sqrt(mean((x-y)^2))}

lst_bs <- list()

### base forecast ######
lst_mae$base <- t(sapply(seq_along(basef), 
                         function(x){ point_metric(basef[[x]]$x_test, res[[x]]$y, mae)},
                         simplify="array"))

lst_rmse$base <- t(sapply(seq_along(basef), 
                          function(x){ point_metric(basef[[x]]$x_test, res[[x]]$y, rmse)},
                          simplify="array"))
lst_bs$base <- lapply(seq_along(basef), function(x){brier_score(basef[[x]]$x_test, res[[x]]$y)})


### bottom up #######
lst_mae$bu <- t(sapply(seq_along(basef), 
                         function(x){ point_metric(marginal2Joint(basef[[x]]$x_test, res[[x]]$y$meta, method = "bu"), res[[x]]$y, mae)},
                         simplify="array"))

lst_rmse$bu <- t(sapply(seq_along(basef), 
                         function(x){ point_metric(marginal2Joint(basef[[x]]$x_test, res[[x]]$y$meta, method = "bu"), res[[x]]$y, rmse)},
                         simplify="array"))
lst_bs$bu <- lapply(seq_along(basef), function(x){brier_score(marginal2Joint(basef[[x]]$x_test, res[[x]]$y$meta, method = "bu"), res[[x]]$y)})


### top down ####
lst_mae$td <- t(sapply(tdres, function(x){point_metric(x$dist, x$y, mae)},simplify="array"))
lst_rmse$td <- t(sapply(tdres, function(x){point_metric(x$dist, x$y, rmse)},simplify="array"))

lst_bs$td <- lapply(tdres, function(x){brier_score(x$dist, x$y)})

### rec ######
lst_mae$sfr <- t(sapply(res, function(x){point_metric(x$dist, x$y, mae)},simplify="array"))
lst_rmse$sfr <- t(sapply(res, function(x){point_metric(x$dist, x$y, rmse)},simplify="array"))

lst_bs$sfr <- lapply(res, function(x){brier_score(x$dist, x$y)})


### summarize ######
ns <- c("base", "bu", "td", "sfr")
output_mae <- sapply(ns, function(x){colMeans(lst_mae[[x]])})
output_rmse <- sapply(ns, function(x){colMeans(lst_rmse[[x]])})
output_bs_series <- sapply(ns, function(x){rowMeans(sapply(lst_bs[[x]], function(x){x$series}))})
output_bs_hierarchy <- sapply(ns, function(x){mean(sapply(lst_bs[[x]], function(x){sum(x$hierarchy)}))})


### MCB Test #####
library(tsutils)
library(dplyr)
mcb_plot <- function(dat, metric, level){
  if (level == "total") dat[,"td"] = dat[, "base"]
  if (level == "bottom") dat[,"bu"] = dat[, "base"]
  level_name <- level
  if (level != "hierarchy") level <- paste0(level, " level")
  metric_name <- metric
  if (metric == "Brier Score") metric_name <- "BS"
  pdf(sprintf("figures/temporal_mcb_%s_%s.pdf", metric_name, level_name),
      width = 5, height = 4, pointsize = 10)
  nemenyi(dat, plottype="vmcb", 
          labels=c("Base", "Bottom-Up", "Top-Down", "SDFR"),
          main = sprintf("MCB Test on %s of %s", metric, level))
  dev.off()
}


sapply(ns, function(x){lst_mae[[x]][,1]}, simplify = "array") %>%
  mcb_plot("MAE", "total")
sapply(ns, function(x){as.vector(lst_mae[[x]][,2:8])}, simplify = "array") %>%
  mcb_plot("MAE", "bottom")

sapply(ns, function(x){lst_rmse[[x]][,1]}, simplify = "array") %>%
  mcb_plot("RMSE", "total")
sapply(ns, function(x){as.vector(lst_rmse[[x]][,2:8])}, simplify = "array") %>%
  mcb_plot("RMSE", "bottom")


sapply(ns, function(x){sapply(lst_bs[[x]], function(x){sum(x$hierarchy)})}, simplify = "array") %>%
  mcb_plot("Brier Score", "hierarchy")
sapply(ns, function(x){sapply(lst_bs[[x]], function(x){x$series[1]})}, simplify = "array") %>%
  mcb_plot("Brier Score", "total")
sapply(ns, function(x){as.vector(sapply(lst_bs[[x]], function(x){x$series[2:8]}))}, simplify = "array") %>%
  mcb_plot("Brier Score", "bottom")


