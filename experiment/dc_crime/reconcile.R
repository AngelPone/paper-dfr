library(DiscreteRecon)

train_basef <- do.call(c, lapply(paste0('_rslurm_dfr_base_train/', list.files('_rslurm_dfr_base_train', "results")),
                                 readRDS))
test_basef <- readRDS(paste0("_rslurm_dfr_base_test/", list.files("_rslurm_dfr_base_test", "results")))

max_domain <- 2

DOMAIN <- rbind(rep(0, 4), rep(max_domain, 4))
SMAT <- rbind(rep(1, 4), diag(4))

dt <- list(train=train_basef, test=test_basef)
# dt <- readRDS("../experiment/dc_crime/data/basef.rds")

dt$train <- Filter(function(x){
  if (!("try-error" %in% class(x))) return(all(x$mean > 0))
  return(FALSE)
}, dt$train)
dt$test <- Filter(function(x){
  if (!("try-error" %in% class(x))) return(all(x$mean > 0))
  return(FALSE)
}, dt$test)


# dt$train <- Filter(function(x){max(x$train) == max_domain}, dt$train)
# dt$test <- Filter(function(x){max(x$train) == max_domain}, dt$test)

train_basef <- list()
train_basef[[1]] <- t(sapply(dt$train, function(x){
  p <- dpois(0:rowSums(DOMAIN)[2], x$mean[1])
  p[length(p)] <- 1 - sum(p[1:(length(p) - 1)])
  names(p) <- 0:rowSums(DOMAIN)[2]
  p
}, simplify = "array"))

train_basef[2:5] <- lapply(2:5, function(x){
  t(sapply(dt$train, function(y){
    p <- dpois(0:(DOMAIN[2, 1]-1), y$mean[x])
    p <- c(p, 1-sum(p))
    names(p) <- 0:DOMAIN[2, 1]
    p
  }, simplify = "array"))
})

test_basef <- list()
test_basef[[1]] <- t(sapply(dt$test, function(x){
  p <- dpois(0:rowSums(DOMAIN)[2], x$mean[1])
  p[length(p)] <- 1 - sum(p[1:(length(p) - 1)])
  names(p) <- 0:rowSums(DOMAIN)[2]
  p
}, simplify = "array"))

test_basef[2:5] <- lapply(2:5, function(x){
  t(sapply(dt$test, function(y){
    p <- dpois(0:(DOMAIN[2, 1]-1), y$mean[x])
    p <- c(p, 1-sum(p))
    names(p) <- 0:DOMAIN[2, 1]
    p
  }, simplify = "array"))
})




train_real <- t(sapply(dt$train, function(x){as.numeric(x$test)}, simplify = "array"))
train_real[train_real >= max_domain] = max_domain
train_real <- dhts(train_real, SMAT, DOMAIN)

test_real <- t(sapply(dt$test, function(x){as.numeric(x$test)}, simplify = "array"))
test_real[test_real >=max_domain] = max_domain
test_real <- dhts(test_real, SMAT, DOMAIN)

rec_mdl <- reconcile_train(train_basef, train_real, 
                           step_wise = FALSE, optimized = TRUE)
td_mdl <- topdown.train(train_real)


fcasts <- reconcile(rec_mdl, test_basef, test_real$meta)
fcasts_td <- reconcile(td_mdl, test_basef)

series_bs <- function(x, f, y){
  if (is.list(f)){
    new_f <- lapply(f, function(y){y[x,,drop=FALSE]})
  }else {
    new_f <- f[x,,drop=FALSE]
    class(new_f) <- class(f)
  }
  
  y$bts <- y$bts[x,,drop=FALSE]
  brier_score(new_f, y)
}


dfr <- lapply(1:NROW(fcasts), series_bs, fcasts, test_real)
td <- lapply(1:NROW(fcasts), series_bs, fcasts_td, test_real)
base <- lapply(1:NROW(fcasts), series_bs, test_basef, test_real)
bu <- lapply(1:NROW(fcasts), series_bs, marginal2Joint(test_basef, test_real$meta, method="bu"),
             test_real)
base <- lapply(base, function(x){
  x$hierarchy <- sum(x$hierarchy)
  x
})



b_res <- sapply(list(dfr, td, base, bu), function(x){as.vector(sapply(x, function(y){y$series[2:5]}))})
colnames(b_res) <- c("DFR", "Top-Down", "Base", "Bottom-Up")

library(tsutils)
b_res[,"Bottom-Up"] = b_res[,"Base"]

pdf(width = 15, height = 5, pointsize = 16)
par(mfrow = c(1,3))
nemenyi(b_res, plottype="vmcb", 
        labels=c("DFR", "Top-Down", "Base", "Bottom-Up"),
        main = sprintf("MCB Test on %s of %s", "Brier Score", "bottom-level series"))

t_res <- sapply(list(dfr, td, base, bu), function(x){as.vector(sapply(x, function(y){y$series[1]}))})
colnames(t_res) <- c("DFR", "Top-Down", "Base", "Bottom-Up")
t_res[,"Top-Down"] = t_res[,"Base"]
nemenyi(t_res, plottype="vmcb", 
        labels=c("DFR", "Top-Down", "Base", "Bottom-Up"),
        main = sprintf("MCB Test on %s of %s", "Brier Score", "total series"))
colMeans(t_res)

hierarchy_res <- sapply(list(dfr, td, base, bu), function(x){as.vector(sapply(x, function(y){y$hierarchy}))})
colnames(hierarchy_res) <- c("DFR", "Top-Down", "Base", "Bottom-Up")
nemenyi(hierarchy_res, plottype="vmcb", 
        labels=c("DFR", "Top-Down", "Base", "Bottom-Up"),
        main = sprintf("MCB Test on %s of %s", "Brier Score", "hierarchy"))
colMeans(hierarchy_res)

dev.off()

brier_score(fcasts, test_real)
brier_score(fcasts_td, test_real)
brier_score(test_basef, test_real)
brier_score(marginal2Joint(test_basef, test_real$meta, "bu"), test_real)
