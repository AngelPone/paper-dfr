library(DiscreteRecon)

idx <- readRDS('data/idx.rds')
data <- readRDS('data/data.rds')

SMAT <- rbind(rep(1, 4), diag(4))

train_basef <- do.call(c, lapply(paste0('_rslurm_dfr_base_train/', list.files('_rslurm_dfr_base_train', "results")),
                                 readRDS))
for (i in seq_along(train_basef)){
  if ("try-error" %in% class(train_basef[[i]])) next
  train_basef[[i]]$id <- idx$train[[i]]
}
test_basef <- readRDS(paste0("_rslurm_dfr_base_test/", list.files("_rslurm_dfr_base_test", "results")))
for (i in seq_along(test_basef)){
  if ("try-error" %in% class(test_basef[[i]])) next
  test_basef[[i]]$id <- idx$test[[i]]
}

# max_domain <- 2

dt <- list(train = train_basef, test = test_basef)
# dt <- readRDS("../experiment/dc_crime/data/basef.rds")

dt$train <- Filter(function(x){
  if (!("try-error" %in% class(x))) return(all(x$mean > 0))
  return(FALSE)
}, dt$train)
dt$test <- Filter(function(x){
  if (!("try-error" %in% class(x))) return(all(x$mean > 0))
  return(FALSE)
}, dt$test)

train_idx <- unname(sapply(dt$train, function(x){x$id}))
test_idx <- unname(sapply(dt$test, function(x){x$id}))

ct <- list()
j <- 1
for (i in unique(train_idx)){
  ct[[j]] <- list(train=dt$train[which(train_idx == i)],
                  test=dt$test[which(test_idx == i)])
  j = j + 1
}

for (i in seq_along(data)){
  ct[[i]]$data <- data[[i]]
}

rm(train_basef, test_basef, dt, idx, train_idx, test_idx, data)

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



tmpf <- function(input){
  if (max(input$data) >= 2){
    DOMAIN <- rbind(0, rep(2, 4))
    max_domain <- 2
  } else {
    DOMAIN <- rbind(0, rep(1, 4))
    max_domain <- 1
  }
  
  train_basef <- list()
  train_basef[[1]] <- t(sapply(input$train, function(y){
    p <- dpois(0:rowSums(DOMAIN)[2], y$mean[1])
    p[length(p)] <- 1 - sum(p[1:(length(p) - 1)])
    names(p) <- 0:rowSums(DOMAIN)[2]
    p
  }, simplify = "array"))
  
  train_basef[2:5] <- lapply(2:5, function(y){
    t(sapply(input$train, function(k){
      p <- dpois(0:(DOMAIN[2, 1]-1), k$mean[y])
      p <- c(p, 1-sum(p))
      names(p) <- 0:DOMAIN[2, 1]
      p
    }, simplify = "array"))
  })
  
  test_basef <- list()
  test_basef[[1]] <- t(sapply(input$test, function(y){
    p <- dpois(0:rowSums(DOMAIN)[2], y$mean[1])
    p[length(p)] <- 1 - sum(p[1:(length(p) - 1)])
    names(p) <- 0:rowSums(DOMAIN)[2]
    p
  }, simplify = "array"))
  
  test_basef[2:5] <- lapply(2:5, function(k){
    t(sapply(input$test, function(y){
      p <- dpois(0:(DOMAIN[2, 1]-1), y$mean[k])
      p <- c(p, 1-sum(p))
      names(p) <- 0:DOMAIN[2, 1]
      p
    }, simplify = "array"))
  })
  
  
  train_real <- t(sapply(input$train, function(y){as.numeric(y$test)}, simplify = "array"))
  train_real[train_real >= max_domain] = max_domain
  train_real <- dhts(train_real, SMAT, DOMAIN)
  
  test_real <- t(sapply(input$test, function(y){as.numeric(y$test)}, simplify = "array"))
  test_real[test_real >= max_domain] = max_domain
  test_real <- dhts(test_real, SMAT, DOMAIN)
  
  rec_mdl <- reconcile_train(train_basef, train_real,
                             step_wise = FALSE, optimized = TRUE)
  td_mdl <- topdown.train(train_real)
  
  
  fcasts <- reconcile(rec_mdl, test_basef, test_real$meta)
  fcasts_td <- reconcile(td_mdl, test_basef)
  fcasts_bu <- marginal2Joint(test_basef, test_real$meta, "bu")
  
  dfr <- lapply(1:NROW(fcasts), series_bs, fcasts, test_real)
  td <- lapply(1:NROW(fcasts), series_bs, fcasts_td, test_real)
  base <- lapply(1:NROW(fcasts), series_bs, test_basef, test_real)
  bu <- lapply(1:NROW(fcasts), series_bs, fcasts_bu, test_real)
  base <- lapply(base, function(x){
    x$hierarchy <- sum(x$hierarchy)
    x
  })
  
  b_res <- sapply(list(dfr, td, base, bu), function(x){as.vector(sapply(x, function(y){y$series[2:5]}))})
  colnames(b_res) <- c("DFR", "Top-Down", "Base", "Bottom-Up")
  b_res[,"Bottom-Up"] = b_res[,"Base"]
  
  t_res <- sapply(list(dfr, td, base, bu), function(x){as.vector(sapply(x, function(y){y$series[1]}))})
  colnames(t_res) <- c("DFR", "Top-Down", "Base", "Bottom-Up")
  t_res[,"Top-Down"] = t_res[,"Base"]
  
  hierarchy_res <- sapply(list(dfr, td, base, bu), function(x){as.vector(sapply(x, function(y){y$hierarchy}))})
  colnames(hierarchy_res) <- c("DFR", "Top-Down", "Base", "Bottom-Up")
  
  list(b_res, t_res, hierarchy_res)
}

library(foreach)
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

res <- foreach(d = iterators::iter(ct), .errorhandling = "pass",
               .packages = c("DiscreteRecon")) %dopar%
  {tmpf(d)}

saveRDS(res, "ct_results.rds")


