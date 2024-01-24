library(tsintermittent)
library(foreach)
crostonNB <- function(x, h=1) {
  x <- x[min(which(x>0)):length(x)]
  crost_mdl <- crost(x, h=h, type = "sba")
  crost_var <- var(x - crost_mdl$frc.in, na.rm = TRUE)
  size <- mean((crost_mdl$frc.in)^2 / (crost_var - crost_mdl$frc.in), 
            na.rm = TRUE)
  size <- ifelse(size < 0, crost_mdl$frc.out, size)
  list(size = size, point = crost_mdl$frc.out)
}


process <- function(hierarchy) {
  bts <- hierarchy$series
  total <- rowSums(bts)
  total <- total[min(which(total > 0)):length(total)]
  bts <- lapply(iterators::iter(bts, by = "column"), function(x) {
    x[min(which(x>0)):length(x)]
  })
  
  min_size <- min(sapply(bts, length))
  if (min_size < 730 + 300 + 10) {
    return (NULL)
  } else {
    output <- c(list(total), bts)
    names(output) <- c("Total", hierarchy$names)
    return(output)
  }
}


rolling_train <- function(x, window_length = 730, h=1) {
  l <- length(x)
  x_lst <- lapply((window_length-1):0, function(o){
    list(train = x[1:(l-o-h)], test = x[(l-o-h+1):(l-o)])
  })
  fcasts <- foreach(x = iterators::iter(x_lst), .export = "crostonNB", .packages = "tsintermittent") %dopar% {
    list(fcasts = crostonNB(x$train, h=h), test = x$test)
  }
  fcasts_ <- do.call(rbind, lapply(fcasts, function(x){x$fcasts}))
  test <- do.call(rbind, lapply(fcasts, function(x){x$test}))
  list(fcasts = fcasts_, test=test)
}
library(DiscreteRecon)

dist2prob <- function(dist, max) {
  point <- simplify2array(dist[,"point"])
  size <- simplify2array(dist[,"size"])
  output <- lapply(1:length(point), function(x){
    if (point[x] == size[x]) {
      return(c(dpois(0:(max-1), point[x]), 1 - ppois(max-1, point[x])))
    } else {
      return(c(dnbinom(0:(max-1), size=size[x], mu=point[x]), 
               1-pnbinom(max-1, size=size[x], mu=point[x])))
    }
  })
  output <- do.call(rbind, output)
  colnames(output) <- 0:max
  output
}


recon_f <- function(hist, fcasts) {
  m <- length(fcasts) - 1
  train_size <- NROW(fcasts[[1]]$fcasts) - 28
  s_mat <- rbind(rep(1, m), diag(m))
  domain <- rbind(rep(0, m), sapply(hist[2:(1+m)], function(x) { max(x) }))
  ht <- dhier(s_mat, domain)
  hist <- do.call(cbind, lapply(fcasts[2:(1+m)], function(x) {x$test}))
  maxs <- c(sum(domain[2,]), domain[2,])
  train_fcasts <- lapply(seq_along(fcasts), function(x) {
    dist2prob(fcasts[[x]]$fcasts[1:train_size, ], maxs[x])
  })
  test_fcasts <- fcasts <- lapply(seq_along(fcasts), function(x) {
    dist2prob(fcasts[[x]]$fcasts[(train_size+1):(train_size+28),], maxs[x])
  })
  
  mdl <- dfr(ht, method = "sdfr",obs_train = hist[1:train_size, ], bf_train = train_fcasts)
  reconcile(mdl, test_fcasts)
}


evaluate <- function(hist, reconf, fcasts, window_length=730, h=1) {
  m <- length(hist) - 1
  s_mat <- rbind(rep(1, m), diag(m))
  train_size <- NROW(fcasts[[1]]$fcasts) - 28
  test_size <- NROW(reconf)
  domain <- rbind(rep(0, m), sapply(hist[2:(1+m)], function(x) { max(x) }))
  ht <- dhier(s_mat, domain)
  obs <- do.call(cbind, lapply(hist[2:(1+m)], function(x) {x[(length(x) - (train_size + test_size) + 1):length(x)]}))
  bs <- brier_score(reconf, obs[(train_size+1):(train_size+test_size),], ht)
  
  maxs <- c(sum(domain[2,]), domain[2,])
  test_fcasts <- fcasts <- lapply(seq_along(fcasts), function(x) {
    dist2prob(fcasts[[x]]$fcasts[(train_size+1):(train_size+test_size),], maxs[x])
  })
  
  mdl_bu <- dfr(ht, method = "bu")
  reconf_bu <- reconcile(mdl_bu, test_fcasts)
  bs_bu <- brier_score(reconf_bu, obs[(train_size+1):(train_size+test_size),], ht)
  
  mdl_td <- dfr(ht, method = "td", obs_train = obs[1:train_size,])
  reconf_td <- reconcile(mdl_td, test_fcasts)
  bs_td <- brier_score(reconf_td, obs[(train_size+1):(train_size+test_size),], ht)
  
  mdl_emp <- dfr(ht, method = "emp", obs_train = obs[1:train_size,])
  reconf_emp <- reconcile(mdl_emp, h = test_size)
  bs_emp <- brier_score(reconf_emp, obs[(train_size+1):(train_size+test_size),], ht)
  
  bs_base <- brier_score(test_fcasts, obs[(train_size+1):(train_size+test_size),], ht)
  
  list(bu = bs_bu, td = bs_td, base=bs_base, sdfr=bs, emp=bs_emp)
}

summary <- function(bs) {
  output <- lapply(c("base", "bu", "td", "sdfr", "emp"), function(method) {
    output <- do.call(rbind, lapply(iterators::iter(bs), function(idx) {
      total <- idx[[method]]$series[1]
      bottom <- mean(idx[[method]]$series[2:length(idx[[method]]$series)])
      hierarchy <- sum(idx[[method]]$hierarchy)
      c(total, bottom, hierarchy)
    })
    )
    colnames(output) <- c("Total", "Bottom", "Hierarchy")
    colMeans(output)
  })
  
  output <- do.call(cbind, output)
  colnames(output) <- c("Base", "Bottom-up", "Top-down", "SDFR", "Empirical")
  output
}
library(DiscreteRecon)
emp_dist <- function(obs, window_length=730) {
  m <- length(obs) - 1
  test_size <- 28
  s_mat <- rbind(rep(1, m), diag(m))
  domain <- rbind(rep(0, m), sapply(hist[2:(1+m)], function(x) { max(x) }))
  ht <- dhier(s_mat, domain)
  hist <- do.call(cbind, lapply(fcasts[2:(1+m)], function(x) {x$test}))
  
  mdl <- dfr(ht, method = "emp", obs_train = hist[1:(window_length - test_size+1)])
  
}





