library(DiscreteRecon)

data <- readRDS('experiment/dc_crime/data/basef.rds')
data <- Filter(function(x){!("error" %in% class(x))}, data)
SMAT <- rbind(rep(1, 4), diag(4))

tmpf <- function(input){
  if (max(input$y) >= 2){
    DOMAIN <- rbind(0, rep(2, 4))
    max_domain <- 2
  } else {
    DOMAIN <- rbind(0, rep(1, 4))
    max_domain <- 1
  }
  T_ <- NROW(input$y)
  ht <- dhier(SMAT, DOMAIN)
  
  basef <- lapply(1:5, function(x){
    max_d <- ifelse(x == 1, max_domain * 4, max_domain)
    bf <- lapply(input$fcasts[,x], function(lambda) {
      if (lambda < 0) { lambda <- 0.00001 }
      c(dpois(0:(max_d-1), lambda), 1 - ppois(max_d-1, lambda))
    })
    do.call(rbind, bf)
  })
  
  train_basef <- lapply(basef, function(x) { 
    output <- x[1:(T_ - 16),]
    colnames(output) <- 0:(NCOL(output)-1)
    output
  })

  test_basef <- lapply(basef, function(x) { 
    output <- x[(T_ - 15):T_,]
    colnames(output) <- 0:(NCOL(output)-1)
    output
  })

  input$y[input$y > max_domain] <- max_domain 
  train_y <- input$y[1:(T_ - 16),]
  test_y <- input$y[(T_-15):T_,]
  

  recf_lst <- list()
  for (method in c("bu", "td", "dfr", "emp")) {
    mdl <- dfr(ht, method, train_y, train_basef)
    recf_lst[[method]] <- reconcile(mdl, test_basef)
  }
  
  recf_lst
}

library(foreach)
cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

res <- foreach(d = iterators::iter(data), .errorhandling = "pass",
               .packages = c("DiscreteRecon")) %dopar%
  {tmpf(d)}

saveRDS(res, "experiment/dc_crime/data/ct_results.rds")
