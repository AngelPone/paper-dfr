library(DiscreteRecon)

data <- readRDS('data/data.rds')

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
    lapply(input$fcasts[,x], function(lambda) {
      c(dpois(0:(max_d-1), lambda), 1 - ppois(max_d-1, lambda))
    })
  })
  
  train_basef <- lapply(basef, function(x) { x[1:(T_ - 16),]})
  test_basef <- lapply(basef, function(x) { x[(T_ - 15):T_,]})
  
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
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

res <- foreach(d = iterators::iter(ct), .errorhandling = "pass",
               .packages = c("DiscreteRecon")) %dopar%
  {tmpf(d)}

saveRDS(res, "data/ct_results.rds")


