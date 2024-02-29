setwd("experiment/dc_crime")
library(DiscreteRecon)

input <- readRDS('data/basef.rds')[[1]]
SMAT <- rbind(rep(1, 4), diag(4))
DOMAIN <- rbind(0, rep(2, 4))
T_ <- NROW(input$y)
ht <- dhier(SMAT, DOMAIN)
max_domain <- 2
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


input$y[input$y > max_domain] <- max_domain 
train_y <- input$y[1:(T_ - 16),]

total_time <- 0
for (i in 1:50) {
  print(i)
  start <- Sys.time()
  mdl <- dfr(ht, "dfr", train_y, train_basef)
  end <- Sys.time()
  total_time <- end - start + total_time
}

total_time <- as.double(total_time, units = "secs")
output_path <- "../../manuscript/figures/time.csv"
results <- c(r = NROW(mdl$A),
             q = NCOL(mdl$A),
             n = sum(mdl$A != 0) - sum(mdl$A == 1),
             time = total_time / 50) 
if (file.exists(output_path)) {
  file <- read.csv(output_path)
  file$crime <- results
  write.csv(file, output_path, row.names = FALSE)
} else {
  write.csv(data.frame(crime = results), 
            output_path,
            row.names = FALSE)
  
}

