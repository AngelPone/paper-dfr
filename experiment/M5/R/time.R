setwd("experiments/M5")

source("R/basef.R")
dts <- readRDS("results/bf.rds") %>% filter(salesmax == 4)

all4_dts <- which(sapply(dts$hierarchy, function(x) { sum(sapply(x[2:length(x)], max) == 4) == 4 }))
dts <- dts[all4_dts,]


fcasts <- dts$fcasts[[1]]
hierarchy <- dts$hierarchy[[1]]
m <- length(fcasts) - 1
train_size <- NROW(fcasts[[1]]$fcasts) - 28
s_mat <- rbind(rep(1, m), diag(m))
domain <- rbind(rep(0, m), rep(4, m))
ht <- dhier(s_mat, domain)
maxs <- c(16, rep(4, m))
hist <- do.call(cbind, lapply(fcasts[2:(1 + m)], function(x) {
  x$test
}))
train_fcasts <- lapply(seq_along(fcasts), function(x) {
  dist2prob(fcasts[[x]]$fcasts[1:train_size, ], maxs[x])
})
test_fcasts <- fcasts <- lapply(seq_along(fcasts), function(x) {
  dist2prob(fcasts[[x]]$fcasts[(train_size + 1):(train_size + 28), ], maxs[x])
})

total_time <- 0
for (i in 1:50) {
  print(i)
  start <- Sys.time()
  mdl <- dfr(ht, method = "sdfr", obs_train = hist[1:train_size, ], bf_train = train_fcasts)
  end <- Sys.time()
  total_time <- total_time + end - start
}

total_n <- 0
for (i in seq_along(mdl$A)) {
  A <- mdl$A[[i]]$model$A
  total_n <- total_n + sum(A != 0) - sum(A == 1)
}
total_time <- as.double(total_time, units = "secs")

results <- c(r = NROW(ht$coherent_domain),
             q = NROW(ht$incoherent_domain),
             n = total_n,
             time = total_time)

if (file.exists(output_path)) {
  file <- read.csv(output_path)
  file$M5 <- results
  write.csv(file, output_path, row.names = FALSE)
} else {
  write.csv(data.frame(M5 = results), 
            output_path,
            row.names = FALSE)
}


