setwd("experiments/simulation")
source('temporal/config.R')
library(DiscreteRecon)

basef <- readRDS('results/basef.rds')
ht <- dhier(smatrix, domain)

total_time <- 0
for (i in 1:50) {
  print(i)
  bf <- basef[[i]]
  bf_train <- bf$x_train
  obs_train <- bf$y_train
  start <- Sys.time()
  mdl <- dfr(ht, method = "sdfr", obs_train = obs_train, bf_train = bf_train)
  end <- Sys.time()
  total_time <- total_time + end - start
}


output_path <- "../../manuscript/figures/time.csv"
total_n <- 0
for (i in seq_along(mdl$A)) {
  A <- mdl$A[[i]]$model$A
  total_n <- total_n + sum(A != 0) - sum(A == 1)
}
total_time <- as.double(total_time, units = "secs")

results <- c(r = NROW(ht$coherent_domain),
             q = NROW(ht$incoherent_domain),
             n = total_n,
             time = total_time / 50)

if (file.exists(output_path)) {
  file <- read.csv(output_path)
  file$simulation_te <- results
  write.csv(file, output_path, row.names = FALSE)
} else {
  write.csv(data.frame(simulation_te = results), 
            output_path,
            row.names = FALSE)
}
