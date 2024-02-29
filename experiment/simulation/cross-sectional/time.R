setwd("experiment/simulation")
dts <- readRDS("results/cs-res.rds")
library(DiscreteRecon)

history_length = 150
window_length = 300
test_length = 30
domain <- rbind(c(0, 0), c(1, 1))
s_mat <- rbind(c(1, 1), diag(2))
n = history_length+window_length+test_length

total_time <- 0
for (i in 1:100) {
  dt <- dts[[i]]
  obs <- dt$series[n - ((window_length + test_length):1) + 1,]
  
  obs_train <- obs[1:window_length,]
  obs_test <- obs[(window_length+1):(window_length + test_length),]
  
  bf_train <- dt$fcasts$base_train
  bf_test <- dt$fcasts$base_test
  
  ht <- dhier(s_mat, domain)
  
  start <- Sys.time()
  mdl <- dfr(ht, "dfr", bf_train = bf_train, obs_train = obs_train)
  end <- Sys.time()
  total_time <- total_time + end - start
}

output_path <- "../../manuscript/figures/time.csv"
results <- c(r = NROW(mdl$A),
             q = NCOL(mdl$A),
             n = sum(mdl$A != 0) - sum(mdl$A == 1),
             time = total_time / 100)
if (file.exists(output_path)) {
  file <- read.csv(output_path)
  file$simulation_cs <- results
  write.csv(file, output_path, row.names = FALSE)
} else {
  write.csv(data.frame(simulation_cs = results), 
            output_path,
            row.names = FALSE)
  
}

