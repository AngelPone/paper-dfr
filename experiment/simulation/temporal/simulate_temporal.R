# Simulate Daily Time Series
library(dplyr)
library(foreach)
setwd('simulation/temporal')
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

# config
lb_window = 300
train_size = 500
test_size = 21
h = 7
simulation_size = 1000

time_window = train_size + lb_window + test_size + 6

# random_params <- function(){
#   a <- abs(rnorm(1))
#   a <- c(a, a * 0.7^(1:5))
#   a/sum(a) * 0.9
# }




# simulate function to simulate daily time series
simfunc <- function(time_window, n){
  a = tsintermittent::simID(n, time_window, idi = 2, level = 1.5)
  a[a>1] = 1
  a
}

# function to construct train and test model



# simulate series
series <- simfunc(time_window, simulation_size)

tt <- list()
for (i in 1:100){
  print(i)
  tt[[i]] = train_test_prepare(reduction2lst(series[,i]))
}

saveRDS(tt, 'sim_temporal_base.rds')


