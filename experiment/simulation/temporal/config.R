# file including config of simulation experiments in temporal setting.
library(dplyr)
library(foreach)


time_window = 50 + 76 + 3
nsim = 1000
lb_window = 50 * 7
n_samples = time_window * 7 - lb_window - 6
test_size = 21
train_size =  n_samples - test_size


set.seed(42)

# domain of hierarchy
domain <- rbind(rep(0, 7), rep(1, 7))
smatrix <- rbind(rep(1, 7), diag(7))

## utility functions
logit <- function(x) exp(x)/(1+exp(x))
temporalAgg <- function(x){
  mult <- floor(length(x)/7)*7
  x <- x[(length(x) - mult + 1) : length(x)]
  sapply(split(ifelse(x>=1, 1, 0), ceiling((1:mult)/7)), sum)
}

