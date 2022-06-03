library(gratis)
library(feasts)
library(tsibble)

source("config.R")

# simulation setting
time_window = 50 + 76 + 3
nsim = 100
lb_window = 50 * 7
n_samples = time_window * 7 - lb_window - 6
train_size =  n_samples - 21
test_size = 21

generate_bts <- function(x){
  d <- rep(0, 7)
  d[order(rbeta(7, 4, 1:7))[1:x]] <- 1
  d
}

simulate_temporal2 <- function(time_window, times){
  total <- arima_model(4, p=3, d=0, D=0) %>% 
    generate(time_window, times) %>%
    group_by_key() %>%
    mutate(value = (value - min(value))/(max(value) - min(value)) * 2.5 +  2) %>%
    ungroup() %>%
    mutate(value = rpois(nrow(.), value)) %>%
    mutate(obs=ifelse(value > 7, 7, value))
  
  tmp <- function(x){
    matrix(sapply(as.list(x), generate_bts, simplify = TRUE), nrow = 1)[1,]
  }
  bts <- group_by_key(total) %>%
    group_map(~ tmp(.x$obs))
}

series <- simulate_temporal2(time_window, nsim)

basef_temporal2 <- list()
for (i in 73:nsim){
  basef_temporal2[[i]] <- train_test_prepare(reduction2lst(series[[i]]))
}
saveRDS(basef_temporal2, "basef_temporal_2.rds")


