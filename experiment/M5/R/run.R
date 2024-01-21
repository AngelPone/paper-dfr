# basef
setwd("experiment/M5")
source("R/basef.R")

dtarget <- readRDS("data.rds")
library(dplyr)
library(foreach)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

dtarget <- dtarget %>% rowwise() %>%
  mutate(hierarchy = list(process(hierarchy))) %>%
  filter(!is.null(hierarchy))

dtarget$fcasts <- vector("list", NROW(dtarget))
for (hierarchy_idx in seq_along(dtarget$hierarchy)) {
  print(timestamp())
  print(sprintf("base forecasts of %s in State %s", dtarget$item_id[[hierarchy_idx]], dtarget$state_id[[hierarchy_idx]]))
  ht_fcasts <- list()
  hierarchy <- dtarget$hierarchy[[hierarchy_idx]]
  for (series in names(hierarchy)){
    ht_fcasts[[series]] <- rolling_train(hierarchy[[series]], window_length = 730, h = 1)
  }
  dtarget$fcasts[[hierarchy_idx]] <- ht_fcasts
  saveRDS(dtarget, "bf.rds")
}
library(dplyr)
dtarget <- readRDS("bf.rds") %>% filter(salesmax <= 3)
dtarget$recf <- vector("list", NROW(dtarget))
for (hierarchy_idx in seq_along(dtarget$hierarchy)) {
  print(hierarchy_idx)
  dtarget$recf[[hierarchy_idx]] <- recon_f(dtarget$hierarchy[[hierarchy_idx]], dtarget$fcasts[[hierarchy_idx]])
}
saveRDS(dtarget, "reconf-3.rds")
