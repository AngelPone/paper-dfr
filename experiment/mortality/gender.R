library(dplyr)
library(tidyr)
library(ggplot2)
library(tscount)

a <- read.table("data/au_51.txt", header=TRUE) %>%
  drop_na() %>%
  filter(Age %in% c("20-24")) %>%
  mutate(Year = as.numeric(Year), Male=as.numeric(Male)*1000, Female=as.numeric(Female)*1000) %>%
  select(Year, Female, Male)

a <- a %>% arrange(Year) %>%
  mutate(roll_11 = lag(Female, 11), roll_1 = lag(Female, 1), expected_f = roll_1 - (roll_11 - roll_1)/10) %>%
  mutate(roll_11 = lag(Male, 11), roll_1 = lag(Male, 1), expected_m = roll_1 - (roll_11 - roll_1)/10) %>%
  mutate(yf = ifelse(Female < expected_f, 1, 0)) %>%
  mutate(ym = ifelse(Male < expected_m, 1, 0)) %>%
  select(Year, yf, ym) %>%
  drop_na()

build_bm <- function(x, y){
  y <- x[[y]]
  train_data <- lapply(21:length(y), function(year){
    list(train=y[1:(year-1)], test=y[year])
  })
  
  train_model <- lapply(train_data, function(dt){
    mdl <- tsglm(dt$train, model=list(past_obs=1:3, past_mean=5),
                 link="log", distr="poisson")
    pred <- predict(mdl, n.ahead=1)$pred
    list(distr=c(dpois(0, pred), 1-dpois(0, pred)), y=dt$test)
  })
  
  train_model
}

build_bm_total <- function(x){
  y <- x$yf + x$ym
  train_data <- lapply(21:length(y), function(year){
    list(train=y[1:(year-1)], test=y[year])
  })
  
  train_model <- lapply(train_data, function(dt){
    mdl <- tsglm(dt$train, model=list(past_obs=1:3, past_mean=5),
                 link="log", distr="poisson")
    pred <- predict(mdl, n.ahead=1)$pred
    list(distr=c(dpois(0:1, pred), 1-sum(dpois(0:1, pred))), y=dt$test)
  })
  
  train_model
}

f <- build_bm(a, "yf")
m <- build_bm(a, "ym")

total <- build_bm_total(a)



li2li <- function(x){
  basef <- do.call(rbind, lapply(x, function(g){g$distr}))
  colnames(basef) <- 0:(dim(basef)[2]-1)
  basef
}


y <- cbind(a$yf, a$ym)
basef <- list(li2li(total), li2li(f), li2li(m))

train_basef <- lapply(basef, function(x){x[1:59,]})
train_y <- y[1:59,]
test_basef <- lapply(basef, function(x){x[60:69,]})
test_y <- y[60:69,]



library(DiscreteRecon)
S <- rbind(rep(1, 2), diag(2))
domain <- rbind(rep(0, 2), rep(1, 2))

train_y <- dhts(train_y, S, domain)
test_y <- dhts(test_y, S, domain)

recon_mdl <- reconcile_train(train_basef, train_y, step_wise = FALSE)

pred_dfr <- reconcile(recon_mdl, test_basef, meta = train_y$meta)

# dfr
brier_score(pred_dfr, test_y)

# base
brier_score(test_basef, test_y)


# bu
pred_bu <- marginal2Joint(test_basef, test_y$meta, "bu")
brier_score(pred_bu, test_y)


# td
td_mdl <- topdown.train(train_y)
pred_td <- reconcile(td_mdl, test_basef)
brier_score(pred_td, test_y)
