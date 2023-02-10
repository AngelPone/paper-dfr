library(dplyr)
library(tidyr)
library(ggplot2)
library(tscount)
library(foreach)
library(DiscreteRecon)

death <- read.table("data/Deaths_1x1.txt", header = TRUE) %>%
  mutate(Age = ifelse(Age == "110+", "110", Age)) %>%
  mutate(Total = as.numeric(Total), Age = as.numeric(Age)) %>%
  filter(Age >= 55)
  
exposure <- read.table("data/Exposures.txt", header=TRUE) %>%
  mutate(Age = ifelse(Age == "110+", "110", Age)) %>%
  mutate(Total = as.numeric(Total), Age = as.numeric(Age)) %>%
  filter(Age >= 55)


death_rate <- left_join(death, exposure, by = c("Age", "Year")) %>% 
  rename(death = Total.x, exposure = Total.y) %>%
  select(Age, Year, death, exposure) %>%
  mutate(age_group = ifelse(Age < 65, "55-64",
                     ifelse(Age < 75, "65-74",
                     ifelse(Age < 85, "75-84", "85+")))) %>%
  group_by(Year, age_group) %>%
  summarise(death = sum(death), exposure = sum(exposure), .groups = "drop") %>%
  mutate(death_rate = death/exposure)
  

dt <- death_rate %>% group_by(age_group) %>%
  arrange(Year) %>%
  mutate(roll_11 = lag(death_rate, 11), roll_1 = lag(death_rate, 1), threshold = (roll_1 - (roll_11 - roll_11) / 10)) %>%
  mutate(y = ifelse(death_rate > threshold * 1.005, 1, 0)) %>%
  drop_na() %>%
  select(Year, age_group, y)


dt %>% filter(age_group == "65-74") %>%
  ggplot(aes(x = Year, y = y, color=age_group)) + 
  geom_point()

# a <- a %>% group_by(Age) %>% arrange(Year) %>%
#   mutate(roll_11 = lag(Total, 11), roll_1 = lag(Total, 1), expected = roll_1 - (roll_11 - roll_1)/10) %>%
#   mutate(y = ifelse(Total < expected, 1, 0)) %>%
#   select(Year, Age, Total, y) %>%
#   drop_na()
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)


build_bm <- function(x){
  y <- x$y
  train_data <- lapply(21:length(y), function(year){
    list(train=y[1:(year-1)], test=y[year])
  })
  
  train_model <- foreach(dt = iterators::iter(train_data), .errorhandling="pass", 
                         .packages = c("tscount")) %dopar% 
    {
      mdl <- tsglm(dt$train, model=list(past_obs=1:3, past_mean=5),
                   link="log", distr="poisson")
      pred <- predict(mdl, n.ahead=1)$pred
      list(distr=c(dpois(0, pred), 1-dpois(0, pred)), y=dt$test)
    }
  train_model
}

build_bm_total <- function(x){
  y <- x$y
  train_data <- lapply(21:length(y), function(year){
    list(train=y[1:(year-1)], test=y[year])
  })
  
  foreach(dt = iterators::iter(train_data), .errorhandling="pass", 
                         .packages = c("tscount")) %dopar% 
    {
      mdl <- tsglm(dt$train, model=list(past_obs=1:3, past_mean=5),
                   link="log", distr="poisson")
      pred <- predict(mdl, n.ahead=1)$pred
      list(distr=c(dpois(0:3, pred), 1-sum(dpois(0:3, pred))), y=dt$test)
    }
}

prepare <- function(b, total){
  func <- function(y){
    bf <- do.call(rbind, lapply(y, function(g){g$distr}))
    colnames(bf) <- 0:(dim(bf)[2]-1)
    bf
  }
  basef_b <- lapply(b$bm, func)
  basef_t <- func(total)
  true_b <- do.call(cbind, lapply(b$bm, function(g){sapply(g, function(z){z$y})}))
  list(basef = c(list(basef_t), basef_b), y = true_b)
}

b <- dt %>% group_nest() %>%
  mutate(bm=purrr::map(data, build_bm))

total <- dt %>% ungroup() %>% group_by(Year) %>% summarise(y=sum(y))
total <- build_bm_total(total)



final_dt <- prepare(b, total)
basef <- final_dt$basef
y <- final_dt$y

train_basef <- lapply(basef, function(x){x[1:59,]})
train_y <- y[1:59,]
test_basef <- lapply(basef, function(x){x[60:69,]})
test_y <- y[60:69,]




S <- rbind(rep(1, 4), diag(4))
domain <- rbind(rep(0, 4), rep(1, 4))

train_y <- dhts(train_y, S, domain)
test_y <- dhts(test_y, S, domain)



# base
brier_score(test_basef, test_y)


# bu
pred_bu <- marginal2Joint(test_basef, test_y$meta, "bu")
brier_score(pred_bu, test_y)


# td
td_mdl <- topdown.train(train_y)
pred_td <- reconcile(td_mdl, test_basef)
brier_score(pred_td, test_y)
# dfr
recon_mdl <- reconcile_train(train_basef, train_y, step_wise = FALSE)
pred_dfr <- reconcile(recon_mdl, test_basef, meta = train_y$meta)
brier_score(pred_dfr, test_y)


