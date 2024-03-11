rm(list = ls())
setwd("experiment/M5")
library(dplyr)
dt <- read.csv("/path/to/sales_train_evaluation.csv")

salesmax <- dt %>% select(starts_with("d_")) %>% t() %>%
  apply(2, max)


dtmeta <- dt %>% select(!starts_with("d_")) %>%
  mutate(salesmax = salesmax)



# It is infeasible to use a whole department
dtmeta %>% select(dept_id, item_id) %>% unique() %>%
  group_by(dept_id) %>% count()

# It is feasible to use item sold within 3 or 4 stores located in one state.
dtmeta %>% select(item_id, id, state_id) %>% group_by(item_id, state_id) %>% count()


dtmeta %>% select(item_id, state_id, salesmax) %>%
  group_by(item_id, state_id) %>% 
  summarise(salesmax = max(salesmax)) %>% 
  NROW()

# filter low-sale items, Only 1556/9147 hierarchies are left.
dtarget <- dtmeta %>% select(item_id, state_id, salesmax) %>%
  group_by(item_id, state_id) %>%
  summarise(salesmax = max(salesmax), .groups = "drop") %>%
  filter(salesmax <= 5) %>%
  arrange(salesmax) %>%
  select(item_id, state_id)


dtarget <- dt %>% right_join(dtarget, by = c("item_id", "state_id")) %>%
  select(item_id, state_id, store_id, starts_with("d_"))


dtarget <- dtarget %>% tidyr::nest(hierarchy = -c("state_id", "item_id")) %>%
  mutate_at("hierarchy", purrr::map, function(x){
    names <- x$store_id
    series <- t(unname(as.matrix(x[,2:NCOL(x)])))
    list(names=names, series=series)
  })

dtarget <- dtarget %>% rowwise() %>%
  mutate(salesmax = max(hierarchy$series))

saveRDS(dtarget, "data.rds")
