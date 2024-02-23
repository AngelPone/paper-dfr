library(DiscreteRecon)
recf <- readRDS('experiment/simulation/results/temporal_recf.rds')

s_mat <- rbind(rep(1, 7), diag(7))
domain <- rbind(rep(0, 7), rep(1, 7))

lst_bs <- list()
hierarchy <- dhier(s_mat, domain)


ns <- c("base", "bu", "td", "sdfr","emp")

### Calculate brier score ####
for (method in ns){
  lst_bs[[method]] <- lapply(iterators::iter(recf), 
  function(x) {
    fcasts <- x$fcasts[[method]]
    y <- x$y
    bs <- brier_score(fcasts, y, hierarchy)
    c(sum(bs$hierarchy), bs$series)
  })
  lst_bs[[method]] <- do.call(rbind, lst_bs[[method]])
}


### MCB Test #####
library(tsutils)

pdf(sprintf("manuscript/figures/temporal_mcb.pdf"),
    width = 12, height = 4, pointsize = 16)
par(mfrow=c(1,3), mar = c(5.1, 6, 4.1, 2.1))
dat <- sapply(ns, function(x){ lst_bs[[x]][,2]}, simplify="array")
dat[,"td"] <- dat[,"base"]
nemenyi(dat, 
        plottype="vmcb", 
        labels=c("Base", "DBU", "DTD", "SDFR", "Empirical"),
        main = sprintf("MCB Test for total series"))
dat <- sapply(ns, function(x){as.vector(lst_bs[[x]][,3:9])}, simplify = "array")
dat[,"bu"] <- dat[,"base"]
nemenyi(dat, 
        plottype="vmcb", 
        labels=c("Base", "DBU", "DTD", "SDFR", "Empirical"),
        main = sprintf("MCB Test for bottom series"))
nemenyi(sapply(ns, function(x){lst_bs[[x]][,1]}, simplify = "array"), 
        plottype="vmcb", 
        labels=c("Base", "DBU", "DTD", "SDFR", "Empirical"),
        main = sprintf("MCB Test for hierarchy"))
dev.off()


tb <- sapply(ns, function(x){
  colMeans(lst_bs[[x]]) * 100
})

rownames(tb) <- c("Y", paste0("Y", 1:8))

write.csv(format(tb, digits=2, nsmall=2), "experiment/simulation/results/temporal_bs.csv")
# sapply(ns, function(x){sapply(lst_bs[[x]], function(x){sum(x$hierarchy)})}, simplify = "array") %>%
#   mcb_plot("Brier Score", "hierarchy")
# sapply(ns, function(x){sapply(lst_bs[[x]], function(x){x$series[1]})}, simplify = "array") %>%
#   mcb_plot("Brier Score", "total")
# sapply(ns, function(x){as.vector(sapply(lst_bs[[x]], function(x){x$series[2:8]}))}, simplify = "array") %>%
#   mcb_plot("Brier Score", "bottom")


