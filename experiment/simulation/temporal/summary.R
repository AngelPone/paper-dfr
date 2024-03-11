library(DiscreteRecon)
recf <- readRDS('experiment/simulation/results/temporal_recf.rds')

s_mat <- rbind(rep(1, 7), diag(7))
domain <- rbind(rep(0, 7), rep(1, 7))

lst_bs <- list()
hierarchy <- dhier(s_mat, domain)


ns <- c("base", "bu", "td", "sdfr","emp")

### Calculate brier score ####
for (method in c("emp")){
  res <- NULL
  for (i in 1:1000) {
    x <- recf[[i]]
    fcasts <- x$fcasts[[method]]
    y <- x$y
    if (method == "emp") {
      fcasts <- matrix(fcasts, byrow = TRUE, ncol = 128, nrow = dim(y)[1])
    }
    bs <- brier_score(fcasts, y, hierarchy)
    bs <- c(sum(bs$hierarchy), bs$series)
    res <- rbind(res, bs)
  }
  lst_bs[[method]] <- res
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

write.csv(format(tb, digits=2, nsmall=2), "manuscript/figures/simulation-temporal_bs.csv")


