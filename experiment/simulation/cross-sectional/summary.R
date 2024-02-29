res <- readRDS("experiment/simulation/results/cs-res.rds")

# filter out errors
res <- Filter(function(x){!is.null(x$metric)}, res)

# summarise
accs_list <- list()
for (m in c("base", "bu", "td", "dfr", "emp")){
  accs_list[[m]] <- sapply(res, function(x){
    x$metric[[m]]
  })
}
accs_sum <- data.frame(lapply(accs_list, rowMeans))
row.names(accs_sum) <- c("y_3", "y_1", "y_2", "Y")

write.csv(format(accs_sum*100, digits=2, nsmall=2), "manuscript/figures/simulation-cs-bs.csv")


library(tsutils)

plot_data <- list()
for (i in 1:4){
  plot_data[[i]] <- sapply(accs_list, function(x){
    x[i,]
  })
}

# prevent inaccuracy caused by floating error, the accuracy should be same
plot_data[[1]][,"td"] <- plot_data[[1]][,"base"]
plot_data[[2]][,"bu"] <- plot_data[[2]][,"base"]
plot_data[[3]][,"bu"] <- plot_data[[3]][,"base"]

pdf(file="manuscript/figures/sim_cross_mcb.pdf", width = 12, height = 4,
    pointsize = 16)
par(mfrow=c(1,3))
nemenyi(plot_data[[1]], plottype = "vmcb", 
        labels = c("Base", "DBU", "DTD", "DFR", "Empirical"), 
        main = "MCB Test for total series")
nemenyi(rbind(plot_data[[2]], plot_data[[3]]), plottype = "vmcb", 
        labels = c("Base", "DBU", "DTD", "DFR", "Empirical"), 
        main = "MCB Test for bottom series")
nemenyi(plot_data[[4]], plottype = "vmcb", 
        labels = c("Base", "DBU", "DTD", "DFR", "Empirical"), 
        main = "MCB Test for hierarchy")
dev.off()

