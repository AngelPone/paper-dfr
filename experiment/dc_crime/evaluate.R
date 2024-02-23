library(DiscreteRecon)

res <- readRDS("experiment/dc_crime/data/ct_results.rds")
data <- readRDS("experiment/dc_crime/data/basef.rds")

SMAT <- rbind(rep(1, 4), diag(4))

bs <- lapply(seq_along(res), function(idx) {
  input <- data[[idx]]
  fcasts <- res[[idx]]
  if (max(input$y) >= 2) {
    DOMAIN <- rbind(0, rep(2, 4))
    max_domain <- 2
  } else {
    DOMAIN <- rbind(0, rep(1, 4))
    max_domain <- 1
  }
  T_ <- NROW(input$y)
  ht <- dhier(SMAT, DOMAIN)

  basef <- lapply(1:5, function(x) {
    max_d <- ifelse(x == 1, max_domain * 4, max_domain)
    bf <- lapply(input$fcasts[, x], function(lambda) {
      if (lambda < 0) {
        lambda <- 0.00001
      }
      c(dpois(0:(max_d - 1), lambda), 1 - ppois(max_d - 1, lambda))
    })
    do.call(rbind, bf)
  })

  test_basef <- lapply(basef, function(x) {
    output <- x[(T_ - 15):T_, ]
    colnames(output) <- 0:(NCOL(output) - 1)
    output
  })

  fcasts$base <- test_basef
  fcasts$emp <- matrix(fcasts$emp, nrow = 16, ncol = NCOL(fcasts$emp), byrow=TRUE)
  input$y[input$y > max_domain] <- max_domain
  test_y <- input$y[(T_ - 15):T_, ]
  lapply(fcasts, function(m) {
    brier_score(m, test_y, ht)
  })
})




get_fcasts <- function(getf) {
  methods <- c("dfr", "base", "bu", "td", "emp")
  output <- lapply(methods, function(m){
    do.call(c, lapply(bs, function(b) {
      getf(b[[m]])
    }))
  })
  output <- do.call(cbind, output)
  colnames(output) <- c("DFR", "Base", "DBU", "DTD", "Empirical")
  output
}


b_res <- get_fcasts(function(x){ x$series[2:5]} )

library(tsutils)
b_res[, "DBU"] <- b_res[, "Base"]

pdf("manuscript/figures/dc_crime_mcb.pdf",width = 15, height = 5, pointsize = 16)
par(mfrow = c(1, 3))
nemenyi(b_res,
  plottype = "vmcb",
  main = sprintf("MCB Test on %s of %s", "Brier Score", "bottom-level series")
)

t_res <- get_fcasts(function(x){ x$series[1]} )
t_res[, "DTD"] <- t_res[, "Base"]
nemenyi(t_res,
  plottype = "vmcb",
  main = sprintf("MCB Test on %s of %s", "Brier Score", "total series")
)
colMeans(t_res)

hierarchy_res <- get_fcasts(function(x){sum(x$hierarchy)})

nemenyi(hierarchy_res,
  plottype = "vmcb",
  main = sprintf("MCB Test on %s of %s", "Brier Score", "hierarchy")
)
colMeans(hierarchy_res)

dev.off()


library(dplyr)

get_summarised_bs <- function(f) {
  lapply(list(t_res, b_res, hierarchy_res), function(res){
    apply(res, 2, f) * 100
  }) %>%
    do.call(rbind, .) %>%
    round(digits = 2) %>%
    `row.names<-`(c("Total", "Bottom", "Hierarchy"))
}

get_summarised_bs(mean) %>%
  write.csv("manuscript/figures/dc_crime_mean.csv")

