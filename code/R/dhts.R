#' discrete hierarchical time series constructor.
#' 
#' Create a 'dhts' object 
#' @param domain the domain of bottom-level series, 2 * m matrix.
#' Each column corresponds to a series, and the first row corresponds to
#' the minimums, and the second row corresponds the maximums.
#' @param s_mat the summing matrix
#' @param bts bottom-level time series, T * m
#' @return a dhts object
#' @examples 
#' # create a 3-nodes discrete hts.
#' bts <- matrix(sample(0:1, 100, replace=TRUE), ncol=2)
#' s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
#' domain <- matrix(c(1, 0, 1, 0), 2)
#' dhts(bts, s_mat, domain)
#' @export
dhts <- function(bts, s_mat, domain){
  stopifnot(dim(bts)[2] == dim(s_mat)[2],
            dim(domain)[2] == dim(s_mat)[2],
            dim(domain)[1] == 2)
  val <- list()
  val$bts <- bts
  val$domain <- domain
  val$s_mat <- s_mat
  class(val) <- "dhts"
  val
}

#' aggregate discrete hierarchical time series.
agg.dhts <- function(dhts){
  t(dhts$s_mat %*% t(dhts$bts))
}



