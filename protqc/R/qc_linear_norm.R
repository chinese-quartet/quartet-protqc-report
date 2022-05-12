#' Linear normalization: 1~10
#'
#' @param x A value to be normalized
#' @param x_max Queried maximum in historical values
#' @param x_min Queried minimum in historical values
#' @export

qc_linear_norm <- function(x, x_min, x_max, decreasing = F) {
  if(length(x) == 1){
    if(x <= x_min) x <- x_min else if(x >= x_max) x <- x_max
  }
  if(decreasing == F) {
    x_norm <- sapply(x, function(a) 1 + (a - x_min)*9 / (x_max - x_min))
  }else {
    x_norm <- sapply(x, function(a) 10 - (a - x_min)*9 / (x_max - x_min))
  }
  return(x_norm <- round(x_norm, digits = 3))
}
