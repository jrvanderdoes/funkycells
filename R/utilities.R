#' Get Mode of vector
#'
#' @param x Vector of some items
#'
#' @return The mode(s) of the vector
#' @export
#'
#' @examples
#' .getMode(c(1,2,3,1))
#' .getMode(c('A','B','C','A','B'))
#' .getMode(c(1,2,3,'A','A'))
#' .getMode(c(1,2,3,'A','A',2,2,2))
.getMode <- function(x) {
  a <- table(x)
  mode <- names(a)[a == max(a)]
  if(is.numeric(x))
    mode <- as.numeric(mode)
  mode
  #uniqx <- unique(na.omit(x))
  #uniqx[which.max(tabulate(match(x, uniqx)))]
}
