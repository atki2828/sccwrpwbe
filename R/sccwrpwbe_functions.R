#' @title Minmax Scale
#'
#' @description Scales values of numeric vector to be between 0 and 1
#'
#' @param \code{x}     numeric vector
#'
#' @return \code{x} values from input vector scaled to be between 0 and 1
#' @export
#'
#' @examples
#' x <- runif(10, min =0 ,max=100)
#' minmax(x)
#'
#'
minmax <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}






