#' Plot a cartogram object
#'
#' @param x a cartogram object
#'
#' @param   \\dots arguments passed to or from other methods.
#' @return No return value, called for side effects
#'
#' @export
#' @examples
#' \donttest{
#'   data(usa)
#'   carto <- cartogramR(usa, "electors64")
#'   plot(carto)
#' }
#'
#' @import sf
#' @md
plot.cartogramR <- function(x, ...) {
 if  (!inherits(x, "cartogramR")) stop(paste(deparse(substitute(x)), "must be a cartogramR object"))
 x <- x$cartogram
 plot(x)
 return(invisible(NULL))
}
