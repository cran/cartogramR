#' Coerce a cartogramR to a sfc object
#'
#' Coerce a cartogramR to a  sfc object extracting the component cartogram of
#'  the cartogramR object
#'
#' @param x a cartogramR object
#' @param   \\dots arguments passed to or from other methods.
#'
#' @return a sfc object
#'
#' @rdname as.sfc.cartogramR
#' @export
#'
#' @md
as.sfc.cartogramR <- function(x, ...) {
    if  (!inherits(x, "cartogramR")) stop(paste(deparse(substitute(x)), "must be a cartogramR object"))
    return(x$cartogram)
}
