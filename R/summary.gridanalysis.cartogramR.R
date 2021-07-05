## summary.gridanalysis.cartogramR

#' Summary of a gridanalysis.cartogram object  
#'
#' @param object a gridanalysis.cartogramR object
#' @param \\dots arguments passed to or from other methods.
#' @return A vector which indicate the grid size necessary to have more than `steps` 
#'        grid points in each polygon
#'     
#' @rdname summary.gridanalysis.cartogramR
#' @examples
#'   data(usa)
#'   ga <- grid_analysis(data=usa, gridpower2=4:9)
#'   summary(ga)
#'
#' @export
#' @md
summary.gridanalysis.cartogramR <- function(object, ...) {
    if  (!inherits(object, "gridanalysis.cartogramR"))
        stop(paste(deparse(substitute(object)), "must be a gridanalysis.cartogramR object"))
    apply(object, 2, function(y) summary(y, ...))
}
