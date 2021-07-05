## summary
##     par L3 on fait une analyse de la

#' Summary of a dbv.cartogram object  
#'
#' @param object a dbv.cartogramR object
#' @param \\dots arguments passed to or from other methods.
#' @return a data-table which contains by region (`L3`)
#'    - the sample quantiles corresponding to the probability
#'       0.8, 0.85, ...,1
#'    - the total number of vertices divided by the 
#'       perimeter of the region (the sum of all polygons
#'       perimeter of the region, `NbyPerim`)
#' @rdname summary.dbv.cartogramR
#' @examples
#'   data(usa)
#'   dbv <- dist_between_vertices(data=usa)
#'   summary(dbv)
#'
#' @importFrom stats quantile
#' @export
#' @md
summary.dbv.cartogramR <- function(object, ...) {
    L3 <- dbv <- NULL # due to NSE notes in R CMD check
    if  (!inherits(object, "dbv.cartogramR"))
        stop(paste(deparse(substitute(object)), "must be a dbv.cartogramR object"))
    object2 <- copy(object)
    return(setorder(object2, -dbv)[1:10,])
}
