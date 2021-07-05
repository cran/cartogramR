## plot

#' Plot a gridanalysis.cartogram object
#'
#' @param x a gridanalysis.cartogram object
#'     
#' @param nthsmallest plot only  the nthsmallest values among all polygons
#' @param redrawxaxis if `TRUE` redraw ticks and labels of x axe at grid size on
#'         log scale
#' @param type character string (length 1 vector) or vector of 1-character
#'           strings indicating the type of plot for each polygons,
#'           see [graphics::matplot] for all possible `type`s.
#' @param xlab titles for x axe, as in [graphics::matplot].
#' @param ylab titles for y axe, as in [graphics::matplot].
#' @param ylim ranges of y axe, as in [graphics::matplot].
#' @param \\dots arguments passed to or from other methods.
#' @return No return value, called for side effects
#'     
#' @export
#' @examples
#' \donttest{
#'   data(usa)
#'   precarto <- precartogramR(usa, method="gsm", pf=1.2, verbose=TRUE)
#'   plot(precarto)
#' }
#'
#' @md
plot.gridanalysis.cartogramR <- function(x, nthsmallest=5, redrawxaxis= TRUE,
                                         type="b", xlab= NULL, ylab= NULL,
                                         ylim=c(0, 20), ...) {
 if  (!inherits(x, "gridanalysis.cartogramR"))
     stop(paste(deparse(substitute(x)), "must be a gridanalysis.cartogramR object"))
 if  (!is.numeric(nthsmallest) || !is.vector(nthsmallest) || (length(nthsmallest)>1))
     stop("nthsmallest must be a numeric vector of length one")
 LL <- as.numeric(substr(colnames(x), 2, nchar(colnames(x))))
 sel <- order(x[, ncol(x)])[1:nthsmallest]
   if (is.null(xlab)) 
        xlab <- "Grid size"
    if (is.null(ylab)) 
        ylab <- "Min. of counts by polygon"
 if (redrawxaxis) {
     graphics::matplot(LL, t(x[sel, ,drop=FALSE]), type=type, xlab = xlab, 
                       ylab = ylab, ylim=ylim, log="x", xaxt="n", ...)
     graphics::axis(1, at= LL, labels= as.character(LL)) 
    
 } else{
     graphics::matplot(LL, t(x[sel, ,drop=FALSE]), xlab = NULL, 
                       ylab = NULL, xlatype=type, ylim=ylim, ...)
     }
 return(invisible(NULL))
}
