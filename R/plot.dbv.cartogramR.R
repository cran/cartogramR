#' Plot a dbv.cartogram object
#'
#' @param x a dbv.cartogram object
#' @param which if a subset of the plots is required, specify a subset of the
#'          numbers `1:2`
#' @param ask logical; if `TRUE`, the user is asked before each plot, see
#'          [par]`(ask=.)`
#' @param key logical; if `TRUE`, a legend is drawn
#' @param last draw the density of distance between vertices for the `last`
#'    coordinates
#' @param probminx the sample quantiles (of distance between vertices)
#'      corresponding to the probability is used as a minimum
#'      of x-axis for the density plot (used only if `last` is `NULL`)
#' @param \\dots arguments passed to or from other methods.
#' @return No return value, called for side effects
#' @details The first plot is the density of distance between
#'     consecutive vertice by region. Only the upper quantiles
#'     are shown. The second plot is a barplot by region of the
#'     number of vertice divided by the perimeter of the region
#'
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics barplot lines legend
#' @importFrom stats density
#' @export
#' @examples
#'   data(usa)
#'   precarto <- precartogramR(usa, method="dcn")
#'   plot(precarto)
#'
#' @md
plot.dbv.cartogramR <- function(x, which=1:2, ask=TRUE, key=TRUE, last=10, probminx=0.9, ...) {
    L3 <- dbv <- NULL # due to NSE notes in R CMD check
    if  (!inherits(x, "dbv.cartogramR"))
        stop(paste(deparse(substitute(x)), "must be a dbv.cartogramR object"))
    x2 <- data.table::copy(x)
    show <- rep(FALSE, 2)
    show[which] <- TRUE
    if (sum(show)>1) {
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
    }
    if (show[1L]) {
        maxi <- xx <- yy <- GRPE <- NULL # due to NSE notes in R CMD check
        dd <- function(x) {
            prov <- density(x)
            return(list(xx=prov$x,yy=prov$y))
        }
        ddd <- x2[, dd(dbv), by=L3]
        if (is.null(last)) {
            if ((probminx>=1)|(probminx<0)) stop("probminx should be in [0, 1)")
            minx <- x2[, quantile(dbv, probs=probminx, na.rm=TRUE)]
        } else {
            data.table::setorder(x2, -dbv)
            minx <- x2[last, dbv]
        }
        coords <- ddd[xx>=minx, ]
        coords[, ":="("GRPE"=.GRP), by=list(L3)]
        maxx <- coords[, max(xx, na.rm=TRUE)]
        maxy <- coords[, max(yy, na.rm=TRUE)]
        plot(c(minx, maxx), c(0, maxy), type="n", xlab="distance between consecutive vertices", ylab="density", ...)
        coords[, lines(xx, yy, col=GRPE[1], ...), by=list(L3)]
        if (key) {
            L3sel <- coords[, list(maxi=unique(L3))][,maxi]
            colL3sel <- coords[, list(maxi=unique(GRPE))][,maxi]
            ncols <- max(ceiling(length(L3sel)/20),2)
            legend("topright", legend=L3sel, lty = 1,
                   col = colL3sel, bty = "n", ncol = ncols)
        }
    }
    if (show[2L]) {
        eff <- NULL # due to NSE notes in R CMD check
        barplot(x2[, list(eff=.N/sum(dbv)), by=list(L3)][, eff], ...)
    }
}
