#' Analyse some of the grid options
#'
#' @param data a sf object to be used in cartogram.
#' @param gridpower2 a vector of exponent (to be raised at the power of 2) that
#'    gives the log2(size) of the grid (default to `8:11`)
#' @param pf Determines space between map and boundary (default to 1.5)
#' @param verbose a boolean object to set on verbose mode (default to `TRUE`)
#'
#' @return a `gridanalysis.cartogramR` object which is a matrix 
#' @examples
#'   data(usa)
#'   ga <- grid_analysis(data=usa, gridpower2=4:8, verbose=TRUE)
#'   summary(ga)
#' @export
grid_analysis <- function(data, gridpower2=8:11, pf=1.5, verbose=FALSE) {
    y_geom <- sf::st_geometry(data)
    if (!is.numeric(pf) || (length(pf)>1) ) stop("pf is not numeric")
    if (pf<1) stop("pf must greater or equal to 1")
    if (!is.numeric(gridpower2)) stop("gridpower2 is not numeric")
    if (!is.vector(gridpower2)) stop("gridpower2 must be a vector")
    if (min(gridpower2)<=0) stop("minimum of gridpower2 must be greater than 0")
    maxgridsize <- 2^max(gridpower2)
    if (max(gridpower2)>15) warning(paste("the maximum grid size requested is: ", maxgridsize,"x", maxgridsize, sep=""))
    bbox <- sf::st_bbox(y_geom)
    LL <- 2^gridpower2
    howmuch <- matrix(0, nrow=nrow(data), ncol=length(gridpower2))
    colnames(howmuch) <- paste("L",2^gridpower2,sep="")
    for (j in seq_along(gridpower2)) {
        if (verbose) cat(paste("starting L=2^",gridpower2[j],"\n",sep=""))
        if (verbose) cat(paste("  - making the grid ",2^gridpower2[j],"x",2^gridpower2[j],"\n",sep=""))
        grille <- .Call(carto_gridanalysis,pf , as.integer(LL[j]), bbox, as.integer(verbose))
        sf::st_crs(grille) <- sf::st_crs(y_geom)
        if (verbose) cat("  - intersection polygon/grid using sf (slow step)\n")
        pointsinpoly <- sf::st_intersects(y_geom,grille)
        howmuch[,j] <- unlist(lapply(pointsinpoly,length))
    }
    class(howmuch) <- c("gridanalysis.cartogramR","matrix","array")
    return(howmuch)
}
