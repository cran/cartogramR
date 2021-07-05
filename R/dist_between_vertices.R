#' Analyse some of the grid options
#'
#' @param data a sf object to be used in cartogram.
#'
#' @return a `dbv.cartogramR` object which is a data-table 
#'         which contains distance between vertices (`dbv`) and polygons names
#'         (`L1`, `L2`, `L3`) inherited from [sf::st_coordinates]
#' @examples
#'   data(usa)
#'   dbv <- dist_between_vertices(data=usa)
#'   summary(dbv)
#' @export
dist_between_vertices <- function(data) {
    L1 <- L2 <- L3 <- X <- Y <- NULL # due to NSE notes in R CMD check
    y_geom <- sf::st_geometry(data)
    coord <- data.table::data.table(data.frame(sf::st_coordinates(y_geom)))
    dista2  <- function(x,y) return(sqrt(diff(x)^2+diff(y)^2))
    dbv <- coord[, list(dbv=dista2(X,Y)), by=list(L3,L2,L1)]
    class(dbv) <- c("dbv.cartogramR","data.table","data.frame")
    return(dbv)
}
