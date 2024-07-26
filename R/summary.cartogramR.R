#' Summary of a cartogram object
#'
#' @param object a cartogramR object
#' @param   \\dots arguments passed to or from other methods. The following
#'   arguments are available:
#'     - digits integer, used for number formatting with signif if
#'       not specified (i.e., `[missing](.)`, `[signif]()` will not be
#'       called anymore (since \\R >= 3.4.0, where the default has been
#'       changed to only round in the `print` and `format` methods).
#'     - quantile.type integer code used in `quantile(*, type=quantile.type)`.
#'     - center character string code used in [residuals.cartogramR].
#'     - initial_data; the initial sf object given as input of cartogramR.
#'       Only needed for symmetric differences residuals.
#' @return A summary.cartogramR object: a list with the following components:
#'  - qrr, the summary of absolute relative residuals
#'  - qres, the summary of absolute residuals
#'  - qsymdiff, the summary of all pairwise symmetric difference beween two
#'    scaled (multi)polygons representative of two regions. These residuals
#'    are calculated only if `initial_data` argument is provided.
#'
#' @rdname summary.cartogramR
#' @export
#' @examples
#' \donttest{
#'   data(usa)
#'   carto <- cartogramR(usa, "electors64")
#'   summary(carto)
#' }
#'
#' @md
summary.cartogramR <- function(object, ...) {
  if  (!inherits(object, "cartogramR")) stop(paste(deparse(substitute(object)), "must be a cartogramR object"))
  default_options <- list(digits=NULL, quantile.type = 7, center = "best")
  argsup <- eval(substitute(alist(...)))
  req_options <- names(argsup)
  default_options[req_options] <- argsup
  rr <- residuals.cartogramR(object, type="relative error")
  qrr <- stats::quantile(abs(rr), names = FALSE, type = default_options$quantile.type)
  qrr <- c(qrr[1L:3L], mean(abs(rr)), qrr[4L:5L])
  if(!is.null(default_options$digits)) qrr <- signif(qrr, default_options$digits)
  names(qrr) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  res <- residuals.cartogramR(object, type="absolute error")
  qres <- stats::quantile(abs(res), names = FALSE, type = default_options$quantile.type)
  qres <- c(qres[1L:3L], mean(abs(res)), qres[4L:5L])
  if(!is.null(default_options$digits))  qres <- signif(qres, default_options$digits)
  names(qres) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  initial_data <- eval(parse(text=attr(object,"initial_data_name")),
                             envir = parent.frame())
  if (!inherits(initial_data, "sf")) stop(parse("Initial data object",
                                                attr(object,"initial_data_name"),
                                          "not found in parent frame"))
  if (default_options$center!="best") {
    delta <- residuals.cartogramR(object, type="symmetric difference", center=default_options$center)
  } else {
    if (attr(object, "method")%in%c("gsm", "gn"))    {
      Delta <- matrix(0,nrow=length(res),3)
      Delta[,3] <- residuals.cartogramR(object, type="symmetric difference", center="deformed_center")
    } else {
      Delta <- matrix(0,nrow=length(res),2)
    }
    Delta[,1] <- residuals.cartogramR(object, type="symmetric difference", center="point_on_surface")
    Delta[,2] <- residuals.cartogramR(object, type="symmetric difference", center="centroid")
    delta <- apply(Delta,1,min)
  }
  qdelta <- stats::quantile(delta, names = FALSE, type = default_options$quantile.type)
  qdelta <- c(qdelta[1L:3L], mean(delta), qdelta[4L:5L])
  if(!is.null(default_options$digits))  qdelta <- signif(qdelta, default_options$digits)
  names(qdelta) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  ans <- list(qrr=qrr, qres=qres, qsymdiff=qdelta)
  class(ans) <- "summary.cartogramR"
  return(ans)
}
