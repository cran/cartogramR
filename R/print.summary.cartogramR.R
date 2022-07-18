#' Print a summary of a cartogram object
#'
#' @param x a summary.cartogramR object
#' @param   \\dots arguments passed to or from other methods. The following
#'   argument is available at this level :  digits, the number of significant
#'   digits to use when printing.
#' @return x.
#'
#' @rdname print.summary.cartogramR
#' @export
#' @examples
#' \donttest{
#'   data(usa)
#'   carto <- cartogramR(usa, "electors64")
#'   summary(carto)
#' }
#'
#' @md
print.summary.cartogramR <- function(x, ...) {
    if  (!inherits(x, "summary.cartogramR")) stop(paste(deparse(substitute(x)), "must be a summary.cartogramR object"))
    default_options <- list(digits=max(3L, getOption("digits") - 3L))
    argsup <- eval(substitute(alist(...)))
    req_options <- names(argsup)
    select <- names(req_options) %in% names(default_options)
    default_options[req_options] <- argsup
    cat("* Objective:\n")
    cat("  - Relative error:\n")
    if (!all(select)) {
       ll <- c(list(x=x$qrr), default_options, req_options[!select])
       do.call(print, ll)
    } else {
       print(x$qrr, digits = default_options$digits)
    }
    cat("  - Absolute error:\n")
    if (!all(select)) {
       ll <- c(list(x=x$qres), default_options, req_options[!select])
       do.call(print, ll)
    } else {
       print(x$qres, digits = default_options$digits)
    }
    cat("\n* Distortion of regions\n")
    cat("  Symmetrized difference error\n")
    if (!all(select)) {
       ll <- c(list(x=x$qsymdiff), default_options, req_options[!select])
       do.call(print, ll)
    } else {
       print(x$qsymdiff, digits = default_options$digits)
    }
    invisible(x)
}
