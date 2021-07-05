#' Coerce an object to a sfc object
#'
#'
#' @param x  object to be coerced
#' @param   \\dots arguments passed to or from other methods.
#'
#' @return a sfc object
#'
#' @rdname as.sfc
#' @export
#'
#' @md
as.sfc <- function(x, ...) {
  UseMethod("as.sfc")
}
