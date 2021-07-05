#' Coerce an object to a sf object
#'
#'
#' @param x  object to be coerced
#' @param   \\dots arguments passed to or from other methods.
#'
#' @return an sf object
#'
#' @rdname as.sf
#' @export
#'
#' @md
as.sf <- function(x, ...) {
  UseMethod("as.sf")
}
