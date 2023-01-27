#' @name simres-internals
#' @title simres internal and secondary functions
#' @description Listed below are supporting functions for the major methods in simres.
# @param
# @export
#' @keywords internal
NULL


#' @description `eval.curve` evaluates the curve (or expression) `expr` in the environment
#' specified by `enclos`.
#' @param expr name of the function (or call or expression written as a function)
#' of `x` which will evaluate to an object of the same length as `x`.
#' See for instance [`curve`].
#' @param envir the environment in which `expr` is to be evaluated.
#' @param enclos an *enclosure* where looks for objects not found in `envir`.
#' @rdname simres-internals
#' @keywords internal
#' @export
eval.curve <- function(expr, x, enclos = parent.frame()){
  y <- eval(expr, envir = list(x = x), enclos = enclos)
  if (length(y) != length(x))
    stop("'expr' did not evaluate to an object of length 'n'")
  return(y)
}
