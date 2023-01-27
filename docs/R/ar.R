#' Ratio-of-uniforms method for Cauchy distribution generation
#'
#' Uses the ratio-of-uniforms method for simulating a standard Cauchy distribution.
#' @param n number of observations.
#' @returns
#' Returns a numeric vector with the random deviates and an attribute `ngen`
#' with the required number of uniform generations.
#' @seealso [rcauchy()]
#' @examples
#' set.seed(1)
#' nsim <- 10^4
#' rx <- rcauchy.rou(nsim)
#' hist(rx, breaks = "FD", freq = FALSE, main = "", xlim = c(-6, 6))
#' curve(dcauchy, add = TRUE)
#' ngen <- attr(rx, "ngen")
#' {cat("Number of generations = ", ngen)
#'   cat("\nAveraged number of generations = ", ngen/nsim)
#'   cat("\nRejection rate = ", 1-nsim/ngen,"\n")}
#' @export
rcauchy.rou <- function(n) {
  # Cauchy mediante cociente de uniformes
  ngen <- n
  u <- runif(n, 0, 1)
  v <- runif(n, -1, 1)
  x <- v/u
  ind <- u^2 + v^2 > 1 # TRUE si no verifica condición
  # Volver a generar si no verifica condición
  while (le <- sum(ind)){ # mientras le = sum(ind) > 0
    ngen <- ngen + le
    u <- runif(le, 0, 1)
    v <- runif(le, -1, 1)
    x[ind] <- v/u
    ind[ind] <- u^2 + v^2 > 1 # TRUE si no verifica condición
  }
  attr(x, "ngen") <- ngen
  return(x)
}



