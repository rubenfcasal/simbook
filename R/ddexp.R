# ------------------------
#' The double-exponential distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the double-exponential distribution with rate `rate`.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be
#' the number required.
#' @param rate vector of rates.
#' @details
#' If `rate` is not specified, it assumes the default value of 1.
#'
#' The double-exponential distribution with rate \eqn{\lambda} has density
#' \deqn{f(x)  =\frac{\lambda}{2}e^{-\lambda | x |}}
#' for \eqn{x \in \mathbb{R}}.
#' The cummulative distribution is:
#' \deqn{F(x) = \int_{-\infty}^{x}f(t) dt
#' =\left\{
#' \begin{array}{ll}
#' \frac{1}{2}e^{\lambda x} & \text{if } x < 0\\
#' 1-\frac{1}{2}e^{-\lambda x} & \text{if } x \geq 0
#' \end{array}
#' \ \right.}
#' The inverse cumulative distribution function is given by
#' \deqn{F^{-1}(p) = - \operatorname{sign}(p-0.5)\frac{\ln(1 - 2 | p - 0.5 | )}
#' {\lambda}.}
#'
#' @returns
#' `ddexp()` gives the density, `pdexp()` gives the distribution function,
#' `qdexp()` gives the quantile function, and `rdexp()` generates random deviates.
#' @seealso [dexp()]
#' @examples
#' set.seed(54321)
#' rate <- 2
#' rx <- rdexp(10^3, rate)
#' hist(rx, breaks = "FD", freq = FALSE)
#' lines(density(rx), lty = 2, col = "blue")
#' curve(ddexp(x, rate), lwd = 2, add = TRUE)
#' p <- c(0.005, 0.025, 0.05, 0.1)
#' abline(v = qdexp(c(p, 1-p), rate = rate), lty = 3, col = "red")
#' plot(ecdf(rx))
#' curve(pdexp(x, rate), lwd = 2, add = TRUE)
#' @export
# PENDENTE:
# - log = FALSE
# - lower.tail = TRUE, log.p = FALSE
# Densidad doble exponencial
ddexp <- function(x, rate = 1){
  # rate * exp(-rate * abs(x))/2
  0.5 * stats::dexp(abs(x), rate = rate)
}

# ------------------------
#' @rdname ddexp
#' @export
# Distribucion doble exponencial
pdexp <- function(q, rate = 1){
  sgn <- 1 - 2 * (q < 0)
  0.5 * (1 + sgn * stats::pexp(abs(q), rate = rate))
}

# ------------------------
#' @rdname ddexp
#' @export
# Funcion cuantil doble exponencial
qdexp <- function(p, rate = 1){
  sgn <- 1 - 2*(p < 0.5)
  # - sgn * log(1 - 2 * abs(p - 0.5))/rate
  sgn * stats::qexp(2  *abs(p - 0.5), rate = rate)
}


# ------------------------
#' @rdname ddexp
#' @export
# Simulacion n valores de doble exponencial
rdexp <- function(n = 1000, rate = 1) {
    x <- stats::rexp(n, rate = rate)
    sgn <- 1 - 2 * (stats::runif(n) < 0.5)
    return(sgn * x)
}
