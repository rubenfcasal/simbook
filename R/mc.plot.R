#' Exploratory plots of a random sequence
#'
#' `mc.plot()` draws the approximation of the distribution, the
#' convergence plot (by calling `conv.plot()`), a normal QQ plot
#' and a sequential plot.
#' @param x the simulated values.
#' @param level the confidence level required.
#' @param true.value the theoretical value.
#' @param main an overall title for the plot.
#' @param omd a vector of the form `c(x1, x2, y1, y2)` giving the region inside
#' outer margins in normalized device coordinates (i.e. fractions of the device
#' region).
#' @param ... further arguments passed to other functions (e.g. to `conv.plot()`).
#' @return Return, invisibly in the case of plot functions, the approximation by
#' simulation (the arithmetic mean) and the corresponding error range (half
#' width of the confidence interval).
#' @examples
#' set.seed(1)
#' teor <- 0
#' res <- mc.plot(rnorm(1000, mean = teor), true.value = teor)
#' res
#' @export mc.plot
mc.plot <- function(x, level = 0.95, true.value = NULL, main,
                    omd = c(0.05, 0.95, 0.01, 0.95), ...){
  old.par <- par(mfrow = c(2, 2), omd = omd)
  on.exit(par(old.par))
  # Aproximacion distribución
  hist(x, breaks = "FD", freq = FALSE,
     main = "", xlab = "Simulated values")
  lines(density(x))
  abline(v = mean(x))
  abline(v = quantile(x, probs = (1 + level * c(-1, 1))/2), lty = 3)
  if (!is.null(true.value))
    abline(v = true.value, lty = 2, lwd = 2, col = "red")
  # Grafico de converxencia
  result <- conv.plot(x, level = level, ...)
  if (!is.null(true.value))
    abline(h = true.value, lty = 2, lwd = 2, col = "red")
  # Grafico QQ
  qqnorm(x, main = "", xlab = "Quantiles of Standard Normal")
  qqline(x)
  # Grafico secuencial
  plot(x, type = 'l', ylab = "Simulated values")
  par(old.par)
  # Titulo
  if (missing(main))
    main <- paste("Exploratory plots of", deparse(substitute(x)))
  title(main = main)
  return(invisible(result))
}


#' @rdname mc.plot
#' @description `conv.plot()` draws a convergence plot.
#' @param lty a vector of line types (of the form `c(conv, value, error)`).
#' @param lwd a vector of line widths (of the form `c(conv, value, error)`).
#' @param ylim the y limits of the plot.
#' @param xlab,ylab the axis titles.
#' @examples
#' set.seed(1)
#' p <- 0.4
#' res <- conv.plot(rbinom(1000, size = 1, prob = p))
#' abline(h = p, lty = 2, col = "red") # Theoretical value
#' res
#' @export conv.plot
conv.plot <- function(x, level = 0.95, lty = c(conv = 1, value = 1, error = 3),
      lwd = c(conv = 2, value = 1, error = 1), ylim = NULL,
      xlab = "Number of generations", ylab = "Mean and error range", ...) {
# PENDENTE: n.burn = 0
  nsim <- length(x)
  n <- 1:nsim
  est <- cumsum(x)/n
  # (cumsum(x^2) - n*est^2)/(n-1) # Varianzas muestrales
  valor <- est[nsim]
  esterr <- sqrt((cumsum(x^2)/n - est^2)/(n-1)) # Errores estándar
  alpha <- (1 - level)/2
  q <- qnorm(1 - alpha)
  max.error <- q * esterr[nsim]
  if (is.null(ylim)) ylim <- valor + max.error * c(-5, 5)
  plot(est, type = "l", lty = lty[1], lwd = lwd[1], xlab = "Number of generations",
       ylab = "Mean and error range", ylim = ylim, ...)
  abline(h = valor, lty = lty[2], lwd = lwd[2])
  lines(est + q * esterr, lty = lty[3], lwd = lwd[3])
  lines(est - q * esterr, lty = lty[3], lwd = lwd[3])
  return(invisible(list(approx = valor, max.error = max.error)))
}


#' @rdname mc.plot
#' @description `mc.integral()` integrates an one-dimensional function
#' over a bounded interval using classic Monte-Carlo integration and draws
#' the corresponding convergence plot.
#' @param fun an one-dimensional function to be integrated on \[a, b\]
#' @param a,b the limits of integration (must be finite).
#' @param n number of uniform generations.
#' @param plot logical; if `TRUE` a convergence plot is draw.
#' @examples
#' fun <- function(x) ifelse((x > 0) & (x < 1), 4 * x^4, 0)
#' curve(fun, 0, 1)
#' abline(h = 0, lty = 2)
#' abline(v = c(0, 1), lty = 2)
#' set.seed(1)
#' mc.integral(fun, 0, 1, 1000)
#' abline(h = 4/5, lty = 2, col = "red") # Theoretical value
#' set.seed(1)
#' mc.integral(fun, 0, 1, 5000, plot = FALSE)
#' @export mc.integral
mc.integral <- function(fun, a, b, n, level = 0.95, plot = TRUE, ...) {
  fx <- sapply(runif(n, a, b), fun) * (b - a)
  result <- if (plot) conv.plot(fx, level = level, ...) else {
    q <- qnorm((1 + level)/2)
    list(approx = mean(fx), max.error = q * sd(fx)/sqrt(n))
  }
  return(result)
}

