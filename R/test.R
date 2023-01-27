# ------------------------------------------------------------------------------
#' Chi-squared goodness-of-fit test for continuous data
#'
#' Performs the chi-squared goodness-of-fit test for a continuous distribution
#' by grouping data into *bins* (also called *intervals*, *classes* or *cells*)
#' with equal probabilities under the null hypothesis.
#' @param x numeric vector containing the observed values.
#' @param distribution character string naming a continuous distribution,
#' such as "norm" or "unif" (the cumulative distribution `p<distribution>()` and
#' quantile `q<distribution>()` functions must exist).
#' @param nclass number of bins.
#' @param output logical; if `TRUE` an histogram is plotted and a table with
#' the results for each class is printed.
#' @param nestpar number of estimated parameters (composite null hypothesis).
#' @param ...  parameters of the distribution (specified by `distribution`).
#' @return A list with class `"htest"` containing the following components:
#' \item{statistic}{the value the chi-squared test statistic.}
#' \item{parameter}{the degrees of freedom of the approximate chi-squared
#' distribution of the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{method}{a character string indicating the type of test performed.}
#' \item{data.name}{a character string with the actual `x` argument name.}
#' \item{classes}{a character vector with the class labels.}
#' \item{observed}{the observed counts.}
#' \item{expected}{the expected counts under the null hypothesis.}
#' \item{residuals}{the Pearson residuals, `(observed - expected) / sqrt(expected)`.}
#' @examples
#' nx <- 30
#' x <- rnorm(nx)
#' chisq.cont.test(x, distribution = "norm", nestpar = 2,
#'                 mean = mean(x), sd = sqrt((nx - 1) / nx) * sd(x))
#' @seealso [`chisq.test`] [`freq.test`]
#' @export
chisq.cont.test <- function(x, distribution = "norm", nclass = floor(length(x)/5),
                            output = TRUE, nestpar = 0, ...) {
  # Función distribución
  q.distrib <- eval(parse(text = paste("q", distribution, sep = "")))
  # Puntos de corte
  q <- q.distrib((1:(nclass - 1))/nclass, ...)
  tol <- sqrt(.Machine$double.eps)
  xbreaks <- c(min(x) - tol, q, max(x) + tol)
  # Gráficos y frecuencias
  if (output) {
    xhist <- hist(x, breaks = xbreaks, freq = FALSE,
                  lty = 2, border = "grey50")
    # Función densidad
    d.distrib <- eval(parse(text = paste("d", distribution, sep = "")))
    curve(d.distrib(x, ...), add = TRUE)
  } else {
    xhist <- hist(x, breaks = xbreaks, plot = FALSE)
  }
  # Cálculo estadístico y p-valor
  O <- xhist$counts  # Equivalente a table(cut(x, xbreaks)) pero más eficiente
  E <- length(x)/nclass
  DNAME <- deparse(substitute(x))
  METHOD <- "Pearson's Chi-squared test"
  STATISTIC <- sum((O - E)^2/E)
  names(STATISTIC) <- "X-squared"
  PARAMETER <- nclass - nestpar - 1
  names(PARAMETER) <- "df"
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  # Preparar resultados
  classes <- format(xbreaks)
  classes <- paste("(", classes[-(nclass + 1)], ",", classes[-1], "]",
                   sep = "")
  RESULTS <- list(classes = classes, observed = O, expected = E,
                  residuals = (O - E)/sqrt(E))
  if (output) {
    cat("\nPearson's Chi-squared test table\n")
    print(as.data.frame(RESULTS))
  }
  if (any(E < 5))
    warning("Chi-squared approximation may be incorrect")
  structure(c(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL,
                   method = METHOD, data.name = DNAME), RESULTS), class = "htest")
}

# ------------------------------------------------------------------------------
#' Frequency test
#'
#' Performs the chi-squared goodness-of-fit test for the particular case
#' of an uniform distribution in the unit interval
#' (by grouping data into *bins*, `y = floor(nclass*x) + 1`,
#' and calling [`chisq.test`]`(y)`).
#' Used for testing random number generators.
#' @param x numeric vector containing the observed values
#' (or the pseudorandom numbers).
#' @param nclass number of bins (to partition the unit interval).
#' @return A list with class `"htest"` returned by [`chisq.test()`].
#' @seealso [`chisq.cont.test`] [`chisq.test`]
#' @examples
#' set.rng(321, "lcg", a = 5, c = 1, m = 512)    # set.seed(321)
#' nsim <- 500
#' res <- freq.test(rng(nsim))
#' res # Pearson's Chi-squared test
#' # Pearson's Chi-squared test table
#' sapply(res[c("observed", "expected", "residuals", "stdres")], as.vector)
#' # Plot observed and expected counts.
#' plot(res$observed)
#' abline(h = res$expected[1], lty = 2)
#' @export
freq.test <- function(x, nclass = floor(length(x)/10)){
  y <- floor(nclass*x) + 1
  # Test chi-cuadrado
  f <- table(factor(y, levels = seq_len(nclass)))
  result <- chisq.test(f)
  result$data.name <- deparse(substitute(x))
  result
}

# ------------------------------------------------------------------------------
#' Test repetition
#'
#' Applies a hypothesis test to simulated samples.
#' @param n sample size.
#' @param nsim number of simulations.
#' @param rand.gen optional: function to generate the samples.
#' @param test function (or function name) which performs an one-sample test.
#' @param ... arguments to be passed to other functions (for instance to `test()`)
#' or methods.
#' @return `rephtest()` returns a list with class `"rhtest"` containing the following components:
#' \item{statistics}{the values of the test statistic.}
#' \item{p.values}{the p-values for the test.}
#' with attributes:
#' \item{method}{a character string indicating the type of test performed.}
#' \item{names.stat}{a character string indicating the distribution
#' of the test statistic.}
#' \item{parameter}{the parameters of the distribution of the test statistic.}
#' @examples
#' set.rng(543210, "lcg", a = 2^16 + 3, c = 0, m = 2^31)  # set.seed(543210)
#' res <- rephtest(n = 30, test = chisq.cont.test, rand.gen = rng,
#'                  distribution = "unif", output = FALSE, nestpar = 0)
#' str(res)
#' summary(res)
#' plot(res)
#' @seealso [`chisq.cont.test`] [`freq.test`]
#' @export
rephtest <- function(n = 30, nsim = 1000, test,
                      rand.gen = runif, ...){
# n = 30; nsim = 1000; test = chisq.cont.test; rand.gen = runif
# distribution = "unif"; nclass = 100; output = FALSE; nestpar = 0
# min = 0; max = 1
  # if (!length(arguments <- list(...)))
  #   if (test == chisq.cont.test)
  #     arguments <- list(distribution = "unif", output = FALSE, nestpar = 0,
  #                       min = 0, max = 1)
  arguments <- list(...)
  # Resultados
  estadistico <- numeric(nsim)
  pvalor <- numeric(nsim)
  # Realizar contrastes
  for(isim in 1:nsim) {
    arguments$x <- rand.gen(n)    # Generar
    tmp <- do.call(test, arguments) # Realizar contrastes
    # stopifnot(inherits(tmp, "htest"))
    if(!inherits(tmp, "htest"))
      stopifnot(c("statistic", "p.value") %in% names(tmp))
    estadistico[isim] <- tmp$statistic
    pvalor[isim] <- tmp$p.value
  }
  res <- structure(list(statistics = estadistico, p.values = pvalor),
                   class = "rhtest")
  if(inherits(tmp, "htest"))
    attributes(res) <- c(attributes(res), with(tmp, list(method = method,
                                      names.stat = attr(statistic, "names"),
                                      parameter = parameter)))
  return(res)
}

# ------------------------
#' @rdname rephtest
#' @param object an object for which a summary is desired.
#' @param alpha numeric vector of probabilities (significance levels).
#' @export
summary.rhtest <- function(object, alpha = c(0.01, 0.05, 0.1, 0.25, 0.5), ...){
  # object <- res; alpha = c(0.01, 0.05, 0.01)
  res <- sapply(alpha, function(a) mean(object$p.values < a) )
  names(res) <- paste0(format(alpha*100, trim = TRUE, ...), "%")
  oldClass(res) <- c("summaryrhtest", oldClass(res))
  res
}

# ------------------------
#' @rdname simres-internals
#' @keywords internal
#' @export
print.summaryrhtest <- function(x, ...){
  cat("Proportion of rejections:\n")
  print(unclass(x), ...)
}

# ------------------------
#' @rdname rephtest
#' @param x an object with class `"rhtest"`.
#' @param y if a subset of the plots is required,
#' specify a subset of the numbers \code{1:3}.
#' @param ask	logical; if `TRUE`, the user is asked before each plot,
#' see [`graphics::par`]`(ask=.)`.
#' @export
plot.rhtest <- function(x, y = 1:3, ask = length(y) > 1 && dev.interactive(),
                        ...){
  show <- rep(FALSE, 3)
  show[y] <- TRUE
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if(show[1]) { # statistic
    hist(x$statistics, breaks = "FD", freq = FALSE,
         main = "Histogram of statistics", xlab = "Statistic")
  }
  if(show[2]) { # p.value
    hist(x$p.values, breaks = "FD", freq = FALSE,
         main = "Histogram of p-values", xlab = "p-value")
    abline(h = 1) # curve(dunif(x,0,1), add=TRUE)
  }
  if(show[3]) { # rejections
    ecd <- ecdf(x$p.values)
    curve(ecd, type = "s", lwd = 2, main = "Rejections",
      xlab = "Significance level (alpha)", ylab = "Proportion of rejections")
    abline(a = 0, b = 1, lty = 2)   # curve(punif(x, 0, 1), add = TRUE)
  }
}




