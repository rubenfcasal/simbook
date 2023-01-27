# -------------------------
#' Inverse transform method for discrete distribution generation
#'
#' Uses the inverse transform method with sequential search for simulating a
#' discrete/categorical probability distribution (that takes on only a finite
#' number of values) from its probability mass function.
#' @param x numeric vector giving the possible values of the discrete random variable.
#' @param prob numeric vector giving the probabilities corresponding to `x`.
#' @param n number of observations to generate.
#' @param as.factor logical; if `TRUE`, the returned vector is encoded as a factor
#' with levels `x`.
#' @returns
#' Returns a numeric vector, or a factor if `as.factor = TRUE`, with the random
#' deviates and an attribute `ncomp` with the required number of comparisons in
#' the sequential search.
#' @seealso [sample()], [rpmf.table()], [rpmf.alias()]
#' @examples
#' set.seed(1)
#' # Simulation of a binomial distribution
#' n <- 10
#' p <- 0.5
#' nsim <- 10^5
#' x <- 0:n
#' pmf <- dbinom(x, n, p)
#' rx <- rpmf(x, pmf, nsim)
#' # Relative frequency plot
#' plot(table(rx)/nsim, ylab = "Relative frequency", xlab = "Value")
#' abline(v = mean(rx))
#' # Theoretical values
#' points(x, pmf, pch = 4, col = "blue")
#' abline(v = p*n, lty = 2, col = "blue")
#' # Number of comparisons
#' ncomp <- attr(rx, "ncomp")
#' ncomp/nsim # mean number
#' sum((1:length(x))*pmf) # theoretical expected number
#' @export
rpmf <- function(x, prob = 1/length(x), n = 1000, as.factor = FALSE) {
  # Numero de comparaciones
  ncomp <- 0
  # Inicializar FD
  Fx <- cumsum(prob)
  # Simular
  X <- numeric(n)
  U <- runif(n)
  for(j in 1:n) {
    i <- 1
    while (Fx[i] < U[j]) i <- i + 1
    X[j] <- x[i]
    ncomp <- ncomp + i
  }
  if(as.factor) X <- factor(X, levels = x)
  attr(X, "ncomp") <- ncomp
  return(X)
}


# -------------------------
#' Guide table method for discrete distribution generation
#'
#' Uses the guide table aided inversion method (also know as indexed search method)
#' for simulating a discrete/categorical probability distribution (that takes on
#' only a finite number of values) from its probability mass function.
#' @inheritParams rpmf
#' @param m size of the table with starting points for the search.
#' @returns
#' Returns a numeric vector, or a factor if `as.factor = TRUE`, with the random
#' deviates and an attribute `ncomp` with the required number of comparisons in
#' the sequential search.
#' @seealso [sample()], [rpmf()], [rpmf.alias()]
#' @examples
#' set.seed(1)
#' # Simulation of a binomial distribution
#' n <- 10
#' p <- 0.5
#' nsim <- 10^5
#' x <- 0:n
#' pmf <- dbinom(x, n, p)
#' rx <- rpmf.table(x, pmf, n-1, nsim)
#' # Relative frequency plot
#' plot(table(rx)/nsim, ylab = "Relative frequency", xlab = "Value")
#' abline(v = mean(rx))
#' # Theoretical values
#' points(x, pmf, pch = 4, col = "blue")
#' abline(v = p*n, lty = 2, col = "blue")
#' # Number of comparisons
#' ncomp <- attr(rx, "ncomp")
#' ncomp/nsim # mean number
#' sum((1:length(x))*pmf) # theoretical expected number with sequential search
#' @export
rpmf.table <- function(x, prob = 1/length(x), m, n = 1000, as.factor = FALSE) {
  # Inicializar tabla y FD
  Fx <- cumsum(prob)
  g <- rep(1,m)
  i <- 1
  for(j in 2:m) {
    while (Fx[i] < (j-1)/m) i <- i + 1
    g[j] <- i
  }
  ncomp <- i - 1
  # Generar valores
  X <- numeric(n)
  U <- runif(n)
  for(j in 1:n) {
    i <- i0 <- g[floor(U[j] * m) + 1]
    while (Fx[i] < U[j]) i <- i + 1
    ncomp <- ncomp + i - i0
    X[j] <- x[i]
  }
  if(as.factor) X <- factor(X, levels = x)
  attr(X, "ncomp") <- ncomp
  return(X)
}


# -------------------------
#' Alias method for discrete distribution generation
#'
#' Uses alias method (with the Robin Hood setup algorithm)
#' for simulating a discrete/categorical probability distribution (that takes on
#' only a finite number of values) from its probability mass function.
#' @inheritParams rpmf
#' @returns
#' Returns a numeric vector, or a factor if `as.factor = TRUE`, with the random
#' deviates.
#' @seealso [sample()], [rpmf()], [rpmf.table()]
#' @examples
#' set.seed(1)
#' # Simulation of a binomial distribution
#' n <- 10
#' p <- 0.5
#' nsim <- 10^5
#' x <- 0:n
#' pmf <- dbinom(x, n, p)
#' rx <- rpmf.alias(x, pmf, nsim)
#' # Relative frequency plot
#' plot(table(rx)/nsim, ylab = "Relative frequency", xlab = "Value")
#' abline(v = mean(rx))
#' # Theoretical values
#' points(x, pmf, pch = 4, col = "blue")
#' abline(v = p*n, lty = 2, col = "blue")
#' @export
rpmf.alias <- function(x, prob = 1/length(x), n = 1000, as.factor = FALSE) {
  # Inicializar tablas
  a <- numeric(length(x))
  q <- prob*length(x)
  low <- q < 1
  high <- which(!low)
  low <- which(low)
  while (length(high) && length(low)) {
    l <- low[1]
    h <- high[1]
    a[l] <- h
    q[h] <- q[h] - (1 - q[l])
    if (q[h] < 1) {
      high <- high[-1]
      low[1] <- h
    } else low <- low[-1]
  } # while
  # Generar valores
  V <- runif(n)
  i <- floor(runif(n)*length(x)) + 1
  X <- x[ ifelse( V < q[i], i, a[i]) ]
  if(as.factor) X <- factor(X, levels = x)
  return(X)
}
