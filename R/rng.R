#' Pseudorandom number generation
#'
#' `.rng` is a list containing the state of the uniform pseudorandom number
#' generator (with components `seed`, `method` and `parameters`).
#' It is advisable to use `set.rng()` to set it.
#' @aliases .rng
#' @param seed initial seed.
#' @param type string specifying the generator.
#' @param ... generator parameters.
#' @returns `set.rng()` stores the state of the generator in `.rng` in the global
#' environment.
# (and returns it invisibly).
#' @export
#' @examples
#' print(set.rng())
#' rng(10)
#' .rng
set.rng <- function(seed = as.numeric(Sys.time()), type = c("lcg", "vng"), ...) {
  type <- match.arg(type)
  if (!length(parameters <- list(...)))
    if (type == "lcg") parameters <- list(a = 7^5, c = 0, m = 2^31 - 1)
    if (type == "vng") parameters <- list(k = 4)
  assign(".rng", list(seed = seed, type = type, parameters = parameters),
      envir = globalenv())
  # return(invisible(
  #   .rng <<- list(seed = seed, type = type, parameters = parameters)))
}


#' @rdname set.rng
#' @description `rng()` returns a sequence of uniform pseudorandom numbers
#' (using the generator selected with `set.rng()`).
#' @param n number of generations.
#' @returns `rng()`, `rlcg()` and `rvmg()` return a numeric vector.
#' @export
rng <- function(n) {
  if (!exists(".rng", envir = globalenv())) set.rng()
  switch(.rng$type,
     lcg = with(.rng, rlcg(n, seed, parameters$a, parameters$c, parameters$m)),
     vmg = with(.rng, rvng(n, seed, parameters$k)))
}

# Generador congruencial lineal
#
#' @rdname set.rng
#' @description `rlcg()` returns a sequence of uniform pseudorandom numbers
#' using the linear congruential generator (LCG).
# @param n number of generations.
# @param seed seed.
#' @param a multiplier.
#' @param c increment.
#' @param m modulus.
# @returns A numeric vector.
# @seealso [set.rng()], [rng()]
#' @export
rlcg <- function(n, seed = as.numeric(Sys.time()), a = 7^5, c = 0, m = 2^31 - 1) {
  u <- numeric(n)
  for(i in 1:n) {
    seed <- (a * seed + c) %% m
    u[i] <- seed/m # (seed + 1)/(m + 1)
  }
  # Almacenar semilla y parámetros
  assign(".rng", list(seed = seed, type = "lcg",
          parameters = list(a = a, c = c, m = m)), envir = globalenv())
  # .rng <<- list(seed = seed, type = "lcg", parameters = list(a = a, c = c, m = m))
  # Para continuar con semilla y parámetros:
  #   with(.rng, rlcg(n, seed, parameters$a, parameters$c, parameters$m))
  # Devolver valores
  return(u)
}

# Generador de Von Neumann
#
#' @rdname set.rng
#' @description `rvng()` returns a sequence of uniform pseudorandom numbers
#' using the Von Neumann middle-square method.
# @param n number of generations.
# @param seed seed.
#' @param k number of digits.
# @returns A numeric vector.
# @seealso [set.rng()], [rng()]
#' @export
rvng <- function(n, seed = as.numeric(Sys.time()), k = 4) {
  seed <- seed %% 10^k
  aux <- 10^(2*k-k/2)
  aux2 <- 10^(k/2)
  u <- numeric(n)
  for(i in 1:n) {
    z <- seed^2
    seed <- trunc((z - trunc(z/aux)*aux)/aux2)
    u[i] <- seed/10^k
  }
  # Almacenar semilla y parámetros
  assign(".rng", list(seed = seed, type = "vm", parameters = list(k = k)),
      envir = globalenv())
  # .rng <<- list(seed = seed, type = "vm", parameters = list(k = k))
  # Para continuar con semilla y parámetros:
  #   with(.rng, rvng(n, seed, parameters$k))
  # Devolver valores
  return(u)
}
