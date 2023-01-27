# pkgdown::build_site()

#' simres: Simulation and resampling techniques
#'
#' Functions and datasets used in the books:
#' * [Técnicas de Simulación y Remuestreo](https://rubenfcasal.github.io/simbook)
#' * [Tecnicas de Remuestreo](https://rubenfcasal.github.io/book_remuestreo)
#'
#' For more information visit <https://rubenfcasal.github.io/simres.html>.
#' @keywords simulation bootstrap Monte-Carlo
#' @name simres-package
#' @aliases simres
#' @docType package
#' @import graphics
#' @import stats
#' @importFrom utils flush.console
#' @importFrom grDevices dev.interactive devAskNewPage dev.flush dev.hold
#' @references
#' Cao R., Fernández-Casal R. (2021). *[Técnicas de Remuestreo](https://rubenfcasal.github.io/book_remuestreo)*,  ([github](https://github.com/rubenfcasal/book_remuestreo)).
#'
#' Fernández-Casal R., Cao R., Costa J. (2023). *[Técnicas de Simulación y Remuestreo](https://rubenfcasal.github.io/simbook2)*, segunda edición, ([github](https://github.com/rubenfcasal/simbook2)).
#'
NULL



if(getRversion() >= "2.15.1")
    utils::globalVariables(c(".rng"))


#' Microorganism lifetimes
#'
#' Lifespan of a sample of a particular population of microorganisms.
#' @name lifetimes
#' @docType data
#' @format A vector with 15 observations of the time to death event.
#' @source Unknown.
#' @keywords datasets
#' @examples
#' hist(lifetimes)
NULL
# lifetimes <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 0.611, 0.712,
#    1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
# usethis::use_data(lifetimes)


#--------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  #--------------------------------------------------------------------
  #   pkg.info <- utils::packageDescription(pkgname, libname, fields = c("Title", "Version", "Date"))
  pkg.info <- drop( read.dcf( file = system.file("DESCRIPTION", package = "simres"),
                              fields = c("Title", "Version", "Date") ))
  packageStartupMessage(
    paste0(" simres: ", pkg.info["Title"], ",\n"),
    paste0(" version ", pkg.info["Version"], " (built on ", pkg.info["Date"], ").\n"),
    " Copyright (C) R. Fernandez-Casal 2022.\n",
    " Type `help(simres)` for an overview of the package or\n",
    " visit https://rubenfcasal.github.io/simres.\n")
}
