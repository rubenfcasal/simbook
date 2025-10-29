# ················································
# Herramientas para crear documentos rmarkdown
# y opciones por defecto chunks
#
# Fecha modificación: 2023/09/28
# ················································
# NOTA: Ctrl + Shift + O para Document Outline

if (!require(simres))
  remotes::install_github("rubenfcasal/simres", INSTALL_opts = "--with-keep.source")

# Funciones auxiliares rmarkdown -----------------
# ················································

# Output bookdown
# ················································
is_latex <- function(...) knitr:::is_latex_output(...)
is_html <- function(...) knitr:::is_html_output(...)

# Rmd code:
# ················································
inline <- function(x = "") paste0("`` `r ", x, "` ``")
inline2 <- function(x = "") paste0("`r ", x, "`")


owidth <- 75 # 70
fowidth <- function(d) sprintf("%i%%", owidth - d)
rowidth <- function(p) sprintf("%i%%", round(0.01*p*owidth))

# Opciones knitr
# ················································

# Opciones salida de texto
options(max.print = 60, width = 73, # 67, # (por defecto 80)
        str = strOptions(strict.width = "cut"), # str()
        digits = 5)

# Opciones gráficas
knitr::knit_hooks$set(small.mar = function(before, ...) {
  if (before){
    par(mar = c(bottom = 4, left = 4, top = 1, right = 1) + 0.1)
  } else
    par(mar = c(bottom = 5, left = 4, top = 4, right = 2) + 0.1)

})


# Establecer opciones chunks
knitr::opts_chunk$set(
  fig.dim = c(7, 5), fig.align = "center", fig.pos = "!htbp", # fig.pos = "!htb"
  out.width = fowidth(0), # out.lines = 50,
  cache = FALSE, cache.path = 'cache/',
  echo = TRUE, warning = FALSE, message = FALSE,
  comment = if(is_html()) "##" else " ##",
  fig.show = if(is_html()) "asis" else "hold",
  small.mar = !is_html() # , crop = TRUE, # Recortar figuras
)


# Directorio figuras
# ················································
fig.path <- "images/"
# fig.path <- ""

.regerar <- FALSE


# Citas Figuras ----------------------------------
# ················································

# (ver Figura/figuras ...)
cite_fig <- function(..., text = "(ver ") {
    x <- as.character(rlang::enexprs(...)) # as.character(rlang::ensyms(...)) # list(...)
    x <- gsub(" - ", "-", x) # Evita que a-b se traduzca en a - b
    paste0(text, if(length(x)>1) "figuras: " else "Figura ",
      paste0("\\@ref(fig:", x, ")", collapse = ", "),
    ")")
}

# [Figura/s: ...]
cite_fig2 <- function(..., text = "") {
    x <- as.character(rlang::enexprs(...)) # as.character(rlang::ensyms(...)) # list(...)
    x <- gsub(" - ", "-", x) # Evita que a-b se traduzca en a - b
    paste0(text, if(length(x)>1) "[Figuras: " else "[Figura ",
      paste0("\\@ref(fig:", x, ")", collapse = ", "),
    "]")
}

# (ver Figura/figuras ...) si latex/pdf
latexfig <- function(..., output = is_latex())
    if (output) cite_fig(..., text = " ") else ""

# [Figura/s: ...] si latex/pdf
latexfig2 <- function(..., output = is_latex())
    if (output) cite_fig2(..., text = " ") else ""

# cite_fig(fig1)
# "(ver Figura \\@ref(fig:fig1))"
# cite_fig("fig1")
# "(ver Figura \\@ref(fig:fig1))"
# cite_fig(fig1, fig2)
# "(ver figuras: \\@ref(fig:fig1), \\@ref(fig:fig2))"
# cite_fig("fig1", "fig2")
# "(ver figuras: \\@ref(fig:fig1), \\@ref(fig:fig2))"


# Citas paquetes y funciones ---------------------
# ················································

## Citas paquetes --------------------------------
# ················································


# Cita paquete CRAN
cite_cran <- function(pkg) {
    pkg <- as.character(substitute(pkg))
    paste0("[`", pkg, "`](https://CRAN.R-project.org/package=", pkg, ")")
}


# https://rubenfcasal.github.io/

cite_github <- function(pkg = simres){
  pkg <- as.character(substitute(pkg))
  cite_pkg_(pkg, paste0("https://rubenfcasal.github.io/", pkg))
}

cite_fgithub <- function(fun, pkg = simres, lnk = NULL, full = FALSE)  {
  fun <- as.character(substitute(fun))
  pkg <- as.character(substitute(pkg))
  paste0(if(full) paste0("[`", pkg, "::") else "[`", fun,
         paste0("()`](https://rubenfcasal.github.io/", pkg, "/reference/"),
         if (!is.null(lnk)) lnk else fun,
         ".html)",
         collapse = ", ")
  # downlit::autolink_url
}


# Pendiente: múltiples paquetes

url_pkg <- function(pkg){
  url <- downlit::href_package(pkg)
  if (is.na(url)) url <- paste0("https://CRAN.R-project.org/package=", pkg)
  url
}

cite_pkg_ <- function(pkg, url = sapply(pkg, url_pkg)) {
    paste0("[`", pkg, "`](", url, ")",  collapse = ", ")
}

# cite_pkg_(c("dplyr", "tidyr"))
# "[`dplyr`](https://dplyr.tidyverse.org), [`tidyr`](https://tidyr.tidyverse.org)"
# cite_pkg_(c("dplyr", "tidyr"), c("https://dplyr.tidyverse.org", "https://tidyr.tidyverse.org"))
# "[`dplyr`](https://dplyr.tidyverse.org), [`tidyr`](https://tidyr.tidyverse.org)"
# Cuidado: cite_pkg_("Rcmdr")
# "[`Rcmdr`](https://www.r-project.org)" Error descripción Rcmdr
# path <- system.file("DESCRIPTION", package = "Rcmdr")
# read.dcf(path, fields = "URL")[[1]]
# ················································

cite_pkg <- function(...) {
    cite_pkg_(as.character(rlang::ensyms(...)))
}

# cite_pkg(dplyr)
# "[`dplyr`](https://dplyr.tidyverse.org)"
# cite_pkg(dplyr, tidyr)
# "[`dplyr`](https://dplyr.tidyverse.org), [`tidyr`](https://tidyr.tidyverse.org)"
# ················································


## Citas funciones -------------------------------
# ················································

cite_fun_ <- function(fun, pkg, url, full = FALSE) {
    fun_full <- if (!missing(pkg))
        paste(pkg, fun, sep = "::") else fun
    if (missing(url)) url <- downlit::autolink_url(fun_full)
    if (full) fun <- fun_full
    paste0("[`", fun, "`](", url, ")", collapse = ", ")
}

# cite_fun_("subset()")
# "[`subset()`](https://rdrr.io/r/base/subset.html)"
# cite_fun_("group_by()", "dplyr")
# "[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)"
# Pendente: cite_fun_(c("summarise", "group_by"), "dplyr")
# ················································

cite_fun <- function(fun, pkg, ...) {
    fun <- paste0(deparse(substitute(fun)),"()")
    if (missing(pkg)) return(cite_fun_(fun, ...))
    cite_fun_(fun, deparse(substitute(pkg)), ...)

}

cite_funs <- function(..., pkg, url, full = FALSE) {
    fun <- paste0(as.character(rlang::ensyms(...)), "()")
    fun_full <- if (!missing(pkg))
        paste(deparse(substitute(pkg)), fun, sep = "::") else fun
    if (missing(url)) url <- sapply(fun_full, downlit::autolink_url)
    if (full) fun <- fun_full
    paste0("[`", fun, "`](", url, ")", collapse = ", ")
    # Pendiente: " y " en lugar de ", " antes de última...
}


cite_method <- function(fun, class, pkg, url, full = FALSE) {
    fun <- deparse(substitute(fun))
    fun_full <- if (!missing(pkg))
        paste(deparse(substitute(pkg)), fun, sep = "::") else fun
    fun_full <-  paste0(fun_full, ".", deparse(substitute(class)), "()")
    if (missing(url)) url <- downlit::autolink_url(fun_full)
    fun <- if (full) fun_full else paste0(fun, "()")
    paste0("[`", fun, "`](", url, ")", collapse = ", ")
}

# cite_fun(subset)
# "[`subset()`](https://rdrr.io/r/base/subset.html)"
# cite_fun(group_by, dplyr)
# "[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)"
# cite_fun(dplyr::group_by)
# cite_fun(group_by, dplyr, full = TRUE)
# "[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)"
# cite_funs(group_by, filter, pkg = dplyr)
# "[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html), [`filter()`](https://dplyr.tidyverse.org/reference/filter.html)"
# Pendente:
#   cite_funs(dplyr::group_by) # No lo admite rlang::ensyms()
# ················································


# Recortar figuras
# knitr::knit_hooks$set(crop = knitr::hook_pdfcrop)
# # https://ghostscript.com/releases/gsdnld.html
# Sys.setenv(R_GSCMD="C:/Program Files/gs/gs10.03.0/bin/gswin64.exe")

# PENDENTE: ----
# ················································
# "[ Se ha omitido parte de la salida de texto... ]\n"
# rmd.lines <- function(l = 1) paste0("<br> \vspace{0.5cm}\n")
#     cat(rep("<br>", l), "\n") # 0.5*l


# https://stackoverflow.com/questions/35317587/extract-names-of-dataframes-passed-with-dots
# H. Wickham, Advanced R, [Section 19.3.2](https://adv-r.hadley.nz/quasiquotation.html#capturing-symbols)).
names_from_dots <- function(...) as.character(rlang::ensyms(...))
names_from_dots(swiss, iris)
# [1] "swiss" "iris"


