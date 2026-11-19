#' ---
#' title: "Utilidades Bookdown"
#' author: "Rubén Fernández Casal"
#' date: "octubre 2026"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: yes
#'     toc_depth: 3
#'   bookdown::pdf_document2: default
#' ---
#'
#' # Resumen
#' <!--
# NOTE: Press Ctrl + Shift + O to show document outline in RStudio
# Resumen -----------------------------------------------
#' -->
#'
#' Este fichero contiene funciones auxiliares para la escritura
#' de documentos con [`knitr`](https://yihui.org/knitr/) o
#' [`bookdown`](https://github.com/rstudio/bookdown), con la idea
#' de facilitar la documentación de las herramientas de R en libros y apuntes.
#'
#' Se puede incluir al principio del fichero un chunk de la forma:
#'
#' ````markdown
#' ````r '{r global-options, include=FALSE}'`
#' source("_global_options.R")
#' ```
#' ````
#'
#' Las herramientas principales serían:
#'
#' - Utilidades knitr:
#'    * Detectar el formato de salida (HTML o PDF): `is_latex()`, `is_html()`.
#'    * Insertar expresiones R en línea: `inline()`
#'    * Configurar opciones globales: `owidth`, `fowidth()`, `fig.path`,...
#' - Crear referencias cruzadas a figuras:
#'   `cite_fig()`, `cite_fig2()`
#' - Generar enlaces a paquetes y funciones de R:
#'    * `cite_cran()`, `cite_pkg()`, `cite_github()`, `cite_fgithub()`
#'    * `cite_fun()`, `cite_funs()`, `cite_method()`
#'
#' Se utilizan los paquetes `rlang`, `downlit` y `knitr`.
#'
#' Para más información copilar este documento [spin](https://rdrr.io/cran/knitr/man/spin.html),
#' pulsando *Ctrl + Shift + K*, el icono correspondiente en la barra superior o seleccionando
#' *File > Compile Report...* en RStudio (o ejecutando
#' `browseURL(url = rmarkdown::render(knitr::spin("_global_options.R", knit = FALSE)))`).
#'
#' # Utilidades `knitr`
#' <!--
# Utilidades `knitr` -----------------------------------------------
#' -->
#'
#' ## Detectar el formato de salida
#'
#' Devuelven `TRUE` dependiendo del formato de salida del documento.

is_latex <- function(...) knitr::is_latex_output(...)
is_html  <- function(...) knitr::is_html_output(...)

#'
#' ## Código R en línea
#'
#' Permiten mostrar el texto rmarkdown necesario para incluir código R en línea
#' (dentro de un párrafo).

inline <- function(x = "") paste0("`` `r ", x, "` ``")
inline2 <- function(x = "") paste0("`r ", x, "`")

#'
#' ## Anchura relativa de figuras
#'
#' `owidth` es el porcentaje por defecto de anchura de figuras respecto al
#' ancho total, opción [`out.width`](https://yihui.org/knitr/options/#plots) del chunk.

owidth <- if(is_html()) 80 else 75 # 70

#'
#' `fowidth()` y `rowidth()` permiten establecer la opción `out.width` de un
#' chunk a partir del valor por defecto, como reducción o valor relativo.

fowidth <- function(d) sprintf("%i%%", owidth - d)
rowidth <- function(p) sprintf("%i%%", round(p*owidth))

#'
#' ## Añadir opciones chunks
#'
#' `small.mar = TRUE` reduce los márgenes por defecto.

knitr::knit_hooks$set(small.mar = function(before, ...) {
  if (before){
    par(mar = c(bottom = 4, left = 4, top = 1, right = 1) + 0.1)
  } else
    par(mar = c(bottom = 5, left = 4, top = 4, right = 2) + 0.1)

})

# Recortar figuras (`crop = TRUE`).
# knitr::knit_hooks$set(crop = knitr::hook_pdfcrop)
# # https://ghostscript.com/releases/gsdnld.html
# Sys.setenv(R_GSCMD="C:/Program Files/gs/gs10.03.0/bin/gswin64.exe")

#'
#' # Opciones globales
#' <!--
# Opciones globales -----------------------------------------------
#' -->
#'
#' Las siguientes instrucciones configuran el comportamiento por defecto de R y
#' knitr.
#'
#' ## Opciones salida de texto

options(max.print = 60,
        width = 73, # 67, # (por defecto 80)
        str = strOptions(strict.width = "cut"), # str()
        digits = 5
)

#'
#' ## Opciones por defecto chunks

knitr::opts_chunk$set(
  fig.dim = c(7, 5), fig.align = "center",
  fig.pos = "!htbp", # fig.pos = "!htb"
  out.width = fowidth(0), # out.lines = 50,
  cache = FALSE, cache.path = 'cache/', # Cache
  echo = TRUE, warning = FALSE, message = FALSE, # Opciones mostrar
  comment = if(is_html()) "##" else " ##", # Prefijo resultados texto
  fig.show = if(is_html()) "asis" else "hold", # Incluir gráficos al generar/final
  small.mar = !is_html() # Reducir márgenes
  # , crop = TRUE, # Recortar figuras
)

# NOTA: la opción `fig.dim` no es compatible con `fig.width` y `fig.height`
# (usar solo `fig.dim` en las opciones de los chunks).

# Directorio figuras estáticas
fig.path <- "images/"
# fig.path <- ""

# Opción para generar figuras estáticas
.regenerar <- FALSE


#'
#' # Referencias a figuras
#' <!--
# Citas Figuras -----------------------------------------------
#' -->
#'
#' Estas funciones generan referencias compatibles con **bookdown**
#' (admiten evaluación no estándar).
#'

# (ver Figura/figuras ...)
cite_fig <- function(..., text = "(ver ") {
    x <- as.character(rlang::enexprs(...)) # as.character(rlang::ensyms(...)) # list(...)
    x <- gsub(" - ", "-", x) # Evita que a-b se traduzca en a - b
    paste0(text, if(length(x)>1) "figuras: " else "Figura ",
      paste0("\\@ref(fig:", x, ")", collapse = ", "),
    ")")
}

# (Figura/s: ...)
cite_fig2 <- function(..., text = "(") {
    x <- as.character(rlang::enexprs(...)) # as.character(rlang::ensyms(...)) # list(...)
    x <- gsub(" - ", "-", x) # Evita que a-b se traduzca en a - b
    paste0(text, if(length(x)>1) "Figuras: " else "[Figura ",
      paste0("\\@ref(fig:", x, ")", collapse = ", "),
    ")")
}

# _(ver Figura/figuras ...) si latex/pdf
latexfig <- function(..., text = " (ver ", output = is_latex())
    if (output) cite_fig(..., text = text) else ""

# _(Figura/s: ...) si latex/pdf
latexfig2 <- function(..., text = " (", output = is_latex())
    if (output) cite_fig2(..., text = text) else ""

#' Una figura en un chunk con nombre `chunkname` puede citarse posteriormente
#' como `r inline("cite_fig(chunkname)")`.
#' Ejemplos:
#'
#' ```{r}
#' cite_fig(fig1)
#' cite_fig("fig1")
#' cite_fig(fig1, fig2)
#' cite_fig("fig1", "fig2")
#' ```
#'
#' La idea de las funciones para latex es no interferir cuando se utilizan
#' en html (se pueden poner entre la última palabra de una frase y el punto).
#'
#' # Referencias a paquetes y funciones
#' <!--
# Citas paquetes y funciones ---------------------
#' -->
#'
#' NOTA: Si termina en `_()` los argumentos tiene que ser de texto, en caso
#' contrario utiliza evaluación no estándar.
#'
#' ## Citas paquetes
#' <!--
## Citas paquetes --------------------------------
#' -->
#'
#' Estas funciones crean enlaces automáticos a la web (oficial) del paquete.
#'
#' Citar paquete CRAN:

cite_cran <- function(pkg) {
    pkg <- as.character(substitute(pkg))
    paste0("[`", pkg, "`](https://CRAN.R-project.org/package=", pkg, ")")
}

#' Citar paquete o función en <https://rubenfcasal.github.io>:

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

#' Pendiente: emplear herramientas para definir enlaces "locales":
#' <https://downlit.r-lib.org/reference/autolink.html>.
#'
#' Citar paquete, por defecto según referencia en descripción:

# Obtener la referencia a un paquete según descripción o a la web del CRAN
# si no se encuentra
# Cuidado: algunos paquete pueden tener una url no "válida" en la descripción
# path <- system.file("DESCRIPTION", package = "Rcmdr") # Rcmdr ya la corrigió
# read.dcf(path, fields = "URL")[[1]]
# Pendiente: múltiples paquetes
url_pkg <- function(pkg){
  url <- downlit::href_package(pkg)
  if (is.na(url)) url <- paste0("https://CRAN.R-project.org/package=", pkg)
  url
}

cite_pkg_ <- function(pkg, url = NULL) {
  if (is.null(url)) url <- sapply(pkg, url_pkg)
  paste0("[`", pkg, "`](", url, ")",  collapse = ", ")
}


cite_pkg <- function(..., url = NULL) {
  cite_pkg_(as.character(rlang::ensyms(...), url = url))
}

#' Ejemplos:
#'
#' ```{r}
#' cite_pkg_(c("dplyr", "tidyr"))
#' cite_pkg_(c("dplyr", "tidyr"), c("https://dplyr.tidyverse.org", "https://tidyr.tidyverse.org"))
#' cite_pkg(dplyr)
#' cite_pkg(dplyr, tidyr)
#' ```
#'
#' ## Citas funciones
#' <!--
## Citas funciones --------------------------------
#' -->
#'
#' Generan enlaces hacia la documentación de funciones de R
#' (emplea <https://rdrr.io/r/> si el paquete no tiene una web "adecuadamente"
#' generada con [`pkgdown`](https://pkgdown.r-lib.org/articles/linking.html)).
#'
#' `cite_fun_()` usa evaluación estándar y requiere que el argumento `fun`
#' sea un texto con el nombre de la función terminado en paréntesis, y puede
#' ser una referencia completa `"pkg::fun()"`.
#' Si `full = TRUE` y `!missing(pkg)` se emplea la referencia completa.

cite_fun_ <- function(fun, pkg, url, full = FALSE) {
    fun_full <- if (!missing(pkg))
        paste(pkg, fun, sep = "::") else fun
    if (missing(url)) url <- downlit::autolink_url(fun_full)
    if (full) fun <- fun_full
    paste0("[`", fun, "`](", url, ")", collapse = ", ")
}

# Pendente: cite_fun_(c("summarise", "group_by"), "dplyr")

#' Ejemplos:
#'
#' ```{r}
#' cite_fun_("subset()")
#' cite_fun_("group_by()", "dplyr")
#' ```

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

#' Ejemplos:
#'
#' ```{r}
#' cite_fun(subset)
#' cite_fun(group_by, dplyr)
#' cite_fun(group_by, dplyr, full = TRUE)
#' cite_funs(group_by, filter, pkg = dplyr)
#' cite_method(plot, sf, sf)
#' cite_method(sf::plot, sf) # cite_method(plot, sf, sf, full = TRUE)
#' ```

# Pendente:
# cite_funs(dplyr::group_by) # No lo admite rlang::ensyms()
# https://stackoverflow.com/questions/35317587/extract-names-of-dataframes-passed-with-dots
# H. Wickham, Advanced R, [Section 19.3.2](https://adv-r.hadley.nz/quasiquotation.html#capturing-symbols)).
# names_from_dots <- function(...) as.character(rlang::ensyms(...))
# names_from_dots(swiss, iris)
# [1] "swiss" "iris"


# PENDENTE: ----
# ················································
# Opción chunk recortar líneas de texto (`out.lines = 50`)
# "[ Se ha omitido parte de la salida de texto... ]\n"
# rmd.lines <- function(l = 1) paste0("<br> \vspace{0.5cm}\n")
#     cat(rep("<br>", l), "\n") # 0.5*l




