# Este script genera el libro gitbook y utiliza el paquete downlit para
# añadir automáticamente enlaces a la documentación de las funciones
# install.packages(c("downlit", "xml2"))

# Establecer opciones paquete downlit
# Paquetes cargados en la sesión (enlazar con web en lugar de rdrr.io)
# pks <- c("sf", "sp", "stars", "gstat", "geoR", "sftime", "npsp",
#          "ggplot2", "dplyr", "tmap", "mapview",
#          "osmdata", "rnaturalearth", "mapSpain", "meteospain")
# options(downlit.attached = pks)

# Se supone que el paquete está en el directorio actual
if(!file.exists("_bookdown.yml")) stop()
message(rstudioapi::getActiveProject())

# Generar el libro
bookdown::render_book("index.Rmd", "bookdown::gitbook")

# Obtener directorio de salida
config <- yaml::read_yaml("_bookdown.yml")
output_dir <- config$output_dir
if (is.null(output_dir) || !nzchar(output_dir)) output_dir <- "_book"

# Listado archivos HTML
files <- list.files(output_dir, pattern = "\\.html$", full.names = TRUE)
# files <- files[48:56]

# Añadir enlaces
for (f in files) { # f <- "docs/operaciones-datos.html"
  message("Adding documentation links: ", f)
  html <- xml2::read_html(f)
  if (inherits(html, "xml_node")) {
    downlit::downlit_html_node(html)
    xml2::write_html(html, f)
  }
}
