--- 
title: "Técnicas de Simulación y Remuestreo"
author: 
  - "Rubén Fernández Casal (ruben.fcasal@udc.es)"
  - "Ricardo Cao (rcao@udc.es)"
  - "Julián Costa (julian.costa@udc.es)"
date: "Edición: Enero de 2023. Impresión: 2025-03-05"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
lang: es
bibliography: ["packages.bib", "book.bib"]
biblio-style: apalike-es
link-citations: yes
fontsize: 10pt
github-repo: rubenfcasal/simbook
description: "Apuntes de Simulación Estadística y Remuestreo."
---

# Prólogo {-}



<!-- 
PENDENTE: 
- Código simres en capítulos
-->

Este libro, **en proceso de elaboración**, contiene los apuntes de la asignatura de [Simulación Estadística](http://eamo.usc.es/pub/mte/index.php/es/?option=com_content&view=article&id=2201&idm=13&a%C3%B1o=2019) del [Máster en Técnicas Estadísticas](http://eio.usc.es/pub/mte) y material de apoyo a la docencia de la asignatura de [Técnicas de Simulación y Remuestreo](https://guiadocente.udc.es/guia_docent/index.php?centre=614&ensenyament=614G02&assignatura=614G02036&idioma=cast) del [Grado en Ciencia e Ingeniería de Datos](https://estudos.udc.es/es/study/start/614G02V01) de la [UDC](https://www.udc.es). 

Este libro ha sido escrito en [R-Markdown](http://rmarkdown.rstudio.com) empleando el paquete [`bookdown`](https://bookdown.org/yihui/bookdown/) y está disponible en el repositorio [rubenfcasal/simbook](https://github.com/rubenfcasal/simbook) de Github. 
Este libro es una "segunda edición", la anterior (Fernández-Casal R. y Cao R., 2022, *Simulación Estadística*) está disponible en la rama *[primera_edicion](https://github.com/rubenfcasal/simbook/tree/primera_edicion)* de este repositorio. 

Se puede acceder a la versión en línea a través del siguiente enlace:

<https://rubenfcasal.github.io/simbook/index.html>.

donde puede descargarse en formato [pdf](https://rubenfcasal.github.io/simbook/Simulacion.pdf).

Para poder ejecutar los ejemplos mostrados en el libro es recomendable emplear el paquete [`simres`](https://rubenfcasal.github.io/simres), ***no disponible actualmente en CRAN***, aunque se puede instalar la versión de desarrollo en [GitHub](https://github.com/rubenfcasal/simres)):

```r
# install.packages("remotes")
# remotes::install_github("rubenfcasal/simres")
remotes::install_github("rubenfcasal/simres", INSTALL_opts = "--with-keep.source")
```
Alternativamente se pueden emplear los ficheros de la carpeta *codigo*.

Para instalar los paquetes necesarios se puede emplear los siguientes comandos:

```r
pkgs <- c('tictoc', 'boot', 'randtoolbox', 'MASS', 'DEoptim', 'nortest', 'geoR', 'copula',
          'sm', 'car', 'tseries', 'forecast', 'plot3D', 'rgl', 'rngWELL', 'randtoolbox')
install.packages(setdiff(pkgs, installed.packages()[,"Package"]), 
                 dependencies = TRUE)

# Si aparecen errores debidos a incompatibilidades entre las versiones de los paquetes, 
# probar a ejecutar en lugar de lo anterior:
# install.packages(pkgs, dependencies = TRUE) # Instala todos...
```

Para generar el libro (compilar) serán necesarios paquetes adicionales, 
para lo que se recomendaría consultar el libro de ["Escritura de libros con bookdown" ](https://rubenfcasal.github.io/bookdown_intro) en castellano.

En la Sección [Enlaces](#links) de las Referencias se incluyen recursos adicionales, incluyendo algunos que pueden ser útiles para el aprendizaje de R.

<img src="images/by-nc-nd-88x31.png" width="10%" style="display: block; margin: auto auto auto 0;" />

Este obra está bajo una licencia de [Creative Commons Reconocimiento-NoComercial-SinObraDerivada 4.0 Internacional](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.es_ES) 
(esperamos poder liberarlo bajo una licencia menos restrictiva más adelante...).


