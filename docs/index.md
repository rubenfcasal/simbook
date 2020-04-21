--- 
title: "Simulación Estadística"
author: "Rubén Fernández Casal (ruben.fcasal@udc.es), Ricardo Cao (rcao@udc.es)"
date: "2020-04-21"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: rubenfcasal/simbook
description: "Apuntes de la asignatura de Simulación Estadística del Máster en Técnicas Estadísticas."
---

# Prólogo {-}

Este libro contiene los apuntes de la asignatura de [Simulación Estadística](http://eamo.usc.es/pub/mte/index.php/es/?option=com_content&view=article&id=2201&idm=13&a%C3%B1o=2019) del [Máster en Técnicas Estadísticas](http://eio.usc.es/pub/mte). 

Este libro ha sido escrito en [R-Markdown](http://rmarkdown.rstudio.com) empleando el paquete [`bookdown`](https://bookdown.org/yihui/bookdown/)  y está disponible en el repositorio Github: [rubenfcasal/simbook](https://github.com/rubenfcasal/simbook). 
Se puede acceder a la versión en línea a través del siguiente enlace:

<https://rubenfcasal.github.io/simbook/index.html>.

donde puede descargarse en formato [pdf](https://rubenfcasal.github.io/simbook/Simulacion.pdf).

Para instalar los paquetes necesarios para poder ejecutar los ejemplos mostrados en el libro se puede emplear el siguiente comando:

```r
pkgs <- c('boot', 'randtoolbox', 'MASS', 'DEoptim', 'nortest', 'geoR', 'copula', 'sm',
          'car', 'tseries', 'forecast', 'plot3D', 'rgl')
install.packages(setdiff(pkgs, installed.packages()[,"Package"]), 
                 dependencies = TRUE)

# Si aparecen errores debidos a incompatibilidades entre las versiones de los paquetes, 
# probar a ejecutar en lugar de lo anterior:
# install.packages(pkgs, dependencies = TRUE) # Instala todos...
```

Para generar el libro (compilar) serán necesarios paquetes adicionales, 
para lo que se recomendaría consultar el libro de ["Escritura de libros con bookdown" ](https://rubenfcasal.github.io/bookdown_intro) en castellano.



<img src="images/by-nc-nd-88x31.png" width="44" />

Este obra está bajo una licencia de [Creative Commons Reconocimiento-NoComercial-SinObraDerivada 4.0 Internacional](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.es_ES) 
(esperamos poder liberarlo bajo una licencia menos restrictiva más adelante...).


