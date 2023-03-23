# Simulación Estadística

## Fernández-Casal, R. (rfcasal@udc.es) y Cao, R. (rcao@udc.es)

***IMPORTANTE***: Esta rama se corresponde con la "primera edición" de este libro.
La **nueva edición** esta disponible en la rama principal del repositorio.

Este libro contiene los apuntes de la asignatura de [Simulación Estadística](http://eamo.usc.es/pub/mte/index.php/es/?option=com_content&view=article&id=2201&idm=13&a%C3%B1o=2019) del [Máster en Técnicas Estadísticas](http://eio.usc.es/pub/mte). 

El libro ha sido escrito en [R-Markdown](http://rmarkdown.rstudio.com) empleando el paquete [`bookdown`](https://bookdown.org/yihui/bookdown/).
y está disponible compilado en la carpeta *docs* de esta rama.

Para instalar los paquetes necesarios para poder ejecutar los ejemplos mostrados en el libro se puede emplear el siguiente comando:
```{r eval=FALSE}
pkgs <- c('boot', 'MASS', 'DEoptim', 'nortest', 'geoR', 'copula', 'sm',
          'tseries', 'forecast', 'plot3D', 'rgl')
# install.packages(pkgs, dependencies=TRUE)
install.packages(setdiff(pkgs, installed.packages()[,"Package"]), 
                 dependencies = TRUE)

# Si aparecen errores debidos a incompatibilidades entre las versiones de los paquetes, 
# probar a ejecutar en lugar de lo anterior:
# install.packages(pkgs, dependencies = TRUE) # Instala todos...
```

Para generar el libro (compilar) serán necesarios paquetes adicionales, 
para lo que se recomendaría consultar el libro de ["Escritura de libros con bookdown" ](https://rubenfcasal.github.io/bookdown_intro) en castellano.


Este obra está bajo una licencia de [Creative Commons Reconocimiento-NoComercial-SinObraDerivada 4.0 Internacional ](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.es_ES) 
(esperamos poder liberarlo bajo una licencia menos restrictiva más adelante...).

![](https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png)
