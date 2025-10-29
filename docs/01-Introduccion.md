# Introducción a la simulación {#intro}

<!-- Capítulo \@ref(intro) -->



<!-- 
---
title: "Introducción a la simulación"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    pandoc_args: ["--number-offset", "0,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("01-Introduccion.Rmd")
knitr::purl("01-Introduccion.Rmd", documentation = 2)
knitr::spin("01-Introduccion.R",knit = FALSE)

PENDENTE:
  - cite_fsimres
-->

Cuando pensamos en ciencia pensamos en experimentos y en modelos. 
Se experimenta una y otra vez sobre el fenómeno real que se desea conocer mejor para, con la información así acumulada, construir un modelo teórico, que no es sino una representación simplificada (más o menos acertada) del fenómeno real.
Como el modelo se formula en términos matemáticos, en general es susceptible de un estudio analítico del que poder sacar conclusiones.

La simulación ofrece una alternativa a esa última fase del proceso, y sustituye (en parte o completamente) el estudio analítico por más experimentación, pero esta vez sobre el propio modelo en lugar de sobre la realidad.

Así, se puede definir la *simulación* como una técnica que consiste en realizar experimentos sobre el modelo de un sistema (experimentos de muestreo si la simulación incorpora aleatoriedad), con el objetivo de recopilar información bajo determinadas condiciones. 

## Conceptos básicos {#conceptos}

La experimentación directa sobre la realidad puede tener muchos inconvenientes, entre otros:

-   Coste elevado: por ejemplo cuando las pruebas son destructivas o si es necesario esperar mucho tiempo para observar los resultados.
    
-   Puede no ser ética: por ejemplo la experimentación sobre seres humanos o la dispersión de un contaminante.

-   Puede resultar imposible: por ejemplo cuando se trata de un acontecimiento futuro o una alternativa en el pasado.

Además la realidad puede ser demasiado compleja como para ser estudiada directamente y resultar preferible trabajar con un modelo del sistema real. 
Un modelo no es más que un conjunto de variables junto con ecuaciones matemáticas que las relacionan y restricciones sobre dichas variables. 
Habría dos tipos de modelos:

-   Modelos deterministas: en los que bajo las mismas condiciones (fijados los valores de las variables explicativas) se obtienen siempre los mismos resultados.

-   Modelos estocásticos (con componente aleatoria): tienen en cuenta la incertidumbre asociada al modelo. Tradicionalmente se supone que esta incertidumbre es debida a que no se dispone de toda la información sobre las variables que influyen en el fenómeno en estudio (puede ser debida simplemente a que haya errores de medida), lo que se conoce como *aleatoriedad aparente*:

    > "Nothing in Nature is random... a thing appears random only through the incompleteness of our knowledge."
    >
    > --- @spinoza1667ethics
    
    aunque hoy en día gana peso la idea de la física cuántica de que en el fondo hay una *aleatoriedad intrínseca*^[Como ejemplo, en física cuántica, la [Ecuación de Schrödinger](https://es.wikipedia.org/wiki/Ecuaci%C3%B3n_de_Schr%C3%B6dinger) es un modelo determinista que describe la evolución en el tiempo de la función de onda de un sistema. Sin embargo, como las funciones de onda pueden cambiar de forma aleatoria al realizar una medición, se emplea la [regla de Born](https://es.wikipedia.org/wiki/Regla_de_Born) para modelar las probabilidades de las distintas posibilidades (algo que inicialmente generó rechazo, dió lugar a la famosa frase de Einstein "Dios no juega a los dados", pero experimentos posteriores parecen confirmar). Por tanto en la práctica se emplea un modelo estocástico.].

La modelización es una etapa presente en la mayor parte de los trabajos de investigación, especialmente en las ciencias experimentales.
El modelo debería considerar las variables más relevantes para explicar el fenómeno en estudio y las principales relaciones entre ellas.
La inferencia estadística proporciona herramientas para estimar los parámetros y contrastar la validez de un modelo estocástico a partir de los datos observados.

La idea es emplear el modelo, asumiendo que es válido, para resolver el problema de interés. 
Si se puede obtener la solución de forma analítica, esta suele ser exacta (aunque en ocasiones solo se dispone de soluciones aproximadas, basadas en resultados asintóticos, o que dependen de suposiciones que pueden ser cuestionables) y a menudo la resolución también es rápida.
Cuando la solución no se puede obtener de modo analítico (o si la aproximación disponible no es adecuada) se puede recurrir a la simulación.
De esta forma se pueden obtener resultados para un conjunto más amplio de modelos, que pueden ser mucho más complejos.

Nos centraremos en el caso de la *simulación estocástica*: las conclusiones se obtienen generando repetidamente simulaciones del modelo aleatorio.
Muchas veces se emplea la denominación de *método Monte Carlo*[^01-introduccion-1] como sinónimo de simulación estocástica, pero realmente se trata de métodos especializados que emplean simulación para resolver problemas que pueden no estar relacionados con un modelo estocástico de un sistema real. Por ejemplo, en el Capítulo \@ref(monte-carlo) se tratarán métodos de integración y optimización Monte Carlo.

[^01-introduccion-1]: Estos métodos surgieron a finales de la década de 1940 como resultado del trabajo realizado por Stanislaw Ulam y John von Neumann en el proyecto Manhattan para el desarrollo de la bomba atómica. Al parecer, como se trataba de una investigación secreta, Nicholas Metropolis sugirió emplear el nombre clave de "Monte-Carlo" en referencia al casino de Monte Carlo de Mónaco.
Esta denominación se empleó por primera vez en una publicación en @metropolis1949.

<!-- 
Ejemplo: caballero de Meré 
Experimentación directa sobre la realidad
Modelo de probabilidad
Aproximación por simulación
-->


### Ejemplo {#ealbum}

Supongamos que nos regalan un álbum con $n = 75$ cromos, que se venden sobres con $m = 6$ cromos por 0.8€, y que estamos interesados en el número de sobres que hay que comprar para completar la colección, por ejemplo en su valor medio.

Podemos aproximar la distribución del número de sobres para completar la colección a partir de $nsim=1000$ simulaciones de coleccionistas de cromos:


``` r
# Parámetros
n <- 75 # Número total de cromos
m <- 6  # Número de cromos en cada sobre
repe <- TRUE # Repetición de cromos en cada sobre
# Número de simulaciones
nsim <- 1000
# Resultados simulación
nsobres <- numeric(nsim)
# evol <- vector("list", nsim)
# Fijar semilla
set.seed(1)
# Bucle simulación
for (isim in 1:nsim) {
  # seed <- .Random.seed    # .Random.seed <- seed
  album <- logical(n)
  i <- 0 # Número de sobres
  while(sum(album) < n) {
    i <- i + 1
    album[sample(n,m, replace = repe)] <- TRUE
  }
  nsobres[isim] <- i
}
```

Distribución del número de sobres para completar la colección
(aproximada por simulación):


``` r
hist(nsobres, breaks = "FD", freq = FALSE,
     main = "", xlab = "Número de sobres")
lines(density(nsobres))
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{01-Introduccion_files/figure-latex/distr-ealbum-1} 

}

\caption{Aproximación por simulación de la distribución del número de sobres para completar la colección.}(\#fig:distr-ealbum)
\end{figure}

Aproximación por simulación del número medio de sobres para completar la colección:


``` r
sol <- mean(nsobres)
sol
```

```
 ## [1] 61.775
```

<!-- 
Análisis de la convergencia:


``` r
plot(1:nsim, cumsum(nsobres)/1:nsim, type = "l",
     ylab="Número de sobres", xlab="Número de simulaciones")
abline(h = sol)
```



\begin{center}\includegraphics[width=0.75\linewidth]{01-Introduccion_files/figure-latex/unnamed-chunk-3-1} \end{center}
-->

Número mínimo de sobres para asegurar de que se completa la colección con una probabilidad del 95\%:


``` r
nmin <- quantile(nsobres, probs = 0.95)
ceiling(nmin)
```

```
 ## 95% 
 ##  92
```

``` r
# Reserva de dinero para poder completar la colección el 95% de las veces:
ceiling(nmin)*0.8
```

```
 ##  95% 
 ## 73.6
```

``` r
hist(nsobres, breaks = "FD", freq = FALSE,
     main = "", xlab = "Número de sobres")
lines(density(nsobres))
abline(v = sol)
abline(v = nmin, lty = 2)
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{01-Introduccion_files/figure-latex/distr-ealbum2-1} 

}

\caption{Aproximaciones por simulación de la distribución del número de sobres para completar la colección, de su valor esperado (línea vertical continua) y del cuantil 0.95 (línea vertical discontinua).}(\#fig:distr-ealbum2)
\end{figure}

Por supuesto, la distribución del gasto necesario para completar la colección es esta misma reescalada.

(ref:mc-plot-ealbum) Gráficos exploratorios de las simulaciones del gasto para completar la colección obtenidos con la función `simres::mc.plot()`.


``` r
res <- simres::mc.plot(nsobres*0.8)
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=1\linewidth]{01-Introduccion_files/figure-latex/mc-plot-ealbum-1} 

}

\caption{(ref:mc-plot-ealbum)}(\#fig:mc-plot-ealbum)
\end{figure}

Aproximación del gasto medio:


``` r
res$approx  # sol*0.8
```

```
 ## [1] 49.42
```
 
En el Ejercicio \@ref(exr:album) se propone modificar este código para obtener información adicional sobre la evolución del número de cromos distintos dependiendo de los sobres comprados por un coleccionista.


### Ventajas e inconvenientes de la simulación 

A modo de resumen, la simulación presenta ventajas respecto a la experimentación real o a la resolución analítica en las siguientes situaciones [que coinciden en gran medida con las ya descritas por @shannon1975systems]:

-   Cuando la resolución analítica no puede llevarse a cabo (todavía no se han desarrollado métodos analíticos para resolver el modelo matemático, como es el caso de algunos modelos de colas).

-   Cuando existen medios para resolver analíticamente el problema, pero dicha resolución es demasiado complicada o ardua (comparada con la simplicidad de la solución mediante simulación), o asumen simplificaciones (hipótesis estructurales) que no son muy realistas (y en el mejor de los casos proporcionarían una solución aproximada).

-   Si resulta imposible experimentar sobre el sistema real (en ese momento no existe), resulta demasiado complicado (e.g. misiones espaciales), es demasiado costoso (e.g. requieren pruebas destructivas), o no es ético (e.g. dispersión de contaminantes o pruebas nucleares).

-   En sistemas que evolucionan muy lentamente en el tiempo (la simulación permitiría observar la evolución del proceso acelerando o ralentizando el tiempo).


El principal inconveniente puede ser el tiempo de computación necesario, aunque gracias a la gran potencia de cálculo de los computadores actuales, se puede obtener rápidamente una solución aproximada en la mayor parte de los problemas susceptibles de ser modelizados.
Además siempre están presentes los posibles problemas debidos a emplear un modelo:

-   La construcción de un buen modelo puede ser una tarea muy costosa 
    (compleja, laboriosa y requerir mucho tiempo; 
    e.g. modelos climáticos).

-   Frecuentemente el modelo omite variables o relaciones importantes entre ellas
    (los resultados pueden no ser válidos para el sistema real).

-   Resulta difícil conocer la precisión del modelo formulado.

Otro problema de la simulación es que se obtienen resultados para unos valores concretos de los parámetros del modelo, por lo que en principio
resultaría complicado extrapolar las conclusiones a otras situaciones.


### Aplicaciones de la simulación

La simulación resulta de utilidad en multitud de contextos diferentes.
Los principales campos de aplicación son:

-   Estadística:

    -   Muestreo, remuestreo, *bagging*...
    
    -   Aproximación de distribuciones (de estadísticos, estimadores...)
    
    -   Realización de contrastes, intervalos de confianza...
    
    -   Comparación de estimadores, contrastes...
    
    -   Validación teoría (distribución asintótica...)
    
    -   Inferencia Bayesiana

-   Optimización: Algoritmos genéticos, temple simulado...

-   Análisis numérico: Aproximación de integrales, resolución de ecuaciones...

-   Computación: Diseño, verificación y validación de algoritmos...

-   Criptografía: Protocolos de comunicación segura...

-   Física: Simulación de fenómenos naturales...


En el Capítulo \@ref(monte-carlo) se describen las principales aplicaciones de uso general: integración Monte Carlo (secciones \@ref(int-MC) y \@ref(muestreo-importancia)) y optimización Monte Carlo (Sección \@ref(opt-MC)). 
En la Sección \@ref(aplic-inf), y a partir del Capítulo \@ref(bootstrap) (donde se introducen los métodos de remuestreo bootstrap), nos centraremos en las más empleadas en Estadística.


<!-- 
Muchas veces se utiliza la simulación como primer paso en el desarrollo de métodos de inferencia (es un método sencillo para comenzar a estudiar si un nuevo procedimiento puede resultar de interés para resolver un problema, posteriormente puede ser estudiado analíticamente si aparentemente resulta ser útil).
También pud

-->


## Tipos de números aleatorios

El primer requisito para poder realizar simulación estocástica sería disponer de números aleatorios.
Se distingue entre tres tipos de secuencias:

- *números aleatorios puros* (*true random*): se caracteriza porque no existe ninguna regla o plan que nos permita conocer sus valores.

- *números pseudo-aleatorios*: imitan realizaciones de una variable aleatoria (uniforme).

- *números cuasi-aleatorios*: secuencias deterministas con una distribución más regular en el rango considerado.


### Números aleatorios puros

Normalmente son obtenidos por procesos físicos (loterías, ruletas, ruidos...) y, hasta hace una décadas, se almacenaban en *tablas de dígitos aleatorios*. 
Por ejemplo, en 1955 la Corporación RAND publicó el libro [*A Million Random Digits with 100,000 Normal Deviates*](https://www.rand.org/pubs/monograph_reports/MR1418.html) que contenía números aleatorios generados mediante una ruleta electrónica conectada a una computadora (ver Figura \@ref(fig:randbook)).

(ref:randbook) Líneas 10580-10594, columnas 21-40, del libro *A Million Random Digits with 100,000 Normal Deviates* (Fuente: [Wikimedia Commons](https://es.wikipedia.org/wiki/Un_mill%C3%B3n_de_d%C3%ADgitos_aleatorios_con_100.000_desviaciones_normales)).

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.3\linewidth]{images/RAND} 

}

\caption{(ref:randbook)}(\#fig:randbook)
\end{figure}

El procedimiento que se utilizaba para seleccionar de una tabla, de forma manual, números aleatorios en un rango de $1$ a $m$ era el siguiente:

-   Se selecciona al azar un punto de inicio en la tabla 
    y la dirección que se seguirá.

-   Se agrupan los dígitos de forma que “cubran” el valor de $m$.

-   Se va avanzado en la dirección elegida, seleccionando los valores menores o iguales que $m$ y descartando el resto.

Hoy en día están disponibles generadores de números aleatorios "online". 
Por ejemplo, el dominio [RANDOM.ORG](http://www.random.org/integers) proporciona números "verdaderamente" aleatorios obtenidos a partir de ruido atmosférico (ver paquete [`random`](https://CRAN.R-project.org/package=random) de R).
Aunque para un uso profesional es recomendable emplear generadores implementados mediante hardware (ver e.g. [Wikipedia: Hardware random number generator](https://en.wikipedia.org/wiki/Hardware_random_number_generator)).
Los procesadores recientes de Intel o de AMD incorporan esta funcionalidad, mediante la instrucción `RDSEED` integrada en la CPU[^01-introduccion-2].
Sus principales aplicaciones hoy en día son en criptografía y juegos de azar, donde resulta especialmente importante su impredecibilidad.

[^01-introduccion-2]: Para un uso intensivo, la recomendación es utilizar la instrucción `RDRAND` que emplea un generador de números pseudoaleatorios alimentado por el generador de números aleatorios "verdaderos" (tomando periódicamente la semilla de `RDSEED`). Para más detalles ver, por ejemplo, [Intel Digital Random Number Generator (DRNG) Software Implementation Guide](http://software.intel.com/en-us/articles/intel-digital-random-number-generator-drng-software-implementation-guide).


<!-- 
[HotBits](http://www.fourmilab.ch/hotbits): desintegración radiactiva.
-->

El uso de números aleatorios puros presenta dos grandes inconvenientes. 
El principal para su aplicación en el campo de la Estadística (y en otros casos) es que los valores generados deberían ser independientes e idénticamente distribuidos con distribución conocida, algo que resulta difícil (o imposible) de garantizar.
Siempre está presente la posible aparición de sesgos, principalmente debidos a fallos del sistema o interferencias. 
Por ejemplo, en el caso de la máquina RAND, fallos mecánicos en el sistema de grabación de los datos causaron problemas de aleatoriedad [@hacking1965, p. 118].

El otro inconveniente estaría relacionado con su reproducibilidad, por lo que habría que almacenarlos en tablas si se quieren volver a reproducir los resultados.
A partir de la década de 1960, al disponer de computadoras de mayor velocidad, empezó a resultar más eficiente generar valores mediante software en lugar de leerlos de tablas.


### Números cuasi-aleatorios 
 
<!-- 
Pendiente 
*números cuasi-aleatorios*: ... (se podría pensar que son una única generación de una variable aleatoria)?
-->

Algunos problemas, como la integración numérica (en el Capítulo \@ref(monte-carlo) se tratarán métodos de integración Monte Carlo), no dependen realmente de la aleatoriedad de la secuencia. Para evitar generaciones poco probables, se puede recurrir a secuencias cuasi-aleatorias, también denominadas *sucesiones de baja discrepancia* (hablaríamos entonces de métodos cuasi-Monte Carlo). La idea sería que la proporción de valores en una región cualquiera sea siempre aproximadamente proporcional a la medida de la región (como sucedería en media con la distribución uniforme, aunque no necesariamente para una realización concreta).

Por ejemplo, el paquete [`randtoolbox`](https://CRAN.R-project.org/package=randtoolbox) de R implementa métodos para la generación de secuencias cuasi-aleatorias (ver Figura \@ref(fig:randtoolbox)).


``` r
library(randtoolbox)
n <- 2000
par.old <- par(mfrow = c(1, 3))
plot(halton(n, dim = 2), xlab = 'x1', ylab = 'x2')
plot(sobol(n, dim = 2), xlab = 'x1', ylab = 'x2')
plot(torus(n, dim = 2), xlab = 'x1', ylab = 'x2')
par(par.old)
detach("package:randtoolbox")
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=1\linewidth]{01-Introduccion_files/figure-latex/randtoolbox-1} 

}

\caption{Secuencias cuasi-aleatorias bidimensionales obtenidas con los métodos de Halton (izquierda), Sobol (centro) y Torus (derecha).}(\#fig:randtoolbox)
\end{figure}

En este libro sólo consideraremos los números pseudoaleatorios y por comodidad se eliminará el prefijo "pseudo" en muchos casos.


### Números pseudo-aleatorios

La mayoría de los métodos de simulación se basan en la posibilidad de generar números pseudoaleatorios que imiten las propiedades de valores independientes de la distribución $\mathcal{U}(0,1)$, es decir, que imiten las propiedades de una muestra aleatoria simple^[Aunque hay que distinguir entre secuencia y muestra. En un problema de inferencia, en principio estamos interesados en una característica desconocida de la población. En cambio, en un problema de simulación "la población" es el modelo y lo conocemos por completo (no obstante el problema de simulación puede surgir como solución de un problema de inferencia).] de esta distribución. 

El procedimiento habitual para obtener estas secuencias es emplear un algoritmo recursivo denominado *generador*:

$$x_{i} = f\left( x_{i-1}, x_{i-2}, \cdots, x_{i-k}\right)$$

donde:

-   $k$ es el orden del generador.

-   $\left(  x_{0},x_{1},\cdots,x_{k-1}\right)$ es la *semilla*
  (estado inicial).

El *periodo* o *longitud del ciclo* es la longitud de la secuencia antes de que vuelva a repetirse. Lo denotaremos por $p$.


Los números de la sucesión son predecibles, conociendo el algoritmo y la semilla.
Sin embargo, si no se conociesen, *no se debería poder distinguir* una serie de números pseudoaleatorios *de una sucesión de números verdaderamente aleatoria* (utilizando recursos computacionales razonables).
En caso contrario esta predecibilidad puede dar lugar a serios
problemas (e.g. [http://eprint.iacr.org/2007/419](http://eprint.iacr.org/2007/419)).

Como regla general, por lo menos mientras se está desarrollando un
programa, interesa *fijar la semilla de aleatorización*.

-   Permite la reproducibilidad de los resultados.

-   Facilita la depuración del código.

Todo generador de números pseudoaleatorios mínimamente aceptable debe comportarse como si proporcionase muestras genuinas de datos independientes de una $\mathcal{U}(0,1)$.
Otras propiedades de interés son:

-   Reproducibilidad a partir de la semilla.

-   Periodo suficientemente largo.

-   Eficiencia (rapidez y requerimientos de memoria).

-   Portabilidad.

-   Generación de sub-secuencias (computación en paralelo).

-   Parsimonia.

Es importante asegurarse de que el generador empleado es adecuado:

> "Random numbers should not be generated with a method chosen at random." 
>
>
> --- Knuth, D.E. (TAOCP, 2002)


Se dispone de una gran cantidad de algoritmos. 
Los primeros intentos (cuadrados medios, método de Lehmer...) resultaron infructuosos, pero al poco tiempo ya se propusieron métodos que podían ser ampliamente utilizados (estableciendo adecuadamente sus parámetros).
Entre ellas podríamos destacar:

-   Generadores congruenciales.

-   Registros desfasados.

-   Combinaciones de distintos algoritmos.


La recomendación sería emplear un algoritmo conocido y que haya sido estudiado en profundidad. 
Por ejemplo el generador *Mersenne-Twister* propuesto por @matsumoto1998 y empleado por defecto en R (ver Sección \@ref(ext-cong)).
Además, sería recomendable utilizar alguna de las implementaciones disponibles en múltiples librerías, por ejemplo:

-   GNU Scientific Library (GSL):
    [http://www.gnu.org/software/gsl/manual](http://www.gnu.org/software/gsl/manual/html\_node/Random-Number-Generation.html)

-   StatLib: [http://lib.stat.cmu.edu](http://lib.stat.cmu.edu)

-   Numerical recipes: [http://www.nrbook.com/nr3](http://www.nrbook.com/nr3)

-   UNU.RAN (paquete `Runuran`):
    [http://statmath.wu.ac.at/unuran](http://statmath.wu.ac.at/unuran)

<!-- 
-   [http://random.mat.sbg.ac.at/software](http://random.mat.sbg.ac.at/software)

-   KISS (Keep It Simple Stupid / Small and Simple):
    [http://www.fortran.com/kiss.f90](http://www.fortran.com/kiss.f90)
-->
  
En este libro nos centraremos en los generadores congruenciales, descritos en la Sección \@ref(gen-cong).
Estos métodos son muy simples, aunque con las opciones adecuadas podrían ser utilizados en pequeños estudios de simulación. Sin embargo, su principal interés es que constituyen la base de los generadores avanzados habitualmente considerados.

## Números aleatorios en R {#rrng}

La generación de números pseudoaleatorios en R es una de las mejores
disponibles en paquetes estadísticos. 
Entre las herramientas implementadas en el paquete base de R podemos destacar:

-   [`set.seed()`](https://rdrr.io/r/base/Random.html): permite establecer la semilla (y el generador).

-   [`RNGkind()`](https://rdrr.io/r/base/Random.html): selecciona el generador.

-   `rdistribución(n,...):` genera valores aleatorios de la correspondiente distribución. 
    Por ejemplo, `runif(n, min = 0, max = 1)`, generaría `n` valores de una uniforme. Se puede acceder al listado completo de las funciones disponibles en el paquete `stats` mediante el comando [`?distributions`](https://rdrr.io/r/stats/Distributions.html).

-   [`sample()`](https://rdrr.io/r/base/sample.html): genera muestras aleatorias de variables discretas y permutaciones (se tratará en el Capítulo \@ref(discretas)).

-   [`simulate()`](https://rdrr.io/r/stats/simulate.html): genera realizaciones de la respuesta de un modelo ajustado.

Además están disponibles otros paquetes que implementan distribuciones adicionales (ver [CRAN Task View: Probability Distributions](https://cran.r-project.org/view=Distributions)). 
Entre ellos podríamos destacar los paquetes [`distr`](http://distr.r-forge.r-project.org) (clases S4; con extensiones en otros paquetes) y [`distributions3`](https://alexpghayes.github.io/distributions3/) (clases S3)^[El paquete [`distr6`](https://cran.r-project.org/package=distr6), que implementa clases R6, está actualmente archivado.].

La semilla se almacena en `.Random.seed`:

-   Inicialmente no existe. 
La recomendación es establecerla con `set.seed(entero)`, en caso contrario se generará a partir del reloj del sistema^[y del identificador asignado por el sistema operativo al proceso.] cuando se necesite.

-   Se almacena como un objeto oculto en el entorno de trabajo (o entorno global `.GlobalEnv`). Con las opciones por defecto de R, si al terminar una sesión almacenamos el entorno (en un fichero *.RData*), al iniciar una nueva sesión se restaurará también la semilla (y se podría continuar con las simulaciones).

-   Es un vector de enteros cuya estructura depende del tipo de generador (en la Sección \@ref(ext-cong) se dan algunos detalles sobre la configuración por defecto), por lo que no debería ser modificado manualmente.  
Puede ser recomendable almacenar (el objeto completo) antes de generar simulaciones, e.g. `seed <- .Random.seed`. Esto permite reproducir los resultados y facilita la depuración de posibles errores.

En la mayoría de los ejemplos de este libro se generan todos los valores de una vez,
se guardan y se procesan vectorialmente (normalmente empleando la función [`apply()`](https://rdrr.io/r/base/apply.html)).
En problemas mas complejos, en los que no es necesario almacenar todas las simulaciones,
puede ser preferible emplear un bucle para generar y procesar cada simulación iterativamente.
Por ejemplo podríamos emplear el siguiente esquema:


``` r
# Fijar semilla
set.seed(1)
for (isim in 1:nsim) {
  seed <- .Random.seed
  # Si se produce un error, podremos depurarlo ejecutando:
  #  .Random.seed <- seed
  # ································
  # Generar valores pseudoaleatorios
  # ································
}
```

o alternativamente fijar la semilla en cada iteración, por ejemplo:


``` r
for (isim in 1:nsim) {
  set.seed(isim)
  # ································
  # Generar valores pseudoaleatorios
  # ································
}
```


###  Opciones {#oprrng}

Normalmente no nos va a interesar cambiar las opciones por defecto de R para la generación de números pseudoaleatorios.
Para establecer estas opciones podemos emplear los argumentos `kind = NULL`, `normal.kind = NULL` y  `sample.kind = NULL` en las funciones [`RNGkind()`](https://rdrr.io/r/base/Random.html) o [`set.seed()`](https://rdrr.io/r/base/Random.html).
A continuación se muestran las distintas opciones (resaltando en negrita los valores por defecto):

-   `kind` especifica el generador pseudoaleatorio (uniforme):

    -   “Wichmann-Hill”: Ciclo $6.9536\times10^{12}$

    -   “Marsaglia-Multicarry”: Ciclo mayor de $2^{60}$

    -   “Super-Duper”: Ciclo aprox. $4.6\times10^{18}$ (S-PLUS)

    -   **“Mersenne-Twister”**: Ciclo $2^{19937}-1$ y equidistribution
        en 623 dimensiones.

    -   “Knuth-TAOCP-2002”: Ciclo aprox. $2^{129}$.

    -   “Knuth-TAOCP”

    -   “user-supplied”: permite emplear generadores adicionales.

-   `normal.kind` selecciona el método de generación de normales 
    (se tratará más adelante):
    “Kinderman-Ramage”, “Buggy Kinderman-Ramage”,
    “Ahrens-Dieter”, “Box-Muller”, **“Inversion”** , o “user-supplied”.
    
-   `sample.kind` selecciona el método de generación de uniformes discretas (el empleado por la función `sample()`, que cambió ligeramente^[Para evitar problemas de redondeo con tamaños extremadamente grandes; ver bug [PR#17494](https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17494).] a partir de la versión 3.6.0 de R): "Rounding" (versión anterior a 3.6.0) o **"Rejection"**.    
    
Estas opciones están codificadas (con índices comenzando en 0) en el primer componente de la semilla:

``` r
set.seed(1)
.Random.seed[1]
```

```
 ## [1] 10403
```
Los dos últimos dígitos se corresponden con el generador, las centenas con el método de generación de normales y las decenas de millar con el método uniforme discreto. 


###  Paquetes de R

Otros paquetes de R que pueden ser de interés:

-   [`setRNG`](https://CRAN.R-project.org/package=setRNG) contiene herramientas que facilitan operar con la semilla
    (e.g. dentro de funciones).

-   [`random`](https://CRAN.R-project.org/package=random) permite la descarga de números “true random” desde
    [RANDOM.ORG](https://www.random.org).

-   [`randtoolbox`](https://CRAN.R-project.org/package=randtoolbox) implementa generadores más recientes (`rngWELL`) y
    generación de secuencias cuasi-aleatorias.

-   [`RDieHarder`](https://CRAN.R-project.org/package=RDieHarder) implementa diversos contrastes para el análisis de la
    calidad de un generador y varios generadores.

-   [`Runuran`](https://CRAN.R-project.org/package=Runuran) interfaz para la librería UNU.RAN para la generación
    (automática) de variables aleatorias no uniformes [ver @hormann2004].

-   [`rsprng`](https://CRAN.R-project.org/package=rsprng), [`rstream`](https://CRAN.R-project.org/package=rstream) y [`rlecuyer`](https://CRAN.R-project.org/package=rlecuyer) implementan
    la generación de múltiples secuencias (para programación paralela).




<!-- 
PENDIENTE: Paquetes de simulación 
-   `gls`, `rngwell19937`, `randaes`, `SuppDists`, `lhs`, `mc2d`,
    `fOptions`, ...
-->

### Tiempo de CPU


La velocidad del generador suele ser una característica importante (también medir los tiempos, de cada iteración y de cada procedimento, en estudios de simulación). 
Para evaluar el rendimiento están disponibles en R distintas herramientas:

-   [`proc.time()`](https://rdrr.io/r/base/proc.time.html): permite obtener tiempo de computación real 
    y de CPU.
    
    ```
    tini <- proc.time()
    # Código a evaluar
    tiempo <- proc.time() - tini
    ```

-   `system.time(expresión)`: muestra el tiempo de computación (real y
    de CPU) de expresión.

Por ejemplo, podríamos emplear las siguientes funciones para ir midiendo los tiempos de CPU durante una simulación:


``` r
CPUtimeini <- function() {
  .tiempo.ini <<- proc.time()
  .tiempo.last <<- .tiempo.ini
}

CPUtimeprint <- function() {
  tmp <- proc.time()
  cat("Tiempo última operación:\n")
  print(tmp-.tiempo.last)
  cat("Tiempo total operación:\n")
  print(tmp-.tiempo.ini)
  .tiempo.last <<- tmp
}
```

Llamando a `CPUtimeini()` donde se quiere empezar a contar, y a `CPUtimeprint()` para imprimir el tiempo total y el tiempo desde la última llamada a una de estas funciones. 
Ejemplo:


``` r
funtest <- function(n) mad(runif(n)) 
CPUtimeini()
result1 <- funtest(10^6)
CPUtimeprint()
```

```
 ## Tiempo última operación:
 ##    user  system elapsed 
 ##    0.01    0.02    0.03 
 ## Tiempo total operación:
 ##    user  system elapsed 
 ##    0.01    0.02    0.03
```

``` r
result2 <- funtest(10^3)
CPUtimeprint()
```

```
 ## Tiempo última operación:
 ##    user  system elapsed 
 ##       0       0       0 
 ## Tiempo total operación:
 ##    user  system elapsed 
 ##    0.01    0.02    0.03
```

La función [`cpu.time()`](https://rubenfcasal.github.io/simres/reference/cpu.time.html) del paquete [`simres`](https://rubenfcasal.github.io/simres) implementa una aproximación similar:

-   `cpu.time(restart = TRUE)`: inicia el temporizador y almacena el tiempo de inicio.

-   `cpu.time()`: calcula el tiempo (real y de CPU) total (desde tiempo de inicio) y parcial (desde la última llamada a esta función).


Si se desea ir almacenando automáticamente los tiempos de computación, puede ser recomendable emplear las funciones [`tic()`](https://rdrr.io/pkg/tictoc/man/tic.html) y [`toc()`](https://rdrr.io/pkg/tictoc/man/tic.html) del paquete [`tictoc`](https://github.com/jabiru/tictoc):

- `tic("mensaje")`: inicia el temporizador y almacena el tiempo de inicio junto con el mensaje en una pila. 

- [`toc()`](https://rdrr.io/pkg/tictoc/man/tic.html): calcula el tiempo transcurrido desde la llamada correspondiente a [`tic()`](https://rdrr.io/pkg/tictoc/man/tic.html).



``` r
library(tictoc)
## Timing nested code
tic("outer")
   result1 <- funtest(10^6)
   tic("middle")
      result2 <- funtest(10^3)
      tic("inner")
         result3 <- funtest(10^2)
      toc() # inner
```

```
 ## inner: 0 sec elapsed
```

``` r
   toc() # middle
```

```
 ## middle: 0 sec elapsed
```

``` r
toc() # outer
```

```
 ## outer: 0.04 sec elapsed
```

``` r
## Timing in a loop and analyzing the results later using tic.log().
tic.clearlog()
for (i in 1:10)
{
   tic(i)
   result <- funtest(10^4)
   toc(log = TRUE, quiet = TRUE)
}
# log.txt <- tic.log(format = TRUE)
# log.lst <- tic.log(format = FALSE)
log.times <- do.call(rbind.data.frame, tic.log(format = FALSE))
str(log.times)
```

```
 ## 'data.frame':	10 obs. of  4 variables:
 ##  $ tic         : num  6.59 6.59 6.59 6.59 6.59 6.59 6.59 6.59 6.59 6.59
 ##  $ toc         : num  6.59 6.59 6.59 6.59 6.59 6.59 6.59 6.59 6.59 6.59
 ##  $ msg         : chr  "1" "2" "3" "4" ...
 ##  $ callback_msg: chr  "1: 0 sec elapsed" "2: 0 sec elapsed" "3: 0 sec "..
```

``` r
tic.clearlog()

# timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
log.times$timings <- with(log.times, toc - tic)
summary(log.times$timings)
```

```
 ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 ##       0       0       0       0       0       0
```

Hay que tener en cuenta que, por construcción, aunque se realicen en la mismas condiciones (en el mismo equipo), los tiempos de CPU en R pueden variar "ligeramente" entre ejecuciones.
Si se quieren estudiar tiempos de computación de forma más precisa, se recomendaría promediar los tiempos de varias ejecuciones.
Para ello se pueden emplear las herramientas del paquete [`microbenchmark`](https://github.com/joshuaulrich/microbenchmark/).
No obstante, para los fines de este libro no será necesaria tanta precisión. 

Finalmente, si los tiempos de computación no fuesen asumibles, para identificar los cuellos de botella y mejorar el código para optimizar la velocidad, podríamos emplear la función [`Rprof()`](https://rdrr.io/r/utils/Rprof.html). 
Esta función permite evaluar el rendimiento muestreando la pila en intervalos para determinar en que funciones se emplea el tiempo de computación.
Después de ejecutar `Rprof(fichero)` y el código a evaluar, llamando a `Rprof(NULL)` se desactiva el muestreo y con `summaryRprof(fichero)` se muestran los resultados (para analizarlos puede resultar de utilidad el paquete [`proftools`](https://CRAN.R-project.org/package=proftools)).


## Ejercicios

::: {.exercise #simpi}
<br>
  
Sea $(X,Y)$ es un vector aleatorio con distribución uniforme en el
cuadrado $[-1,1]\times\lbrack-1,1]$ de área 4.

a)  Aproximar mediante simulación $P\left(X + Y \leq 0 \right)$ y
    compararla con la probabilidad teórica (obtenida aplicando la
    regla de Laplace $\frac{\text{área favorable}}{\text{área posible}}$).


b)  Aproximar el valor de $\pi$ mediante simulación a partir de
    $P\left( X^2 +Y^2 \leq 1 \right)$.

:::

Ver solución en Sección \@ref(sol-simpi).



::: {.exercise #bernoulli name="Experimento de Bernoulli"}
<br>

Consideramos el experimento de Bernoulli consistente en el
lanzamiento de una moneda.

a)  Empleando la función `sample`, obtener 1000 simulaciones del
    lanzamiento de una moneda `(0 = cruz, 1 = cara)`, suponiendo que
    no está trucada. Aproximar la probabilidad de cara a partir de
    las simulaciones.
    
b)  En R pueden generarse valores de la distribución de Bernoulli
    mediante la función `rbinom(nsim, size=1, prob)`. Generar un
    gráfico de lineas considerando en el eje $X$ el número de
    lanzamientos (de 1 a 10000) y en el eje $Y$ la frecuencia
    relativa del suceso cara (puede ser recomendable emplear la
    función `cumsum`).
:::


Ver solución en Sección \@ref(sol-bernoulli).

<!-- 
Ejercicio propuesto
Von Neumann (1951) afirmó:

> "... in tossing a coin it is probably easier to make two consecutive tosses independent than to toss heads with probability exactly one-half. 
> If independence of successive tosses is assumed, we can reconstruct a 50-50 chance out of even a badly biased coin by tossing twice. 
> If we get heads-heads or tails-tails, we reject the tosses and try again. 
> If we get heads-tails (or tails-heads), we accept the result as heads (or tails)." 

Comprobar mediante simulación esta afirmación a partir de generaciones de lanzamientos
de una moneda con probabilidad de obtener cara igual a 0.4.
-->


::: {.exercise #circuito name="Simulación de un circuito"}
<br>
  
Simular el paso de corriente a través del circuito mostrado en la Figura \@ref(fig:circuito2), donde se muestran las probabilidades de que pase corriente por cada uno de los interruptores, que se suponen variables aleatorias de Bernoulli independientes.

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.5\linewidth]{images/circuito2} 

}

\caption{Esquema de un circuito eléctrico con interruptores aleatorios.}(\#fig:circuito2)
\end{figure}
:::



::: {.remark}
R maneja internamente los valores lógicos como 1 (`TRUE`) y 0 (`FALSE`).
Recíprocamente, cualquier número puede ser tratado como lógico (al estilo de C).
El entero 0 es equivalente a `FALSE` y cualquier entero distinto de 0 a `TRUE`.
:::


Ver solución en Sección \@ref(sol-circuito).


::: {.exercise #mere name="El problema del Caballero de Méré"}
<br>

En 1651, el Caballero de Méré le planteó a Pascal una pregunta
relacionada con las apuestas y los juegos de azar: ¿es ventajoso
apostar a que en cuatro lanzamientos de un dado se obtiene al menos
un seis? Este problema generó una fructífera correspondencia entre
Pascal y Fermat que se considera, simbólicamente, como el nacimiento
del Cálculo de Probabilidades.

a)  Escribir una función que simule el lanzamiento de $n$ dados. El
    parámetro de entrada es el número de lanzamientos $n$, que toma
    el valor 4 por defecto, y la salida debe ser `TRUE` si se
    obtiene al menos un 6 y `FALSE` en caso contrario.

b)  Utilizar la función anterior para simular $nsim=10000$ jugadas
    de este juego y calcular la proporción de veces que se gana la
    apuesta (obtener al menos un 6 en $n$ lanzamientos), usando
    $n=4$. Comparar el resultado con la probabilidad teórica
    $1-(5/6)^{n}$.
    
    
:::

Ver solución en Sección \@ref(sol-mere).


::: {.exercise #album name="variación del problema del coleccionista (cadena de Markov)"}
<br>

Continuando con el ejemplo de la Sección \@ref(ealbum)
(álbum con $n = 75$ cromos y sobres con $m = 6$). A partir de $nsim=2000$ simulaciones de coleccionistas de cromos, aproximar por simulación la evolución del proceso de compra de un coleccionista (número de cromos distintos dependiendo de los sobres comprados).

:::

Ver solución en Sección \@ref(sol-album).


