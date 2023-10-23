# Generación de números pseudoaleatorios {#gen-pseudo}




<!--
---
title: "Generación de números pseudoaleatorios"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    pandoc_args: ["--number-offset", "1,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("02-Generacion_numeros_aleatorios.Rmd")
knitr::purl("02-Generacion_numeros_aleatorios.Rmd", documentation = 2)
knitr::spin("02-Generacion_numeros_aleatorios.R",knit = FALSE)

PENDIENTE:
- Redactar ejemplo repetición de contrastes
-->

Como ya se comentó, los distintos métodos de simulación requieren disponer de secuencias de números pseudoaleatorios que imiten las propiedades de generaciones independientes de una distribución $\mathcal{U}(0,1)$. 
En primer lugar nos centraremos en el caso de los generadores congruenciales. A pesar de su simplicidad, podrían ser adecuados en muchos casos y constituyen la base de los generadores avanzados habitualmente considerados.
Posteriormente se dará una visión de las diferentes herramientas para estudiar la calidad de un generador de números pseudoaleatorios.

## Generadores congruenciales lineales {#gen-cong}

<!-- 
Pendiente: Incluir nota sobre generadores implementados en ordenadores que trabajan con números enteros o bits
-->

En los generadores congruenciales lineales se considera una combinación lineal de los últimos $k$ enteros generados y se calcula su resto al dividir por un entero fijo $m$. 
En el método congruencial simple (de orden $k = 1$), partiendo de una semilla inicial $x_0$, el algoritmo secuencial es el siguiente:
$$\begin{aligned}
x_{i}  & = (ax_{i-1}+c) \bmod m \\
u_{i}  & = \dfrac{x_{i}}{m} \\
i  & =1,2,\ldots
\end{aligned}$$ 
donde $a$ (*multiplicador*), $c$ (*incremento*) y $m$ (*módulo*) son enteros positivos^[Se supone además que $a$, $c$ y $x_0$ son menores que $m$, ya que, dadas las propiedades algebraicas de la suma y el producto en el conjunto de clases de resto módulo $m$ (que es un anillo), cualquier otra elección de valores mayores o iguales que $m$ tiene un equivalente verificando esta restricción.] fijados de antemano (los parámetros de este generador). Si $c=0$ el generador se denomina congruencial *multiplicativo* (Lehmer, 1951) y en caso contrario se dice que es *mixto* (Rotenburg, 1960).

Obviamente los parámetros y la semilla determinan los valores generados, que también se pueden obtener de forma no recursiva:
$$x_{i}=\left( a^{i}x_0+c\frac{a^{i}-1}{a-1}\right) \bmod m$$

Este método está implementado^[Aunque de forma no muy eficiente. Para evitar problemas computacionales, se recomienda realizar el cálculo de los valores empleando el método de Schrage (ver Bratley *et al.*, 1987; L'Ecuyer, 1988).] en la función `rlcg()` del paquete [`simres`](https://rubenfcasal.github.io/simres), imitando el funcionamiento del generador uniforme de R (ver también `simres::rng()`; fichero [*rng.R*](R/rng.R)):


```r
simres::rlcg
```

```
## function(n, seed = as.numeric(Sys.time()), a = 7^5, c = 0, m = 2^31 - 1) {
##   u <- numeric(n)
##   for(i in 1:n) {
##     seed <- (a * seed + c) %% m
##     u[i] <- seed/m # (seed + 1)/(m + 1)
##   }
##   # Almacenar semilla y parámetros
##   assign(".rng", list(seed = seed, type = "lcg",
##           parameters = list(a = a, c = c, m = m)), envir = globalenv())
##   # .rng <<- list(seed = seed, type = "lcg", parameters = list(a = a, c = c, m = m))
##   # Para continuar con semilla y parámetros:
##   #   with(.rng, rlcg(n, seed, parameters$a, parameters$c, parameters$m))
##   # Devolver valores
##   return(u)
## }
## <bytecode: 0x000001aee7a2a5a0>
## <environment: namespace:simres>
```


Ejemplos de parámetros:

-   $c=0$, $a=2^{16}+3=65539$ y $m=2^{31}$, generador *RANDU* de IBM
    (**no recomendable** como veremos más adelante).

-   $c=0$, $a=7^{5}=16807$ y $m=2^{31}-1$ (primo de Mersenne), Park y Miller (1988)
    *minimal standar*, empleado por las librerías IMSL y NAG.
    
-   $c=0$, $a=48271$ y $m=2^{31}-1$ actualización del *minimal standar* 
    propuesta por Park, Miller y Stockmeyer (1993).
    

A pesar de su simplicidad, una adecuada elección de los parámetros permite obtener de manera eficiente secuencias de números "aparentemente" i.i.d. $\mathcal{U}(0,1)$.
Durante los primeros años, el procedimiento habitual consistía en escoger $m$ de forma que se pudiera realizar eficientemente la operación del módulo, aprovechando la arquitectura del ordenador (por ejemplo $m = 2^{31}$ si se emplean enteros con signo de 32 bits). 
Posteriormente se seleccionaban $c$ y $a$ de forma que el período $p$ fuese lo más largo posible (o suficientemente largo), empleando los resultados mostrados a continuación.


::: {.theorem #hull-dobell name="Hull y Dobell, 1962"}
<br>
Un generador congruencial tiene período máximo ($p=m$) si y solo si:

1.  $c$ y $m$ son primos relativos (i.e. $m.c.d.(c, m) = 1$).

2.  $a-1$ es múltiplo de todos los factores primos de $m$ (i.e.
    $a \equiv 1 \bmod q$, para todo $q$ factor primo de $m$).

3.  Si $m$ es múltiplo de $4$, entonces $a-1$ también lo ha de
    ser (i.e. $m \equiv 0 \bmod 4\Rightarrow a \equiv
    1 \bmod 4$).
 
:::

Algunas consecuencias:

-   Si $m$ primo, $p=m$ si y solo si $a=1$.

-   Un generador multiplicativo no cumple la condición 1 ($m.c.d.(0, m)=m$).


::: {.theorem}
<br>
Un generador multiplicativo tiene período máximo ($p=m-1$) si:

1.  $m$ es primo.

2.  $a$ es una raiz primitiva de $m$ (i.e. el menor entero $q$ tal
    que $a^{q}=1 \bmod m$ es $q=m-1$).
   
:::

Sin embargo, además de preocuparse de la longitud del ciclo, sería mucho más importante que las secuencias generadas se comporten de forma similar a muestras aleatorias simple de una $\mathcal{U}(0,1)$. 
Uno de los principales problemas con este tipo de generadores (y con muchos otros) es que los valores generados pueden mostrar una estructura reticular.
Este es el caso del generador RANDU de IBM muy empleado en la década de los 70 y que muestra una clara estructura reticular cuando se consideran más de dos dimensiones (aunque ya se detectó este problema en 1963 se siguió utilizando hasta los 90).
Para ilustrar este problema podríamos emplear el conjunto de datos `randu` del paquete base `datasets` que contiene 400 tripletas de números sucesivos generados con la implementación de VAX/VMS 1.5 (de 1977, y que no se corrigió hasta la versión 2.0 de 1980). 
Aunque también podemos emplear el siguiente código^[La función `stats::embed(u, 3)` devuelve una matriz en la que la fila i-ésima contiene los valores `u[i+2]`, `u[i+1]` y `u[i]`. Además, en lugar de la función `plot3D::points3D()` se podría utilizar la función `plot3d()` del paquete `rgl`, y rotar la figura (pulsando con el ratón) para ver los hiperplanos `rgl::plot3d(xyz)`.] (ver Figura \@ref(fig:randu)):


```r
library(simres)
system.time(u <- rlcg(n = 9999, 
          seed = 543210, a = 2^16 + 3, c = 0, m = 2^31))
```

```
##    user  system elapsed 
##       0       0       0
```

```r
# xyz <- matrix(u, ncol = 3, byrow = TRUE)
xyz <- stats::embed(u, 3)
library(plot3D)
# points3D(xyz[,1], xyz[,2], xyz[,3], colvar = NULL, phi = 60, 
#          theta = -50, pch = 21, cex = 0.2)
points3D(xyz[,3], xyz[,2], xyz[,1], colvar = NULL, phi = 60, 
         theta = -50, pch = 21, cex = 0.2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/randu-1} 

}

\caption{Grafico de dispersión de tripletas del generador RANDU de IBM (contenidas en 15 planos).}(\#fig:randu)
\end{figure}

En general todos los generadores de este tipo van a presentar estructuras reticulares.
Marsaglia (1968) demostró que las k-uplas de un generadores multiplicativo están contenidas en a lo sumo $\left(k!m\right)^{1/k}$ hiperplanos paralelos (`r `trunc((factorial(3)*2^31)^(1/3))` como máximo con $k=3$ y $m=2^{31}$).
Por tanto habría que seleccionar adecuadamente los parámetros del generador congruencial de forma que la estructura reticular sea imperceptible, teniendo en cuenta el número de datos que se pretende generar (por ejemplo de forma que la distancia mínima entre los puntos sea próxima a la esperada en teoría). 
Para más detalles sobre la estructura reticular ver por ejemplo Ripley (1987, Sección 2.7).

<!-- 
PENDIENTE: 
Ejercicio aleatoriedad dígitos menos significativos (ejemplo sample)
9*a(n-2)-6*a(n-1)+a(n) = 0 mod 2^31 en RANDU con periodo 2^29
$a$ es una raiz primitiva de $m$ en Park y Miller
-->

Se han propuesto diversas pruebas (ver Sección \@ref(calgen)) para
determinar si un generador tiene problemas de este tipo y se han
realizado numerosos estudios para determinadas familias (e.g. Park y
Miller, 1988, estudiaron los multiplicadores adecuados para $m=2^{31}-1$).
En ciertos contextos muy exigentes (por ejemplo en criptografía), se recomienda    considerar un "periodo de seguridad" $\approx \sqrt{p}$ para evitar este tipo de problemas.

Aunque estos generadores tienen limitaciones en su capacidad para producir secuencias muy largas de números i.i.d. $\mathcal{U}(0,1)$, son un elemento básico en generadores más avanzados (incluyendo el empleado por defecto en R) como veremos en la siguiente sección.

    
::: {.example #congru512}
<br>

Consideramos el generador congruencial, de ciclo máximo, definido por: 
$$\begin{aligned}
x_{n+1}  & =(5x_{n}+1)\ \bmod\ 512,\nonumber\\
u_{n+1}  & =\frac{x_{n+1}}{512},\ n=0,1,\dots\nonumber
\end{aligned}$$


a)  Generar 500 valores de este generador, obtener el tiempo de CPU,
    representar su distribución mediante un histograma (en escala
    de densidades) y compararla con la densidad teórica.
   
    
    ```r
    set.rng(321, "lcg", a = 5, c = 1, m = 512)  # Establecer semilla y parámetros
    nsim <- 500
    system.time(u <- rng(nsim)) 
    ```
    
    ```
    ##    user  system elapsed 
    ##       0       0       0
    ```
    
    ```r
    hist(u, freq = FALSE)
    abline(h = 1)                   # Densidad uniforme
    ```
    
    \begin{figure}[!htb]
    
    {\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/ejcona-1} 
    
    }
    
    \caption{Histograma de los valores generados.}(\#fig:ejcona)
    \end{figure}

    En este caso concreto la distribución de los valores generados es aparentemente más uniforme de lo que cabría esperar, lo que induciría a sospechar de la calidad de este generador (ver Ejemplo \@ref(exm:congru512b) en Sección \@ref(calgen)).

b)  Calcular la media de las simulaciones (`mean`) y compararla con
    la teórica.
    
    La aproximación por simulación de la media teórica es:
    
    
    ```r
    mean(u)
    ```
    
    ```
    ## [1] 0.4999609
    ```
    
    La media teórica es 0.5. 
    Error absoluto $\ensuremath{3.90625\times 10^{-5}}$.

c)  Aproximar (mediante simulación) la probabilidad del intervalo
    $(0.4;0.8)$ y compararla con la teórica.

    La probabilidad teórica es 0.8 - 0.4 = 0.4
    
    La aproximación mediante simulación:
    
    
    ```r
    sum((0.4 < u) & (u < 0.8))/nsim
    ```
    
    ```
    ## [1] 0.402
    ```
    
    ```r
    mean((0.4 < u) & (u < 0.8))     # Alternativa
    ```
    
    ```
    ## [1] 0.402
    ```

:::


## Extensiones {#ext-cong}

Se han considerado diversas extensiones del generador congruencial lineal simple:

-   Lineal múltiple: 
    $x_{i}= a_0 + a_1 x_{i-1} + a_2 x_{i-2} + \cdots + a_{k} x_{i-k} \bmod m$,
    con periodo $p\leq m^{k}-1$.

-   No lineal: 
    $x_{i} = f\left(  x_{i-1}, x_{i-2}, \cdots, x_{i-k} \right) \bmod m$. 
    Por ejemplo $x_{i} = a_0 + a_1 x_{i-1} + a_2 x_{i-1}^2 \bmod m$.

-   Matricial: 
    $\boldsymbol{x}_{i} = A_0 + A_1\boldsymbol{x}_{i-1} 
    + A_2\boldsymbol{x}_{i-2} + \cdots 
    + A_{k}\boldsymbol{x}_{i-k} \bmod m$.

Un ejemplo de generador congruencia lineal múltiple es el denominado *generador de Fibonacci retardado* (Fibonacci-lagged generator; Knuth, 1969):
$$x_n = (x_{n-37} + x_{n-100}) \bmod 2^{30},$$
con un período aproximado de $2^{129}$ y que puede ser empleado en R (lo cual no sería en principio recomendable; ver [Knuth Recent News 2002](https://www-cs-faculty.stanford.edu/~knuth/news02.html#rng)) estableciendo `kind` a `"Knuth-TAOCP-2002"` o `"Knuth-TAOCP"` en la llamada a `set.seed()` o `RNGkind()`.
    
El generador *Mersenne-Twister* (Matsumoto y Nishimura, 1998), empleado por defecto en R, de periodo $2^{19937}-1$ y equidistribution en 623 dimensiones, se puede expresar como un generador congruencial matricial lineal.
En cada iteración (*twist*) genera 624 valores (los últimos componentes de la semilla son los 624 enteros de 32 bits correspondientes, el segundo componente es el índice/posición correspondiente al último valor devuelto; el conjunto de enteros solo cambia cada 624 generaciones).


```r
set.seed(1)
u <- runif(1)
seed <- .Random.seed
u <- runif(623)
sum(seed != .Random.seed) 
```

```
## [1] 1
```

```r
# Solo cambia el índice: 
seed[2]; .Random.seed[2]
```

```
## [1] 1
```

```
## [1] 624
```

```r
u <- runif(1)
# Cada 624 generaciones cambia el conjunto de enteros y el índice se inicializa
sum(seed != .Random.seed)
```

```
## [1] 624
```

```r
seed[2]; .Random.seed[2]
```

```
## [1] 1
```

```
## [1] 1
```


Un caso particular del generador lineal múltiple son los denominados *generadores de registros desfasados* (más relacionados con la criptografía).
Se generan bits de forma secuencial considerando $m=2$ y $a_{i} \in \left \{ 0,1\right \}$ y se van combinando $l$ bits para obtener valores en el intervalo $(0, 1)$, por ejemplo $u_i = 0 . x_{it+1} x_{it+2} \ldots x_{it+l}$, siendo $t$ un parámetro denominado *aniquilación* (Tausworthe, 1965). 
Los cálculos se pueden realizar rápidamente mediante operaciones lógicas (los sumandos de la combinación lineal se traducen en un "o" exclusivo XOR), empleando directamente los registros del procesador (ver por ejemplo, Ripley, 1987, Algoritmo 2.1).

Otras alternativas consisten en la combinanción de varios generadores, las más empleadas son:

-   Combinar las salidas: por ejemplo $u_{i}=\sum_{l=1}^L u_{i}^{(l)} \bmod 1$, donde $u_{i}^{(l)}$ es el $i$-ésimo valor obtenido con el generador $l$.

-   Barajar las salidas: por ejemplo se crea una tabla empleando un generador y se utiliza otro para seleccionar el índice del valor que se va a devolver y posteriormente actualizar.

<!-- 
PENDIENTE:  
Ejemplo combinar salidas  generador Wichmann-Hill (1982)  https://en.wikipedia.org/wiki/Wichmann%E2%80%93Hill 
-->


El generador *L'Ecuyer-CMRG* (L'Ecuyer, 1999), empleado como base para la generación de múltiples secuencias en el paquete `parallel`, combina dos generadores concruenciales lineales múltiples de orden $k=3$ (el periodo aproximado es $2^{191}$).


## Análisis de la calidad de un generador {#calgen}

Para verificar si un generador tiene las propiedades estadísticas deseadas hay disponibles una gran cantidad de test de hipótesis y métodos gráficos,
incluyendo métodos genéricos (de bondad de ajuste y aleatoriedad) y contrastes específicos para generadores aleatorios.
Se trata principalmente de contrastar si las muestras generadas son i.i.d. $\mathcal{U}\left(0,1\right)$ (análisis univariante).
Aunque los métodos más avanzados tratan de contrastar si las $d$-uplas:

$$(U_{t+1},U_{t+2},\ldots,U_{t+d}); \ t=(i-1)d, \ i=1,\ldots,m$$

son i.i.d. $\mathcal{U}\left(0,1\right)^{d}$ (uniformes independientes en el hipercubo; análisis multivariante).
En el Apéndice \@ref(gof-aleat) se describen algunos de estos métodos.

En esta sección emplearemos únicamente métodos genéricos, ya que también pueden ser de utilidad para evaluar generadores de variables no uniformes y para la construcción de modelos del sistema real (e.g. para modelar variables que se tratarán como entradas del modelo general). 
Sin embargo, los métodos clásicos pueden no ser muy adecuados para evaluar generadores de números pseudoaleatorios (ver L’Ecuyer y Simard, 2007).
La recomendación sería emplear baterías de contrastes recientes, como las descritas en la Subsección \@ref(baterias).

Hay que destacar algunas diferencias entre el uso de este tipo de métodos en inferencia y en simulación. 
Por ejemplo, si empleamos un constrate de hipótesis del modo habitual, desconfiamos del generador si la muestra (secuencia) no se ajusta a la distribución teórica (p-valor $\leq \alpha$).
En simulación, además, también se sospecha si se ajusta demasiado bien a la distribución teórica (p-valor $\geq1-\alpha$), lo que indicaría que no reproduce adecuadamente la variabilidad.

Uno de los contrastes más conocidos es el test chi-cuadrado de bondad de ajuste (`chisq.test` para el caso discreto). 
Aunque si la variable de interés es continua, habría que discretizarla (con la correspondiente perdida de información). 
Por ejemplo, se podría emplear la función `simres::chisq.cont.test()` (fichero [*test.R*](R/test.R)), que imita a las incluidas en R:


```r
simres::chisq.cont.test
```

```
## function(x, distribution = "norm", nclass = floor(length(x)/5),
##                             output = TRUE, nestpar = 0, ...) {
##   # Función distribución
##   q.distrib <- eval(parse(text = paste("q", distribution, sep = "")))
##   # Puntos de corte
##   q <- q.distrib((1:(nclass - 1))/nclass, ...)
##   tol <- sqrt(.Machine$double.eps)
##   xbreaks <- c(min(x) - tol, q, max(x) + tol)
##   # Gráficos y frecuencias
##   if (output) {
##     xhist <- hist(x, breaks = xbreaks, freq = FALSE,
##                   lty = 2, border = "grey50")
##     # Función densidad
##     d.distrib <- eval(parse(text = paste("d", distribution, sep = "")))
##     curve(d.distrib(x, ...), add = TRUE)
##   } else {
##     xhist <- hist(x, breaks = xbreaks, plot = FALSE)
##   }
##   # Cálculo estadístico y p-valor
##   O <- xhist$counts  # Equivalente a table(cut(x, xbreaks)) pero más eficiente
##   E <- length(x)/nclass
##   DNAME <- deparse(substitute(x))
##   METHOD <- "Pearson's Chi-squared test"
##   STATISTIC <- sum((O - E)^2/E)
##   names(STATISTIC) <- "X-squared"
##   PARAMETER <- nclass - nestpar - 1
##   names(PARAMETER) <- "df"
##   PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
##   # Preparar resultados
##   classes <- format(xbreaks)
##   classes <- paste("(", classes[-(nclass + 1)], ",", classes[-1], "]",
##                    sep = "")
##   RESULTS <- list(classes = classes, observed = O, expected = E,
##                   residuals = (O - E)/sqrt(E))
##   if (output) {
##     cat("\nPearson's Chi-squared test table\n")
##     print(as.data.frame(RESULTS))
##   }
##   if (any(E < 5))
##     warning("Chi-squared approximation may be incorrect")
##   structure(c(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL,
##                    method = METHOD, data.name = DNAME), RESULTS), class = "htest")
## }
## <bytecode: 0x000001aeecba2c50>
## <environment: namespace:simres>
```

::: {.example #congru512b name="análisis de un generador congruencial continuación"}
<br>
  
Continuando con el generador congruencial del Ejemplo \@ref(exm:congru512): 


```r
set.rng(321, "lcg", a = 5, c = 1, m = 512)  # Establecer semilla y parámetros
nsim <- 500
u <- rng(nsim)
```

Al aplicar el test chi-cuadrado obtendríamos:

<!-- PENDIENTE: evitar r-markdown en título figura -->


```r
chisq.cont.test(u, distribution = "unif", 
                nclass = 10, nestpar = 0, min = 0, max = 1)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/chisq-test-unif-1} 

}

\caption{Gráfico resultante de aplicar la función `chisq.cont.test()` comparando el histograma de los valores generados con la densidad uniforme.}(\#fig:chisq-test-unif)
\end{figure}

```
## 
## Pearson's Chi-squared test table
##                          classes observed expected  residuals
## 1  (-1.490116e-08, 1.000000e-01]       51       50  0.1414214
## 2  ( 1.000000e-01, 2.000000e-01]       49       50 -0.1414214
## 3  ( 2.000000e-01, 3.000000e-01]       49       50 -0.1414214
## 4  ( 3.000000e-01, 4.000000e-01]       50       50  0.0000000
## 5  ( 4.000000e-01, 5.000000e-01]       51       50  0.1414214
## 6  ( 5.000000e-01, 6.000000e-01]       51       50  0.1414214
## 7  ( 6.000000e-01, 7.000000e-01]       49       50 -0.1414214
## 8  ( 7.000000e-01, 8.000000e-01]       50       50  0.0000000
## 9  ( 8.000000e-01, 9.000000e-01]       50       50  0.0000000
## 10 ( 9.000000e-01, 9.980469e-01]       50       50  0.0000000
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  u
## X-squared = 0.12, df = 9, p-value = 1
```

Alternativamente, por ejemplo si solo se pretende aplicar el contraste, se podría emplear  la función `simres::freq.test()` (fichero [*test.R*](R/test.R))  para este caso particular (ver Sección \@ref(freq-test)).

Como se muestra en la Figura \@ref(fig:chisq-test-unif) el histograma de la secuencia generada es muy plano (comparado con lo que cabría esperar de una muestra de tamaño 500 de una uniforme), y consecuentemente el p-valor del contraste chi-cuadrado es prácticamente 1, lo que indicaría que este generador no reproduce adecuadamente la variabilidad de una distribución uniforme.   

Otro contraste de bondad de ajuste muy conocido es el test de Kolmogorov-Smirnov, implementado en `ks.test` (ver Sección \@ref(ks-test)). 
Este contraste de hipótesis compara la función de distribución bajo la hipótesis nula con la función de distribución empírica (ver Sección \@ref(empdistr)), representadas en la Figura \@ref(fig:empdistrunif):
    

```r
# Distribución empírica
curve(ecdf(u)(x), type = "s", lwd = 2)
curve(punif(x, 0, 1), add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/empdistrunif-1} 

}

\caption{Comparación de la distribución empírica de la secuencia generada con la función de distribución uniforme.}(\#fig:empdistrunif)
\end{figure}
Podemos realizar el contraste con el siguiente código:

```r
# Test de Kolmogorov-Smirnov
ks.test(u, "punif", 0, 1)
```

```
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  u
## D = 0.0033281, p-value = 1
## alternative hypothesis: two-sided
```

En la Sección \@ref(gof) se describen con más detalle estos contrastes de bondad de ajuste.

Adicionalmente podríamos estudiar la aleatoriedad de los valores generados (ver Sección \@ref(diag-aleat)), por ejemplo mediante un gráfico secuencial y el de dispersión retardado.


```r
plot(as.ts(u))
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/plot-sec-1} 

}

\caption{Gráfico secuencial de los valores generados.}(\#fig:plot-sec)
\end{figure}


```r
plot(u[-nsim],u[-1])
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/plot-ret-1} 

}

\caption{Gráfico de dispersión retardado de los valores generados.}(\#fig:plot-ret)
\end{figure}

Si se observa algún tipo de patrón indicaría dependencia (se podría considerar como una versión descriptiva del denominado “Parking lot test”), ver Ejemplo \@ref(exm:ret-gen).

También podemos analizar las autocorrelaciones (las correlaciones de $(u_{i},u_{i+k})$, con $k=1,\ldots,K$): 


```r
acf(u)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/plot-acf-1} 

}

\caption{Autocorrelaciones de los valores generados.}(\#fig:plot-acf)
\end{figure}
    
Por ejemplo, para contrastar si las diez primeras autocorrelaciones son nulas podríamos emplear el test de Ljung-Box:
    

```r
Box.test(u, lag = 10, type = "Ljung")
```

```
## 
## 	Box-Ljung test
## 
## data:  u
## X-squared = 22.533, df = 10, p-value = 0.01261
```

:::


### Repetición de contrastes

Los contrastes se plantean habitualmente desde el punto de vista de la inferencia estadística: se realiza una prueba sobre la única muestra disponible. 
Si se realiza una única prueba, en las condiciones de $H_0$ hay una probabilidad $\alpha$ de rechazarla.
En simulación tiene mucho más sentido realizar un gran número de pruebas:

-   La proporción de rechazos debería aproximarse al valor de
    $\alpha$ (se puede comprobar para distintos valores de $\alpha$).

-   La distribución del estadístico debería ajustarse a la teórica
    bajo $H_0$ (se podría realizar un nuevo contraste de bondad
    de ajuste).

-   Los p-valores obtenidos deberían ajustarse a una
    $\mathcal{U}\left(0,1\right)$ (se podría realizar también un
    contraste de bondad de ajuste).

Este procedimiento es también el habitual para validar un método de
contraste de hipótesis por simulación (ver Sección \@ref(contrastes)).

::: {.example #rep-test-randu}
<br>
  
Continuando con el generador congruencial RANDU, podemos pensar en estudiar la uniformidad de los valores generados empleando repetidamente el test chi-cuadrado:
  


```r
# Valores iniciales
set.rng(543210, "lcg", a = 2^16 + 3, c = 0, m = 2^31)  # Establecer semilla y parámetros
# set.seed(543210)
n <- 500
nsim <- 1000
estadistico <- numeric(nsim)
pvalor <- numeric(nsim)

# Realizar contrastes
for(isim in 1:nsim) {
  u <- rng(n)    # Generar
  # u <- runif(n)
  tmp <- freq.test(u, nclass = 100)
  # tmp <- chisq.cont.test(u, distribution = "unif", nclass = 100, 
  #     output = FALSE, nestpar = 0, min = 0, max = 1)
  estadistico[isim] <- tmp$statistic
  pvalor[isim] <- tmp$p.value
}
```

Por ejemplo, podemos comparar la proporción de rechazos observados con los que cabría esperar con los niveles de significación habituales:


```r
{
cat("Proporción de rechazos al 1% =", mean(pvalor < 0.01), "\n") # sum(pvalor < 0.01)/nsim
cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")   # sum(pvalor < 0.05)/nsim
cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")   # sum(pvalor < 0.1)/nsim
}
```

```
## Proporción de rechazos al 1% = 0.014 
## Proporción de rechazos al 5% = 0.051 
## Proporción de rechazos al 10% = 0.112
```

Las proporciones de rechazo obtenidas deberían comportarse como una aproximación por simulación de los niveles teóricos.
En este caso no se observa nada extraño, por lo que no habría motivos para sospechar de la uniformidad de los valores generados (aparentemente no hay problemas con la uniformidad de este generador).


Adicionalmente, si queremos estudiar la proporción de rechazos (el *tamaño del contraste*) para los posibles valores de $\alpha$, podemos emplear la distribución empírica del p-valor (proporción de veces que resultó menor que un determinado valor):


```r
# Distribución empírica
plot(ecdf(pvalor), do.points = FALSE, lwd = 2, 
     xlab = 'Nivel de significación', ylab = 'Proporción de rechazos')
abline(a = 0, b = 1, lty = 2)   # curve(punif(x, 0, 1), add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/rep-test-ecdf-1} 

}

\caption{Proporción de rechazos con los distintos niveles de significación.}(\#fig:rep-test-ecdf)
\end{figure}

<!-- 
curve(ecdf(pvalor)(x), type = "s", lwd = 2) 
-->


También podemos estudiar la distribución del estadístico del contraste.
En este caso, como la distribución bajo la hipótesis nula está implementada en R, podemos compararla fácilmente con la de los valores generados (debería ser una aproximación por simulación de la distribución teórica):


```r
# Histograma
hist(estadistico, breaks = "FD", freq = FALSE, main = "")
curve(dchisq(x, 99), add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/rep-test-est-1} 

}

\caption{Distribución del estadístico del constraste.}(\#fig:rep-test-est)
\end{figure}

Además de la comparación gráfica, podríamos emplear un test de bondad de ajuste para contrastar si la distribución del estadístico es la teórica bajo la hipótesis nula:


```r
# Test chi-cuadrado (chi-cuadrado sobre chi-cuadrado)
# chisq.cont.test(estadistico, distribution="chisq", nclass=20, nestpar=0, df=99)
# Test de Kolmogorov-Smirnov
ks.test(estadistico, "pchisq", df = 99)
```

```
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  estadistico
## D = 0.023499, p-value = 0.6388
## alternative hypothesis: two-sided
```

En este caso la distribución observada del estadístico es la que cabría esperar de una muestra de este tamaño de la distribución teórica, por tanto, según este criterio, aparentemente no habría problemas con la uniformidad de este generador (hay que recordar que estamos utilizando contrastes de hipótesis como herramienta para ver si hay algún problema con el generador, no tiene mucho sentido hablar de aceptar o rechazar una hipótesis).

En lugar de estudiar la distribución del estadístico de contraste  siempre podemos analizar la distribución del p-valor.
Mientras que la distribución teórica del estadístico depende del contraste y puede ser complicada, la del p-valor es siempre una uniforme.


```r
# Histograma
hist(pvalor, freq = FALSE, main = "")
abline(h=1) # curve(dunif(x,0,1), add=TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/rep-test-pval-1} 

}

\caption{Distribución del p-valor del constraste.}(\#fig:rep-test-pval)
\end{figure}

```r
# Test chi-cuadrado
# chisq.cont.test(pvalor, distribution="unif", nclass=20, nestpar=0, min=0, max=1)
# Test de Kolmogorov-Smirnov
ks.test(pvalor, "punif",  min = 0, max = 1)
```

```
## 
## 	Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  pvalor
## D = 0.023499, p-value = 0.6388
## alternative hypothesis: two-sided
```

Como podemos observar, obtendríamos los mismos resultados que al analizar la distribución del estadístico. 

Alternativamente podríamos emplear la función [`rephtest()`](https://rubenfcasal.github.io/simres/reference/rephtest.html) del paquete `simres` (fichero [*test.R*](R/test.R)):


```r
set.rng(543210, "lcg", a = 2^16 + 3, c = 0, m = 2^31)
# res <- rephtest(n = 30, test = chisq.cont.test, rand.gen = rng,
#          distribution = "unif", output = FALSE, nestpar = 0)
res <- rephtest(n = 30, test = freq.test, rand.gen = rng, nclass = 6)
str(res)
```

```
## List of 2
##  $ statistics: num [1:1000] 5.2 6.8 12.4 0.8 5.6 7.6 6.4 9.6 5.2 3.2 ...
##  $ p.values  : num [1:1000] 0.392 0.2359 0.0297 0.977 0.3471 ...
##  - attr(*, "class")= chr "rhtest"
##  - attr(*, "method")= chr "Chi-squared test for given probabilities"
##  - attr(*, "names.stat")= chr "X-squared"
##  - attr(*, "parameter")= Named num 5
##   ..- attr(*, "names")= chr "df"
```

```r
summary(res)
```

```
## Proportion of rejections:
##    1%    5%   10%   25%   50% 
## 0.013 0.054 0.096 0.255 0.544
```

```r
old.par <- par(mfrow = c(1, 2))
plot(res, 2:3)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.9\linewidth]{02-Generacion_numeros_aleatorios_files/figure-latex/rephtest-1} 

}

\caption{Distribución de los p-valores y proporción de rechazos.}(\#fig:rephtest)
\end{figure}

```r
par(old.par)
```

:::
  

### Baterías de contrastes {#baterias}

Hay numerosos ejemplos de generadores que pasaron diferentes test de uniformidad y aleatoriedad pero que fallaron estrepitosamente al considerar nuevos contrastes diseñados específicamente para generadores aleatorios (ver Marsaglia *et al.*, 1990). 
Por este motivo, el procedimiento habitual en la práctica es aplicar un número más o menos elevado de contrastes (de distinto tipo y difíciles de pasar, e.g. Marsaglia y Tsang, 2002), de forma que si el generador los pasa tendremos mayor confianza en que sus propiedades son las adecuadas.
Este conjunto de pruebas es lo que se denomina batería de contrastes. Una de las primeras se introdujo en Knuth (1969) y de las más recientes podríamos destacar:

-   Diehard tests (The Marsaglia Random Number CDROM, 1995):
    [http://www.stat.fsu.edu/pub/diehard (versión archivada el 2016-01-25)](https://web.archive.org/web/20160125103112/http://stat.fsu.edu/pub/diehard).

-   Dieharder (Brown y Bauer, 2003):
    [Dieharder Page](https://webhome.phy.duke.edu/~rgb/General/dieharder.php),
    paquete [`RDieHarder`](https://github.com/eddelbuettel/rdieharder).
    
-   TestU01 (L'Ecuyer y Simard, 2007): 
    [http://simul.iro.umontreal.ca/testu01/tu01.html](http://simul.iro.umontreal.ca/testu01/tu01.html).

-   NIST test suite (National Institute of Standards and Technology, USA, 2010): 
    [http://csrc.nist.gov/groups/ST/toolkit/rng](http://csrc.nist.gov/groups/ST/toolkit/rng).


Para más detalles, ver por ejemplo^[También puede ser de interés el enlace [Randomness Tests: A Literature Survey](http://www.ciphersbyritter.com/RES/RANDTEST.HTM) y la entidad certificadora (gratuita) en línea [CAcert](http://www.cacert.at/random).]:

-  Marsaglia, G. y Tsang, W.W. (2002). [Some difficult-to-pass tests of randomness](http://www.jstatsoft.org/v07/i03). Journal of Statistical Software, 7(3), 1-9.    
  
-  Demirhan, H. y Bitirim, N. (2016). [CryptRndTest: an R package for testing the cryptographic randomness](https://journal.r-project.org/archive/2016/RJ-2016-016/index.html). The R Journal, 8(1), 233-247.


Estas baterías de contrastes se suelen emplear si el generador va a ser utilizado en criptografía o si es muy importante la impredecibilidad (normalmente con generadores de números "verdaderamente aleatorios" por hardware).
Si el objetivo es únicamente obtener resultados estadísticos (como en nuestro caso) no sería tan importante que el generador no superase alguno de estos test.


Ejercicios
----------


::: {.exercise #RANDVN name="Método de los cuadrados medios"}
<br>
  
Uno de los primeros generadores utilizados fue el denominado método de los cuadrados medios propuesto por Von Neumann (1946). 
Con este procedimiento se generan números pseudoaleatorios de 4 dígitos de la siguiente forma:

i.  Se escoge un número de cuatro dígitos $x_0$ (semilla).

ii.   Se eleva al cuadrado ($x_0^2$) y se toman los cuatro dígitos centrales ($x_1$).

iii.   Se genera el número pseudo-aleatorio como $$u_1=\frac{x_1}{10^{4}}.$$

iv.  Volver al paso ii y repetir el proceso.

Para obtener los $k$ (número par) dígitos centrales de $x_{i}^2$
se puede utilizar que:
$$x_{i+1}=\left\lfloor \left(  x_{i}^2-\left\lfloor \dfrac{x_{i}^2}{10^{(2k-\frac{k}2)}}\right\rfloor 10^{(2k-\frac{k}2)}\right)
/10^{\frac{k}2}\right\rfloor$$ 
  
Estudiar las características de este generador a partir de una secuencia de 500 valores. 
Emplear únicamente métodos gráficos.  
  
:::

Este algoritmo está implementado en la función `simres::rvng()` (ver también `simres::rng()`; fichero [*rng.R*](R/rng.R)):


```r
simres::rvng
```

```
## function(n, seed = as.numeric(Sys.time()), k = 4) {
##   seed <- seed %% 10^k
##   aux <- 10^(2*k-k/2)
##   aux2 <- 10^(k/2)
##   u <- numeric(n)
##   for(i in 1:n) {
##     z <- seed^2
##     seed <- trunc((z - trunc(z/aux)*aux)/aux2)
##     u[i] <- seed/10^k
##   }
##   # Almacenar semilla y parámetros
##   assign(".rng", list(seed = seed, type = "vm", parameters = list(k = k)),
##       envir = globalenv())
##   # .rng <<- list(seed = seed, type = "vm", parameters = list(k = k))
##   # Para continuar con semilla y parámetros:
##   #   with(.rng, rvng(n, seed, parameters$k))
##   # Devolver valores
##   return(u)
## }
## <bytecode: 0x000001aeeda1f248>
## <environment: namespace:simres>
```


::: {.exercise #parkmiller}
<br>
  
Considerando el generador congruencial multiplicativo de parámetros $a=7^{5}=16807$, $c=0$ y $m=2^{31}-1$ (*minimal standar* de Park y Miller, 1988). 
¿Se observan los mismos problemas que con el algoritmo RANDU al considerar las tripletas $(x_{k},x_{k+1},x_{k+2})$?

:::
