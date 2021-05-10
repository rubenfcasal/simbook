# Generación de números pseudoaleatorios {#cap3}

<!-- 
PENDIENTE:

- Redactar ejemplo repetición de contrastes
  Etiquetar figuras

-->




Como ya se comentó, los distintos métodos de simulación requieren disponer de secuencias de números pseudoaleatorios que imiten las propiedades de generaciones independientes de una distribución $\mathcal{U}(0,1)$. 
En primer lugar nos centraremos en el caso de los generadores congruenciales. A pesar de su simplicidad, podrían ser adecuados en muchos casos y constituyen la base de los generadores avanzados habitualmente considerados.
Posteriormente se dará una visión de las diferentes herramientas para estudiar la calidad de un generador de números pseudoaleatorios.

## Generadores congruenciales {#gen-cong}

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

Este método está implementado en el fichero *RANDC.R*, tratando de imitar el funcionamiento de un generador de R (aunque de un forma no muy eficiente^[Para evitar problemas computacionales, se recomienda realizar el cálculo de los valores empleando el método de Schrage (ver Bratley *et al.*, 1987; L'Ecuyer, 1988).]):

```r
# --------------------------------------------------
# Generador congruencial de números pseudoaleatorios
# --------------------------------------------------

# initRANDC(semilla,a,c,m)
# -----------------------
#   Selecciona el generador congruencial
#   Por defecto RANDU de IBM con semilla del reloj
#   OJO: No se hace ninguna verificación de los parámetros
initRANDC <- function(semilla=as.numeric(Sys.time()), a=2^16+3, c=0, m=2^31) {
  .semilla <<- as.double(semilla) %% m  #Cálculos en doble precisión
  .a <<- a
  .c <<- c
  .m <<- m
  return(invisible(list(semilla=.semilla,a=.a,c=.c,m=.m))) #print(initRANDC())
}

# RANDC()
# -----------------------
#   Genera un valor pseudoaleatorio con el generador congruencial
#   Actualiza la semilla (si no existe llama a initRANDC)
RANDC <- function() {
    if (!exists(".semilla", envir=globalenv())) initRANDC()
    .semilla <<- (.a * .semilla + .c) %% .m
    return(.semilla/.m)
}

# RANDCN(n)
# -----------------------
#   Genera un vector de valores pseudoaleatorios con el generador congruencial
#   (por defecto de dimensión 1000)
#   Actualiza la semilla (si no existe llama a initRANDC)
RANDCN <- function(n=1000) {
    x <- numeric(n)
    for(i in 1:n) x[i]<-RANDC()
    return(x)
    # return(replicate(n,RANDC()))  # Alternativa más rápida    
}

initRANDC(543210)       # Fijar semilla 543210 para reproductibilidad
```


Ejemplos:

-   $c=0$, $a=2^{16}+3=65539$ y $m=2^{31}$, generador *RANDU* de IBM
    (**no recomendable**).

-   $c=0$, $a=7^{5}=16807$ y $m=2^{31}-1$ (primo de Mersenne), Park y Miller (1988)
    *minimal standar*, empleado por las librerías IMSL y NAG.
    
-   $c=0$, $a=48271$ y $m=2^{31}-1$ actualización del *minimal standar* 
    propuesta por Park, Miller y Stockmeyer (1993).
    
Los parámetros y la semilla determinan los valores generados:
$$x_{i}=\left(  a^{i}x_0+c\frac{a^{i}-1}{a-1}\right) \bmod m$$

A pesar de su simplicidad, una adecuada elección de los parámetros
permite obtener de manera eficiente secuencias de números
“aparentemente” i.i.d. $\mathcal{U}(0,1)$.

El procedimiento habitual solía ser escoger $m$ de forma que la operación del módulo se pudiese realizar de forma muy eficiente, para posteriormente seleccionar $c$ y $a$ de forma que el período fuese lo más largo posible (o suficientemente largo).


\BeginKnitrBlock{theorem}\iffalse{-91-72-117-108-108-32-121-32-68-111-98-101-108-108-44-32-49-57-54-50-93-}\fi{}<div class="theorem"><span class="theorem" id="thm:hull-dobell"><strong>(\#thm:hull-dobell)  \iffalse (Hull y Dobell, 1962) \fi{} </strong></span>
Un generador congruencial tiene período máximo ($p=m$) si y solo si:

1.  $c$ y $m$ son primos relativos (i.e. $m.c.d.(c, m) = 1$).

2.  $a-1$ es múltiplo de todos los factores primos de $m$ (i.e.
    $a \equiv 1 \bmod q$, para todo $q$ factor primo de $m$).

3.  Si $m$ es múltiplo de $4$, entonces $a-1$ también lo ha de
    ser (i.e. $m \equiv 0 \bmod 4\Rightarrow a \equiv
    1 \bmod 4$).
    
.   </div>\EndKnitrBlock{theorem}

Algunas consecuencias:

-   Si $m$ primo, $p=m$ si y solo si $a=1$.

-   Un generador multiplicativo no cumple la condición 1 ($m.c.d.(0, m)=m$).


\BeginKnitrBlock{theorem}<div class="theorem"><span class="theorem" id="thm:unnamed-chunk-3"><strong>(\#thm:unnamed-chunk-3) </strong></span>
Un generador multiplicativo tiene período máximo ($p=m-1$) si:

1.  $m$ es primo.

2.  $a$ es una raiz primitiva de $m$ (i.e. el menor entero $q$ tal
    que $a^{q}=1 \bmod m$ es $q=m-1$).

.    </div>\EndKnitrBlock{theorem}

Además de preocuparse de la longitud del ciclo, las secuencias
generadas deben aparentar muestras i.i.d. $\mathcal{U}(0,1)$. 

Uno de los principales problemas es que los valores generados pueden mostrar una clara estructura reticular.
Este es el caso por ejemplo del generador RANDU de IBM muy empleado en la década de los 70 (ver Figura \@ref(fig:randu))^[Alternativamente se podría utilizar la función `plot3d` del paquete `rgl`, y rotar la figura (pulsando con el ratón) para ver los hiperplanos:
`rgl::plot3d(xyz)`].
Por ejemplo, el conjunto de datos `randu` contiene 400 tripletas de números sucesivos obtenidos con la implementación de VAX/VMS 1.5 (1977).


```r
system.time(u <- RANDCN(9999))  # Generar
```

```
##    user  system elapsed 
##    0.03    0.00    0.03
```

```r
xyz <- matrix(u, ncol = 3, byrow = TRUE)
# xyz <- stats::embed(u, 3)

library(plot3D)
points3D(xyz[,1], xyz[,2], xyz[,3], colvar = NULL, phi = 60, 
         theta = -50, pch = 21, cex = 0.2)
```

<div class="figure" style="text-align: center">
<img src="03-Generacion_numeros_aleatorios_files/figure-html/randu-1.png" alt="Grafico de dispersión de tripletas del generador RANDU de IBM (contenidas en 15 planos)" width="70%" />
<p class="caption">(\#fig:randu)Grafico de dispersión de tripletas del generador RANDU de IBM (contenidas en 15 planos)</p>
</div>

En general todos los generadores de este tipo van a presentar estructuras reticulares.
Marsaglia (1968) demostró que las $k$-uplas de un generadores multiplicativo están contenidas en a lo sumo $\left(k!m\right)^{1/k}$ hiperplanos paralelos (para más detalles sobre la estructura reticular, ver por ejemplo Ripley, 1987, sección 2.7).
Por tanto habría que seleccionar adecuadamente $m$ y $c$ ($a$ solo influiría en la pendiente) de forma que la estructura reticular sea impreceptible teniendo en cuenta el número de datos que se pretende generar (por ejemplo de forma que la distancia mínima entre los puntos sea próxima a la esperada en teoría).

Se han propuesto diversas pruebas (ver Sección \@ref(calgen)) para
determinar si un generador tiene problemas de este tipo y se han
realizado numerosos estudios para determinadas familias (e.g. Park y
Miller, 1988, estudiaron que parámetros son adecuados para $m=2^{31}-1$).

-   En cualquier caso, se recomienda considerar un “periodo de
    seguridad” $\approx \sqrt{p}$ para evitar este tipo de problemas.

-   Aunque estos generadores tiene limitaciones en su capacidad para
    producir secuencias muy largas de números i.i.d. $\mathcal{U}(0,1)$,
    es un elemento básico en generadores más avanzados.


### Otros generadores

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

Un caso particular del generador lineal múltiple son los denominados *generadores de registros desfasados* (más relacionados con la Criptografía).
Se generan bits de forma secuencial considerando $m=2$ y $a_{i} \in \left \{ 0,1\right \}$ y se van combinando $l$ bits para obtener valores en el intervalo $(0, 1)$, por ejemplo $u_i = 0 . x_{it+1} x_{it+2} \ldots x_{it+l}$, siendo $t$ un parámetro denominado *aniquilación* (Tausworthe, 1965). 
Los cálculos se pueden realizar rápidamente mediante operaciones lógicas (los sumandos de la combinación lineal se traducen en un "o" exclusivo XOR), empleando directamente los registros del procesador (ver por ejemplo, Ripley, 1987, Algoritmo 2.1).

Otras alternativas consisten en la combinanción de varios generadores, las más empleadas son:

-   Combinar las salidas: por ejemplo $u_{i}=\sum_{l=1}^L u_{i}^{(l)} \bmod 1$, donde $u_{i}^{(l)}$ es el $i$-ésimo valor obtenido con el generador $l$.

-   Barajar las salidas: por ejemplo se crea una tabla empleando un generador y se utiliza otro para seleccionar el índice del valor que se va a devolver y posteriormente actualizar.

El generador *L'Ecuyer-CMRG* (L'Ecuyer, 1999), empleado como base para la generación de múltiples secuencias en el paquete `parallel`, combina dos generadores concruenciales lineales múltiples de orden $k=3$ (el periodo aproximado es $2^{191}$).
    
\BeginKnitrBlock{exercise}\iffalse{-91-65-110-225-108-105-115-105-115-32-100-101-32-117-110-32-103-101-110-101-114-97-100-111-114-32-99-111-110-103-114-117-101-110-99-105-97-108-93-}\fi{}<div class="exercise"><span class="exercise" id="exr:congru512"><strong>(\#exr:congru512)  \iffalse (Análisis de un generador congruencial) \fi{} </strong></span></div>\EndKnitrBlock{exercise}

Considera el generador congruencial definido por: 
$$\begin{aligned}
x_{n+1}  & =(5x_{n}+1)\ \bmod\ 512,\nonumber\\
u_{n+1}  & =\frac{x_{n+1}}{512},\ n=0,1,\dots\nonumber
\end{aligned}$$
(de ciclo máximo).

NOTA: El algoritmo está implementado en el fichero *RANDC.R* y se muestra en la Sección \@ref(gen-cong).

a)  Generar 500 valores de este generador, obtener el tiempo de CPU,
    representar su distribución mediante un histograma (en escala
    de densidades) y compararla con la densidad teórica.
   
    
    ```r
    initRANDC(321, 5, 1, 512)       # Establecer semilla y parámetros
    nsim <- 500
    system.time(u <- RANDCN(nsim))  # Generar
    ```
    
    ```
    ##    user  system elapsed 
    ##       0       0       0
    ```
    
    ```r
    hist(u, freq = FALSE)
    abline(h = 1)                   # Densidad uniforme
    ```
    
    <div class="figure" style="text-align: center">
    <img src="03-Generacion_numeros_aleatorios_files/figure-html/ejcona-1.png" alt="Histograma de los valores generados" width="70%" />
    <p class="caption">(\#fig:ejcona)Histograma de los valores generados</p>
    </div>

    En este caso concreto la distribución de los valores generados es aparentemente más uniforme de lo que cabría esperar, lo que induciría a sospechar de la calidad de este generador.

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
    Error absoluto $3.90625\times 10^{-5}$.

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


Análisis de la calidad de un generador {#calgen}
--------------------------------------

Para verificar si un generador tiene las propiedades estadísticas deseadas hay disponibles una gran cantidad de test de hipótesis y métodos gráficos,
incluyendo métodos genéricos (de bondad de ajuste y aleatoriedad) y contrastes específicos para generadores aleatorios.
Se trata principalmente de contrastar si las muestras generadas son i.i.d. $\mathcal{U}\left(0,1\right)$ (análisis univariante).
Aunque los métodos más avanzados tratan normalmente de contrastar si las $d$-uplas:

$$(U_{t+1},U_{t+2},\ldots,U_{t+d}); \ t=(i-1)d, \ i=1,\ldots,m$$

son i.i.d. $\mathcal{U}\left(0,1\right)^{d}$ (uniformes independientes en el hipercubo; análisis multivariante).
En el Apéndice \@ref(gof-aleat) se describen algunos de estos métodos.

En esta sección emplearemos únicamente métodos genéricos, ya que también pueden ser de utilidad para evaluar generadores de variables no uniformes y para la construcción de modelos del sistema real (e.g. para modelar variables que se tratarán como entradas del modelo general). 
Sin embargo, los métodos clásicos pueden no ser muy adecuados para evaluar generadores de números pseudoaleatorios (e.g. L’Ecuyer y Simard, 2007).
La recomendación sería emplear baterías de contrastes recientes, como las descritas en la Subsección \@ref(baterias).

Hay que destacar algunas diferencias entre el uso de este tipo de métodos en inferencia y en simulación. 
Por ejemplo, si empleamos un constrate de hipótesis del modo habitual, desconfiamos del generador si la muestra (secuencia) no se ajusta a la distribución teórica ($p$-valor $\leq \alpha$).
En este caso además, también se sospecha si se ajusta demasiado
bien a la distribución teórica ($p$-valor $\geq1-\alpha$),
lo que indicaría que no reproduce adecuadamente la variabilidad.

Uno de los contrastes más conocidos es el test ji-cuadrado de bondad de ajuste
(`chisq.test` para el caso discreto). 
Aunque si la variable de interés es continua, habría que discretizarla 
(con la correspondiente perdida de información). 
Por ejemplo, se podría emplear la siguiente función 
(que imita a las incluídas en `R`):



```r
#-------------------------------------------------------------------------------
# chisq.test.cont(x, distribution, nclasses, output, nestpar,...)
#-----------------------------------------------------------------------
# Realiza el test ji-cuadrado de bondad de ajuste para una distribución 
# continua discretizando en intervalos equiprobables.
# Parámetros:
#   distribution = "norm","unif", etc
#   nclasses = floor(length(x)/5)
#   output = TRUE
#   nestpar = 0 = número de parámetros estimados
#   ... = parámetros distribución
# Ejemplo:
#   chisq.test.cont(x, distribution = "norm", nestpar = 2, 
#                   mean = mean(x), sd = sqrt((nx - 1) / nx) * sd(x))
#-----------------------------------------------------------------------
chisq.test.cont <- function(x, distribution = "norm", 
    nclasses = floor(length(x)/5), output = TRUE, nestpar = 0, ...) {
    # Funciones distribución
    q.distrib <- eval(parse(text = paste("q", distribution, sep = "")))
    d.distrib <- eval(parse(text = paste("d", distribution, sep = "")))
    # Puntos de corte
    q <- q.distrib((1:(nclasses - 1))/nclasses, ...)
    tol <- sqrt(.Machine$double.eps)
    xbreaks <- c(min(x) - tol, q, max(x) + tol)
    # Gráficos y frecuencias
    if (output) {
        xhist <- hist(x, breaks = xbreaks, freq = FALSE, 
                      lty = 2, border = "grey50")
        curve(d.distrib(x, ...), add = TRUE)
    } else {
        xhist <- hist(x, breaks = xbreaks, plot = FALSE)
    }
    # Cálculo estadístico y p-valor
    O <- xhist$counts  # Equivalente a table(cut(x, xbreaks)) pero más eficiente
    E <- length(x)/nclasses
    DNAME <- deparse(substitute(x))
    METHOD <- "Pearson's Chi-squared test"
    STATISTIC <- sum((O - E)^2/E)
    names(STATISTIC) <- "X-squared"
    PARAMETER <- nclasses - nestpar - 1
    names(PARAMETER) <- "df"
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    # Preparar resultados
    classes <- format(xbreaks)
    classes <- paste("(", classes[-(nclasses + 1)], ",", classes[-1], "]", 
        sep = "")
    RESULTS <- list(classes = classes, observed = O, expected = E, residuals = (O - 
        E)/sqrt(E))
    if (output) {
        cat("\nPearson's Chi-squared test table\n")
        print(as.data.frame(RESULTS))
    }
    if (any(E < 5)) 
        warning("Chi-squared approximation may be incorrect")
    structure(c(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, 
        method = METHOD, data.name = DNAME), RESULTS), class = "htest")
}
```
Continuando con el generador congruencial anterior, obtendríamos:


```r
chisq.test.cont(u, distribution = "unif", 
                nclasses = 10, nestpar = 0, min = 0, max = 1)
```

<div class="figure" style="text-align: center">
<img src="03-Generacion_numeros_aleatorios_files/figure-html/chisq-test-unif-1.png" alt="Gráfico resultante de aplicar la función `chisq.test.cont()` comparando el histograma de los valores generados con la densidad uniforme." width="70%" />
<p class="caption">(\#fig:chisq-test-unif)Gráfico resultante de aplicar la función `chisq.test.cont()` comparando el histograma de los valores generados con la densidad uniforme.</p>
</div>

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

Como se muestra en la Figura \@ref(fig:chisq-test-unif) el histograma de la secuencia generada es muy plano (comparado con lo que cabría esperar de una muestra de tamaño 500 de una uniforme), y consecuentemente el $p$-valor del contraste ji-cuadrado es prácticamente 1, lo que indicaría que este generador no reproduce adecuadamente la variabilidad de una distribución uniforme.   

Otro contraste de bondad de ajuste muy conocido es el test de Kolmogorov-Smirnov, implementado en `ks.test` (ver Sección \@ref(ks-test)). 
En la Sección \@ref(gof) se describen con más detalle estos contrastes.

\BeginKnitrBlock{exercise}\iffalse{-91-65-110-225-108-105-115-105-115-32-100-101-32-117-110-32-103-101-110-101-114-97-100-111-114-32-99-111-110-103-114-117-101-110-99-105-97-108-44-32-99-111-110-116-105-110-117-97-99-105-243-110-93-}\fi{}<div class="exercise"><span class="exercise" id="exr:congru512b"><strong>(\#exr:congru512b)  \iffalse (Análisis de un generador congruencial, continuación) \fi{} </strong></span></div>\EndKnitrBlock{exercise}

Continuando con el generador congruencial del Ejercicio \@ref(exr:congru512): 


```r
initRANDC(321, 5, 1, 512)
nsim <- 500
system.time(u <- RANDCN(nsim))
```

a)  Realizar el contraste de Kolmogorov-Smirnov para estudiar el
    ajuste a una $\mathcal{U}(0,1)$.
    
    Este contraste de hipótesis compara la función de distribución bajo la hipótesis nula con la función de distribución empírica (ver Sección \@ref(empdistr)), representadas en la Figura \@ref(fig:empdistrunif):
    
    ```r
    # Distribución empírica
    curve(ecdf(u)(x), type = "s", lwd = 2)
    curve(punif(x, 0, 1), add = TRUE)
    ```
    
    <div class="figure" style="text-align: center">
    <img src="03-Generacion_numeros_aleatorios_files/figure-html/empdistrunif-1.png" alt="Comparación de la distribución empírica de la secuencia generada con la función de distribución uniforme." width="70%" />
    <p class="caption">(\#fig:empdistrunif)Comparación de la distribución empírica de la secuencia generada con la función de distribución uniforme.</p>
    </div>
    Podemos realizar el contraste con el siguiente código:
    
    ```r
    # Test de Kolmogorov-Smirnov
    ks.test(u, "punif", 0, 1)
    ```
    
    ```
    ## 
    ## 	One-sample Kolmogorov-Smirnov test
    ## 
    ## data:  u
    ## D = 0.0033281, p-value = 1
    ## alternative hypothesis: two-sided
    ```
    
b)  Obtener el gráfico secuencial y el de dispersión retardado, ¿se
    observa algún problema?

    Gráfico secuencial:
    
    
    ```r
    plot(as.ts(u))
    ```
    
    <div class="figure" style="text-align: center">
    <img src="03-Generacion_numeros_aleatorios_files/figure-html/plot-sec-1.png" alt="Gráfico secuencial de los valores generados." width="70%" />
    <p class="caption">(\#fig:plot-sec)Gráfico secuencial de los valores generados.</p>
    </div>
    
    Gráfico de dispersión retardado:
    
    
    ```r
    plot(u[-nsim],u[-1])
    ```
    
    <div class="figure" style="text-align: center">
    <img src="03-Generacion_numeros_aleatorios_files/figure-html/plot-ret-1.png" alt="Gráfico de dispersión retardado de los valores generados." width="70%" />
    <p class="caption">(\#fig:plot-ret)Gráfico de dispersión retardado de los valores generados.</p>
    </div>

c)  Estudiar las correlaciones del vector $(u_{i},u_{i+k})$, con
    $k=1,\ldots,10$. Contrastar si son nulas.

    Correlaciones:
    
    
    ```r
    acf(u)
    ```
    
    <div class="figure" style="text-align: center">
    <img src="03-Generacion_numeros_aleatorios_files/figure-html/plot-acf-1.png" alt="Autocorrelaciones de los valores generados." width="70%" />
    <p class="caption">(\#fig:plot-acf)Autocorrelaciones de los valores generados.</p>
    </div>
    
    Test de Ljung-Box:
    
    
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


### Repetición de contrastes

Los contrastes se plantean habitualmente desde el punto de vista de
la inferencia estadística en la práctica: se realiza una prueba
sobre la única muestra disponible. Si se realiza una única prueba, 
en las condiciones de $H_0$ hay
una probabilidad $\alpha$ de rechazarla.
En simulación tiene mucho más sentido realizar un gran número de
pruebas:

-   La proporción de rechazos debería aproximarse al valor de
    $\alpha$(se puede comprobar para distintos valores de $\alpha$).

-   La distribución del estadístico debería ajustarse a la teórica
    bajo $H_0$(se podría realizar un nuevo contraste de bondad
    de ajuste).

-   Los *p*-valores obtenidos deberían ajustarse a una
    $\mathcal{U}\left(0,1\right)$ (se podría realizar también un
    contraste de bondad de ajuste).

Este procedimiento es también el habitual para validar un método de
contraste de hipótesis por simulación.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:rep-test-randu"><strong>(\#exm:rep-test-randu) </strong></span></div>\EndKnitrBlock{example}

Consideramos el generador congruencial RANDU:


```r
# Valores iniciales
initRANDC(543210)   # Fijar semilla para reproductibilidad
# set.seed(543210)
n <- 500
nsim <- 1000
estadistico <- numeric(nsim)
pvalor <- numeric(nsim)

# Realizar contrastes
for(isim in 1:nsim) {
  u <- RANDCN(n)    # Generar
  # u <- runif(n)
  tmp <- chisq.test.cont(u, distribution="unif", 
            nclasses=100, output=FALSE, nestpar=0, min=0, max=1)
  estadistico[isim] <- tmp$statistic
  pvalor[isim] <- tmp$p.value
}
```

Proporción de rechazos:


```r
# cat("\nProporción de rechazos al 1% =", sum(pvalor < 0.01)/nsim, "\n")
cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
```

```
## 
## Proporción de rechazos al 1% = 0.014
```

```r
# cat("Proporción de rechazos al 5% =", sum(pvalor < 0.05)/nsim, "\n")
cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
```

```
## Proporción de rechazos al 5% = 0.051
```

```r
# cat("Proporción de rechazos al 10% =", sum(pvalor < 0.1)/nsim, "\n")
cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")
```

```
## Proporción de rechazos al 10% = 0.112
```

Análisis del estadístico contraste:


```r
# Histograma
hist(estadistico, breaks = "FD", freq=FALSE)
curve(dchisq(x,99), add=TRUE)
```

<img src="03-Generacion_numeros_aleatorios_files/figure-html/unnamed-chunk-12-1.png" width="70%" style="display: block; margin: auto;" />

```r
# Test ji-cuadrado
# chisq.test.cont(estadistico, distribution="chisq", nclasses=20, nestpar=0, df=99)
# Test de Kolmogorov-Smirnov
ks.test(estadistico, "pchisq", df=99)
```

```
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  estadistico
## D = 0.023499, p-value = 0.6388
## alternative hypothesis: two-sided
```

Análisis de los *p*-valores:


```r
# Histograma
hist(pvalor, freq=FALSE)
abline(h=1) # curve(dunif(x,0,1), add=TRUE)
```

<img src="03-Generacion_numeros_aleatorios_files/figure-html/unnamed-chunk-13-1.png" width="70%" style="display: block; margin: auto;" />

```r
# Test ji-cuadrado
# chisq.test.cont(pvalor, distribution="unif", nclasses=20, nestpar=0, min=0, max=1)
# Test de Kolmogorov-Smirnov
ks.test(pvalor, "punif",  min=0, max=1)
```

```
## 
## 	One-sample Kolmogorov-Smirnov test
## 
## data:  pvalor
## D = 0.023499, p-value = 0.6388
## alternative hypothesis: two-sided
```

Adicionalmente, si queremos estudiar la proporción de rechazos (el *tamaño del contraste*) para los posibles valores de $\alpha$, podemos emplear la distribución empírica del $p$-valor (proporción de veces que resultó menor que un determinado valor):


```r
# Distribución empírica
curve(ecdf(pvalor)(x), type = "s", lwd = 2, 
      xlab = 'Nivel de significación', ylab = 'Proporción de rechazos')
abline(a = 0, b = 1, lty = 2)   # curve(punif(x, 0, 1), add = TRUE)
```

<img src="03-Generacion_numeros_aleatorios_files/figure-html/unnamed-chunk-14-1.png" width="70%" style="display: block; margin: auto;" />


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
  
-  Demirhan, H. y Bitirim, N. (2016). [CryptRndTest: an R package for testing the cryptographic randomness](https://journal.r-project.org/archive/2016/RJ-2016-016/index.html). 
  The R Journal, 8(1), 233-247.



Ejercicios
----------


\BeginKnitrBlock{exercise}\iffalse{-91-77-233-116-111-100-111-32-100-101-32-108-111-115-32-99-117-97-100-114-97-100-111-115-32-109-101-100-105-111-115-93-}\fi{}<div class="exercise"><span class="exercise" id="exr:RANDVN"><strong>(\#exr:RANDVN)  \iffalse (Método de los cuadrados medios) \fi{} </strong></span></div>\EndKnitrBlock{exercise}
Uno de los primeros generadores fue el denominado método de los
cuadrados medios propuesto por Von Neumann (1946). Con este
procedimiento se generan números pseudoaleatorios de 4 dígitos de la
siguiente forma:

i.  Se escoge un número de cuatro dígitos $x_0$ (semilla).

ii.   Se eleva al cuadrado ($x_0^2$) y se toman los cuatro dígitos centrales ($x_1$).

iii.   Se genera el número pseudo-aleatorio como $$u_1=\frac{x_1}{10^{4}}.$$

iv.  Volver al paso ii y repetir el proceso.

Para obtener los $k$ (número par) dígitos centrales de $x_{i}^2$
se puede utilizar que:
$$x_{i+1}=\left\lfloor \left(  x_{i}^2-\left\lfloor \dfrac{x_{i}^2}{10^{(2k-\frac{k}2)}}\right\rfloor 10^{(2k-\frac{k}2)}\right)
/10^{\frac{k}2}\right\rfloor$$ 

El algoritmo está implementado en el fichero *RANDVN.R*:

```r
# -------------------------------------------------
# Generador Von Neumann de números pseudoaleatorios
# -------------------------------------------------

# initRANDVN(semilla,n)
# -----------------------
#   Inicia el generador 
#   n número de digitos centrales, 4 por defecto (debe ser un número par)
#   Por defecto semilla del reloj
#   OJO: No se hace ninguna verificación de los parámetros
initRANDVN <- function(semilla = as.numeric(Sys.time()), n = 4) {
  .semilla <<- as.double(semilla) %% 10^n  # Cálculos en doble precisión
  .n <<- n
  .aux <<- 10^(2*n-n/2)
  .aux2 <<- 10^(n/2)
  return(invisible(list(semilla=.semilla,n=.n)))
}

# RANDVN()
# -----------------------
#   Genera un valor pseudoaleatorio con el generador de Von Neumann.
#   Actualiza la semilla (si no existe llama a initRANDVN).
RANDVN <- function() {
    if (!exists(".semilla", envir=globalenv())) initRANDVN()
    z <- .semilla^2
    .semilla <<- trunc((z-trunc(z/.aux)*.aux)/.aux2)
    return(.semilla/10^.n)
}

# RANDVNN(n)
# -----------------------
#   Genera un vector de valores pseudoaleatorios, de dimensión `n` 
#   con elgenerador de Von Neumann.
#   Actualiza la semilla (si no existe llama a initRANDVN).
RANDVNN <- function(n = 1000) {
    x <- numeric(n)
    for(i in 1:n) x[i] <- RANDVN()
    return(x)
    # return(replicate(n,RANDVN()))  # Alternativa más rápida
}
```

Estudiar las características del
generador de cuadrados medios a partir de una secuencia de 500
valores. Emplear únicamente métodos gráficos.

\BeginKnitrBlock{exercise}<div class="exercise"><span class="exercise" id="exr:parkmiller"><strong>(\#exr:parkmiller) </strong></span></div>\EndKnitrBlock{exercise}
Considerando el generador congruencial multiplicativo de parámetros
$a=7^{5}=16807$, $c=0$ y $m=2^{31}-1$. ¿Se observan los mismos problemas 
que con el algoritmo RANDU al considerar las tripletas $(x_{k},x_{k+1},x_{k+2})$?
