# Números aleatorios en R {#rrng}




La generación de números pseudoaleatorios en R es una de las mejores
disponibles en paquetes estadísticos. 
Entre las herramientas implementadas en el paquete base de R podemos destacar:

-   `set.seed(entero)`: permite establecer la semilla (y el generador).

-   `RNGkind()`: selecciona el generador.

-   `rdistribución(n,...):` genera valores aleatorios de la correspondiente distribución. 
    Por ejemplo, `runif(n, min = 0, max = 1)`, generaría `n` valores de una uniforme. Se puede acceder al listado completo de las funciones disponibles en el paquete `stats` mediante el comando `?distributions`.

-   `sample()`: genera muestras aleatorias de variables discretas y permutaciones (se tratará en el Capítulo \@ref(discretas)).

-   `simulate()`: genera realizaciones de la respuesta de un modelo ajustado.

Además están disponibles otros paquetes que implementan distribuciones adicionales (ver [CRAN Task View: Probability Distributions](https://cran.r-project.org/view=Distributions)). 
Entre ellos podríamos destacar los paquetes [`distr`](http://distr.r-forge.r-project.org) (clases S4; con extensiones en otros paquetes) y [`distr6`](https://alan-turing-institute.github.io/distr6/index.html) (clases R6).

La semilla se almacena (en `globalenv`) en `.Random.seed`; es un vector
de enteros cuya dimensión depende del tipo de generador:

-   No debe ser modificado manualmente; se guarda con el entorno de trabajo (por ejemplo, si se guarda al terminar la sesión en un fichero *.RData*, se restaurará la semilla al iniciar una nueva sesión y podremos continuar con las simulaciones).

-   Si no se especifica con `set.seed` (o no existe) se genera a partir del reloj del sistema.

-   Puede ser recomendable almacenarla antes de generar simulaciones, e.g. `seed <- .Random.seed`. 
    Esto permite reproducir los resultados y facilita la depuración de posibles errores.

En la mayoría de los ejemplos de este libro se generan todos los valores de una vez,
se guardan y se procesan vectorialmente (normalmente empleando la función `apply`).
En problemas mas complejos, en los que no es necesario almacenar todas las simulaciones,
puede ser preferible emplear un bucle para generar y procesar cada simulación iterativamente.
Por ejemplo podríamos emplear el siguiente esquema:


```r
# Fijar semilla
set.seed(1)
for (isim in 1:nsim) {
  seed <- .Random.seed
  # Si se produce un error, podremos depurarlo ejecutando:
  #  .Random.seed <- seed
  ...
  # Generar valores pseudoaleatorios
  ...
}
```

o alternativamente fijar la semilla en cada iteración, por ejemplo:


```r
for (isim in 1:nsim) {
  set.seed(isim)
  ...
  # Generar valores pseudoaleatorios
  ...
}
```


##  Opciones

Normalmente no nos va a interesar cambiar las opciones por defecto de R para la generación de números pseudoaleatorios.
Para establecer estas opciones podemos emplear los argumentos `kind = NULL`, `normal.kind = NULL` y  `sample.kind = NULL` en las funciones `RNGkind()` o `set.seed()`.
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

```r
set.seed(1)
.Random.seed[1]
```

```
## [1] 10403
```
Los dos últimos dígitos se corresponden con el generador, las centenas con el método de generación de normales y las decenas de millar con el método uniforme discreto. 


##  Paquetes de R

Otros paquetes de R que pueden ser de interés:

-   `setRNG` contiene herramientas que facilitan operar con la semilla
    (dentro de funciones,...).

-   `random` permite la descarga de números “true random” desde [RANDOM.ORG](https://www.random.org).

-   `randtoolbox` implementa generadores más recientes (`rngWELL`) y
    generación de secuencias cuasi-aleatorias.

-   `RDieHarder` implementa diversos contrastes para el análisis de la
    calidad de un generador y varios generadores.

-   [`Runuran`](http://statmath.wu.ac.at/unuran) interfaz para la librería UNU.RAN para la
    generación (automática) de variables aleatorias no uniformes (ver Hörmann et al., 2004).

-   `rsprng`, `rstream` y `rlecuyer` implementan la generación de múltiples
    secuencias (para programación paralela).

-   `gls`, `rngwell19937`, `randaes`, `SuppDists`, `lhs`, `mc2d`,
    `fOptions`, ...


<!-- 
PENDIENTE: Paquetes de simulación 
-->

## Ejercicios

\BeginKnitrBlock{exercise}
<span class="exercise" id="exr:simpi"><strong>(\#exr:simpi) </strong></span>
\EndKnitrBlock{exercise}

Sea $(X,Y)$ es un vector aleatorio con distribución uniforme en el
cuadrado $[-1,1]\times\lbrack-1,1]$ de área 4.

a)  Aproximar mediante simulación $P\left(X + Y \leq 0 \right)$ y
    compararla con la probabilidad teórica (obtenida aplicando la
    regla de Laplace $\frac{\text{área favorable}}{\text{área posible}}$).

    Generamos `nsim = 10000` valores del proceso bidimensional:

    
    ```r
    set.seed(1)
    nsim <- 10000
    x <- runif(nsim, -1, 1)
    y <- runif(nsim, -1, 1)
    ```

    La probabilidad teórica es 1/2 y la aproximación por simulación es la frecuencia relativa del suceso en los valores generados (para calcularla podemos aprovechar que R maneja internamente los valores lógicos como 1, `TRUE`, y 0, `FALSE`):
    
    
    ```r
    indice <- (x+y < 0)
    sum(indice)/nsim
    ```
    
    ```
    ## [1] 0.4996
    ```
    
    Alternativamente (la frecuencia relativa es un caso particular de la media) se puede obtener de forma más simple como:
    
    
    ```r
    mean(indice)
    ```
    
    ```
    ## [1] 0.4996
    ```


b)  Aproximar el valor de $\pi$ mediante simulación a partir de
    $P\left( X^2 +Y^2 \leq 1 \right)$.

    
    ```r
    set.seed(1)
    n <- 10000
    x <- runif(n, -1, 1)
    y <- runif(n, -1, 1)
    indice <- (x^2+y^2 < 1)
    mean(indice)
    ```
    
    ```
    ## [1] 0.7806
    ```
    
    ```r
    pi/4
    ```
    
    ```
    ## [1] 0.7853982
    ```
    
    ```r
    pi_aprox <- 4*mean(indice)
    pi_aprox
    ```
    
    ```
    ## [1] 3.1224
    ```

    Generamos el correspondiente gráfico (ver Figura \@ref(fig:simpiplot)) (los puntos con color negro tienen distribución uniforme en el círculo unidad; esto está relacionado con el método de aceptación-rechazo, ver Ejemplo \@ref(exm:ar-esfera), o con el denominado método *hit-or-miss*).
    
    
    ```r
    # Colores y símbolos dependiendo de si el índice correspondiente es verdadero:
    color <- ifelse(indice, "black", "red") 
    simbolo <- ifelse(indice, 1, 4)
    plot(x, y, pch = simbolo, col = color, 
         xlim = c(-1, 1), ylim = c(-1, 1), xlab="X", ylab="Y", asp = 1) 
         # asp = 1 para dibujar circulo
    symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
    symbols(0, 0, squares = 2, inches = FALSE, add = TRUE)
    ```
    
    \begin{figure}[!htb]
    
    {\centering \includegraphics[width=0.7\linewidth]{02-Numeros_aleatorios_R_files/figure-latex/simpiplot-1} 
    
    }
    
    \caption{Valores generados con distribución uniforme bidimensional, con colores y símbolos indicando si están dentro del círculo unidad.}(\#fig:simpiplot)
    \end{figure}

---    

\BeginKnitrBlock{exercise}\iffalse{-91-69-120-112-101-114-105-109-101-110-116-111-32-100-101-32-66-101-114-110-111-117-108-108-105-93-}\fi{}
<span class="exercise" id="exr:bernouilli"><strong>(\#exr:bernouilli)  \iffalse (Experimento de Bernoulli) \fi{} </strong></span>
\EndKnitrBlock{exercise}
Consideramos el experimento de Bernoulli consistente en el
lanzamiento de una moneda.

a)  Empleando la función `sample`, obtener 1000 simulaciones del
    lanzamiento de una moneda `(0 = cruz, 1 = cara)`, suponiendo que
    no está trucada. Aproximar la probabilidad de cara a partir de
    las simulaciones.

    
    ```r
    set.seed(1)
    nsim <- 10000
    x <- sample(c(cara = 1, cruz = 0), nsim, replace = TRUE, prob = c(0.5,0.5))
    mean(x)
    ```
    
    ```
    ## [1] 0.4953
    ```
    
    ```r
    barplot(100*table(x)/nsim, ylab = "Porcentaje") # Representar porcentajes 
    ```
    
    \begin{figure}[!htb]
    
    {\centering \includegraphics[width=0.7\linewidth]{02-Numeros_aleatorios_R_files/figure-latex/simberplot-1} 
    
    }
    
    \caption{Frecuencias relativas de los valores generados con distribución Bernoulli (aproximaciones por simulación de las probabilidades teóricas).}(\#fig:simberplot)
    \end{figure}

b)  En R pueden generarse valores de la distribución de Bernoulli
    mediante la función `rbinom(nsim, size=1, prob)`. Generar un
    gráfico de lineas considerando en el eje $X$ el número de
    lanzamientos (de 1 a 10000) y en el eje $Y$ la frecuencia
    relativa del suceso cara (puede ser recomendable emplear la
    función `cumsum`).

    
    ```r
    set.seed(1)
    nsim <- 1000
    p <- 0.4
    x <- rbinom(nsim, size = 1, prob = p) # Simulamos una Bernouilli
    # Alternativa programación: x <- runif(nsim) < p
    mean(x)
    ```
    
    ```
    ## [1] 0.394
    ```
    
    ```r
    n <- 1:nsim
    plot(n, cumsum(x)/n, type="l", ylab="Proporción de caras", 
         xlab="Número de lanzamientos", ylim=c(0,1))
    abline(h=p, lty=2, col="red")
    ```
    
    \begin{figure}[!htb]
    
    {\centering \includegraphics[width=0.7\linewidth]{02-Numeros_aleatorios_R_files/figure-latex/simberconv-1} 
    
    }
    
    \caption{Gráfico de convergencia de la aproximación por simulación a la probabilidad teórica.}(\#fig:simberconv)
    \end{figure}

---

\BeginKnitrBlock{exercise}\iffalse{-91-83-105-109-117-108-97-99-105-243-110-32-100-101-32-117-110-32-99-105-114-99-117-105-116-111-93-}\fi{}
<span class="exercise" id="exr:circuito"><strong>(\#exr:circuito)  \iffalse (Simulación de un circuito) \fi{} </strong></span>
\EndKnitrBlock{exercise}
Simular el paso de corriente a través del siguiente circuito, donde
figuran las probabilidades de que pase corriente por cada uno de los
interruptores:


\begin{center}\includegraphics[width=0.5\linewidth]{images/circuito2} \end{center}

Considerar que cada interruptor es una variable aleatoria de Bernoulli independiente
para simular 1000 valores de cada una de ellas.
    
\BeginKnitrBlock{remark}
\iffalse{} <span class="remark"><em>Nota: </em></span>  \fi{}R maneja internamente los valores lógicos como 1 (`TRUE`) y 0 (`FALSE`).
Recíprocamente, cualquier número puede ser tratado como lógico (al estilo de C).
El entero 0 es equivalente a `FALSE` y cualquier entero distinto de 0 a `TRUE`.
\EndKnitrBlock{remark}


```r
set.seed(1)
nsim <- 10000
x1 <- rbinom(nsim, size=1, prob=0.8)
x2 <- rbinom(nsim, size=1, prob=0.9)
z1 <- x1 | x2   # Operador lógico "O"
x3 <- rbinom(nsim, size=1, prob=0.6)
x4 <- rbinom(nsim, size=1, prob=0.5)
z2 <- x3 | x4
z3 <- z1 | z2
x5 <- rbinom(nsim, size=1, prob=0.7)
fin <- z3 & x5  # Operador lógico "Y"
mean(fin)
```

```
## [1] 0.692
```

---

\BeginKnitrBlock{exercise}\iffalse{-91-69-108-32-112-114-111-98-108-101-109-97-32-100-101-108-32-67-97-98-97-108-108-101-114-111-32-100-101-32-77-233-114-233-93-}\fi{}
<span class="exercise" id="exr:mere"><strong>(\#exr:mere)  \iffalse (El problema del Caballero de Méré) \fi{} </strong></span>
\EndKnitrBlock{exercise}
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

    
    ```r
    deMere <- function(n = 4){
      lanz <- sample(1:6, replace=TRUE, size=n)
      return(6 %in% lanz)
    }
    
    n <- 4
    lanz <- sample(1:6, replace=TRUE, size=n)
    lanz
    ```
    
    ```
    ## [1] 3 5 1 6
    ```
    
    ```r
    6 %in% lanz
    ```
    
    ```
    ## [1] TRUE
    ```


b)  Utilizar la función anterior para simular $nsim=10000$ jugadas
    de este juego y calcular la proporción de veces que se gana la
    apuesta (obtener al menos un 6 en $n$ lanzamientos), usando
    $n=4$. Comparar el resultado con la probabilidad teórica
    $1-(5/6)^{n}$.

    
    ```r
    set.seed(1)
    n <- 4
    nsim <- 10000
    mean(replicate(nsim, deMere(n)))
    ```
    
    ```
    ## [1] 0.5148
    ```
    
    ```r
    1-(5/6)^n
    ```
    
    ```
    ## [1] 0.5177469
    ```

---


\BeginKnitrBlock{exercise}\iffalse{-91-118-97-114-105-97-99-105-243-110-32-100-101-108-32-112-114-111-98-108-101-109-97-32-100-101-108-32-99-111-108-101-99-99-105-111-110-105-115-116-97-93-}\fi{}
<span class="exercise" id="exr:album"><strong>(\#exr:album)  \iffalse (variación del problema del coleccionista) \fi{} </strong></span>
\EndKnitrBlock{exercise}

Supongamos que tenemos un álbum con $n = 75$ cromos y para completarlo hay que comprar sobres con $m = 6$ cromos. A partir de $nsim=1000$ simulaciones de coleccionistas de cromos, obtener una aproximación por simulación a la respuesta de las siguientes cuestiones:

a)  ¿Cuál será en media el número de sobres que hay que comprar para completar la colección? 
    ¿Cuál sería el gasto medio si cada sobre cuesta 0.8€? 
    
b)  ¿Cuál sería el número mínimo de sobres para asegurar de que se completa la colección un 95% de las veces?


## Tiempo de CPU


La velocidad del generador suele ser una característica importante (también medir los tiempos, de cada iteración y de cada procedimento, en estudios de simulación). 
Para evaluar el rendimiento están disponibles en R distintas herramientas:

-   `proc.time()`: permite obtener tiempo de computación real y de
    CPU.
    
    ```
    tini <- proc.time()
    # Código a evaluar
    tiempo <- proc.time() - tini
    ```

-   `system.time(expresión)`: muestra el tiempo de computación (real y
    de CPU) de expresión.


Por ejemplo, podríamos emplear las siguientes funciones para
ir midiendo los tiempos de CPU durante una simulación:


```r
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

Llamando a `CPUtimeini()` donde se quiere empezar a contar, 
y a `CPUtimeprint()` para imprimir el tiempo total 
y el tiempo desde la última llamada a una de estas funciones. 
Ejemplo:


```r
funtest <- function(n) mad(runif(n)) 
CPUtimeini()
result1 <- funtest(10^6)
CPUtimeprint()
```

```
## Tiempo última operación:
##    user  system elapsed 
##    0.15    0.00    0.16 
## Tiempo total operación:
##    user  system elapsed 
##    0.15    0.00    0.16
```

```r
result2 <- funtest(10^3)
CPUtimeprint()
```

```
## Tiempo última operación:
##    user  system elapsed 
##    0.17    0.00    0.17 
## Tiempo total operación:
##    user  system elapsed 
##    0.32    0.00    0.33
```


Hay diversos paquetes que implementan herramientas similares, por ejemplo:

* El paquete `tictoc`:

  - `tic("mensaje")`: inicia el temporizador y almacena el tiempo de inicio junto con el mensaje en una pila. 
  
  - `toc()`: calcula el tiempo transcurrido desde la llamada correspondiente a `tic()`.


* La función [`cpu.time()`](https://rubenfcasal.github.io/npsp/reference/cpu.time.html) del paquete `npsp`:

  -   `cpu.time(restart = TRUE)`: inicia el temporizador y almacena el tiempo de inicio.
  
  -   `cpu.time()`: calcula el tiempo (real y de CPU) total (desde tiempo de inicio) y parcial (desde la última llamada a esta función).



```r
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

```r
   toc() # middle
```

```
## middle: 0.02 sec elapsed
```

```r
toc() # outer
```

```
## outer: 0.1 sec elapsed
```

```r
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
## 'data.frame':	10 obs. of  3 variables:
##  $ tic: num  5.14 5.14 5.14 5.14 5.15 5.15 5.15 5.15 5.15 5.15
##  $ toc: num  5.14 5.14 5.14 5.15 5.15 5.15 5.15 5.15 5.15 5.15
##  $ msg: chr  "1" "2" "3" "4" ...
```

```r
tic.clearlog()

# timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
log.times$timings <- with(log.times, toc - tic)
summary(log.times$timings)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.000   0.000   0.000   0.001   0.000   0.010
```


Hay que tener en cuenta que, por construcción, aunque se realicen en la mismas condiciones (en el mismo equipo), los tiempos de CPU en R pueden variar "ligeramente" entre ejecuciones.
Si se quieren estudiar tiempos de computación de forma más precisa, se recomendaría promediar los tiempos de varias ejecuciones.
Para ello se pueden emplear las herramientas del paquete [`microbenchmark`](https://CRAN.R-project.org/package=microbenchmark).
No obstante, para los fines de este libro no será necesaria tanta precisión. 

Finalmente, si los tiempos de computación no fuesen asumibles, para identificar los cuellos de botella y mejorar el código para optimizar la velocidad, podríamos emplear la función `Rprof(fichero)`. 
Esta función permite evaluar el rendimiento muestreando la pila en intervalos para determinar en que funciones se emplea el tiempo de computación.
Después de ejecutar el código, llamando a `Rprof(NULL)` se desactiva el muestreo y con `summaryRprof(fichero)` se muestran los resultados.






