# Soluciones ejercicios {#soluciones}

<!-- Apéndice \@ref(soluciones) -->



<!-- 
---
title: "Soluciones a algunos ejercicios"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    # pandoc_args: ["--number-offset", "0,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("25-Soluciones.Rmd")
knitr::purl("25-Soluciones.Rmd", documentation = 2)
knitr::spin("25-Soluciones.R",knit = FALSE)
-->

A continuación se muestran soluciones de algunos de los ejercicios no resueltos en el texto.

## Capítulo 1 [Introducción a la simulación](intro.html)

### Ejercicio [1.1](ejercicios.html#exr:simpi) {#sol-simpi}

Enunciado \@ref(exr:simpi):

Sea $(X,Y)$ es un vector aleatorio con distribución uniforme en el
cuadrado $[-1,1]\times\lbrack-1,1]$ de área 4.

a)  Aproximar mediante simulación $P\left(X + Y \leq 0 \right)$ y
    compararla con la probabilidad teórica (obtenida aplicando la
    regla de Laplace $\frac{\text{área favorable}}{\text{área posible}}$).

---

Generamos `nsim = 10000` valores del proceso bidimensional:


``` r
set.seed(1)
nsim <- 10000
x <- runif(nsim, -1, 1)
y <- runif(nsim, -1, 1)
```

La probabilidad teórica es 1/2 y la aproximación por simulación es la frecuencia relativa del suceso en los valores generados (para calcularla podemos aprovechar que R maneja internamente los valores lógicos como 1, `TRUE`, y 0, `FALSE`):


``` r
indice <- (x + y < 0)
sum(indice)/nsim
```

```
 ## [1] 0.4996
```

Alternativamente (la frecuencia relativa es un caso particular de la media) se puede obtener de forma más simple como:


``` r
mean(indice)
```

```
 ## [1] 0.4996
```

---

b)  Aproximar el valor de $\pi$ mediante simulación a partir de
    $P\left( X^2 +Y^2 \leq 1 \right)$.

---


``` r
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

``` r
pi/4
```

```
 ## [1] 0.7854
```

``` r
pi_aprox <- 4*mean(indice)
pi_aprox
```

```
 ## [1] 3.1224
```

Generamos el correspondiente gráfico (ver Figura \@ref(fig:simpiplot)) (los puntos con color negro tienen distribución uniforme en el círculo unidad; esto está relacionado con el método de aceptación-rechazo, ver Ejemplo \@ref(exm:ar-esfera), o con el denominado método *hit-or-miss*).


``` r
# Colores y símbolos dependiendo de si el índice correspondiente es verdadero:
color <- ifelse(indice, "black", "red") 
simbolo <- ifelse(indice, 1, 4)
plot(x, y, pch = simbolo, col = color, 
     xlim = c(-1, 1), ylim = c(-1, 1), xlab="X", ylab="Y", asp = 1) 
     # asp = 1 para dibujar circulo
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
symbols(0, 0, squares = 2, inches = FALSE, add = TRUE)
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/simpiplot-1} 

}

\caption{Valores generados con distribución uniforme bidimensional, con colores y símbolos indicando si están dentro del círculo unidad.}(\#fig:simpiplot)
\end{figure}




### Ejercicio [1.2](ejercicios.html#exr:bernoulli) {#sol-bernoulli}

Enunciado \@ref(exr:bernoulli):

Consideramos el experimento de Bernoulli consistente en el
lanzamiento de una moneda.

a)  Empleando la función `sample`, obtener 1000 simulaciones del
    lanzamiento de una moneda `(0 = cruz, 1 = cara)`, suponiendo que
    no está trucada. Aproximar la probabilidad de cara a partir de
    las simulaciones.
    
---    


``` r
set.seed(1)
nsim <- 10000
x <- sample(c(cara = 1, cruz = 0), nsim, replace = TRUE, prob = c(0.5,0.5))
mean(x)
```

```
 ## [1] 0.4953
```

``` r
barplot(100*table(x)/nsim, ylab = "Porcentaje") # Representar porcentajes 
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/simberplot-1} 

}

\caption{Frecuencias relativas de los valores generados con distribución Bernoulli (aproximaciones por simulación de las probabilidades teóricas).}(\#fig:simberplot)
\end{figure}

---

b)  En R pueden generarse valores de la distribución de Bernoulli
    mediante la función `rbinom(nsim, size=1, prob)`. Generar un
    gráfico de lineas considerando en el eje $X$ el número de
    lanzamientos (de 1 a 10000) y en el eje $Y$ la frecuencia
    relativa del suceso cara (puede ser recomendable emplear la
    función `cumsum`).
    
---


``` r
set.seed(1)
nsim <- 1000
p <- 0.4
x <- rbinom(nsim, size = 1, prob = p) # Simulamos una Bernoulli
# Alternativa programación: x <- runif(nsim) < p
mean(x)
```

```
 ## [1] 0.394
```

``` r
n <- 1:nsim
plot(n, cumsum(x)/n, type="l", ylab="Proporción de caras", 
     xlab="Número de lanzamientos", ylim=c(0,1))
abline(h=p, lty=2, col="red")
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/simberconv-1} 

}

\caption{Gráfico de convergencia de la aproximación por simulación a la probabilidad teórica.}(\#fig:simberconv)
\end{figure}


### Ejercicio [1.3](ejercicios.html#exr:circuito) {#sol-circuito}

Enunciado \@ref(exr:circuito):

Simular el paso de corriente a través del siguiente circuito, donde
figuran las probabilidades de que pase corriente por cada uno de los
interruptores:


\begin{center}\includegraphics[width=0.5\linewidth]{images/circuito2} \end{center}

Considerar que cada interruptor es una variable aleatoria de Bernoulli independiente
para simular 1000 valores de cada una de ellas.
    
::: {.remark}
R maneja internamente los valores lógicos como 1 (`TRUE`) y 0 (`FALSE`).
Recíprocamente, cualquier número puede ser tratado como lógico (al estilo de C).
El entero 0 es equivalente a `FALSE` y cualquier entero distinto de 0 a `TRUE`.
:::

---


``` r
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


### Ejercicio [1.4](ejercicios.html#exr:mere) {#sol-mere}

Enunciado \@ref(exr:mere) (el problema del Caballero de Méré):

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

---    


``` r
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

``` r
6 %in% lanz
```

```
 ## [1] TRUE
```

---    

b)  Utilizar la función anterior para simular $nsim=10000$ jugadas
    de este juego y calcular la proporción de veces que se gana la
    apuesta (obtener al menos un 6 en $n$ lanzamientos), usando
    $n=4$. Comparar el resultado con la probabilidad teórica
    $1-(5/6)^{n}$.

---    


``` r
set.seed(1)
n <- 4
nsim <- 10000
mean(replicate(nsim, deMere(n)))
```

```
 ## [1] 0.5148
```

``` r
1-(5/6)^n
```

```
 ## [1] 0.51775
```





### Ejercicio [1.5](ejercicios.html#exr:album) {#sol-album}

Enunciado \@ref(exr:album) (variación del problema del coleccionista, cadena de Markov):

Continuando con el ejemplo de la Sección \@ref(ealbum)
(álbum con $n = 75$ cromos y sobres con $m = 6$). A partir de $nsim=2000$ simulaciones de coleccionistas de cromos, aproximar por simulación la evolución del proceso de compra de un coleccionista (número de cromos distintos dependiendo de los sobres comprados).

---

Generamos `nsim = 2000` simulaciones de coleccionistas de cromos:


``` r
# Parámetros
n <- 75 # Número total de cromos
m <- 6  # Número de cromos en cada sobre
repe <- TRUE # Repetición de cromos en cada sobre

# Número de simulaciones
nsim <- 2000

# Resultados simulación
nsobres <- numeric(nsim)      # Número de sobres
evol <- vector("list", nsim)  # Evolución del número de cromos
# Por comodidad se podría haber fijado un número máximo de cromos
# evol <- matrix(nrow = max_len, ncol = nsim)

# Fijar semilla
set.seed(1)
# Bucle simulación
for (isim in 1:nsim) {
  # seed <- .Random.seed # .Random.seed <- seed
  # Simular
  album <- logical(n)
  evolucion <- c()
  i <- 0 # Número de sobres
  repeat{
    i <- i + 1
    sobre <- sample(n, m, replace = repe)
    album[sobre] <- TRUE
    ncromos <- sum(album)
    evolucion <- c(evolucion, ncromos)
    if (ncromos == n) {
      nsobres[isim] <- i
      evol[[isim]] <- evolucion
      break
    }
  }
}
# simres::plot.sr(nsobres)
```

`evol` contiene las realizaciones de la cadena de Markov.


``` r
# plot(evol[[1]], type = "l")
```

Combinar realizaciones del proceso (evoluciones del número de cromos):


``` r
# Se extienden a la máxima longitud
max_len <- max(lengths(evol)) # max(sapply(evol, length))
evol <- sapply(evol, function(x) c(x, rep(n, max_len - length(x))))
str(evol)
```

```
 ##  num [1:167, 1:2000] 6 12 16 21 23 25 30 34 37 38 ...
```

Aproximar cuantiles (intervalos de predicción):


``` r
alpha <- 0.05
limits <- apply(evol, 1, quantile, probs = c(alpha, 0.5, 1-alpha))
str(limits)
```

```
 ##  num [1:3, 1:167] 5 6 6 10 11 12 14 16 18 18 ...
 ##  - attr(*, "dimnames")=List of 2
 ##   ..$ : chr [1:3] "5%" "50%" "95%"
 ##   ..$ : NULL
```

Ejemplo, aproximación de los límites (y mediana) para el número de cromos en el álbum después de comprar 20 sobres:


``` r
limits[, 20]
```

```
 ##  5% 50% 95% 
 ##  55  60  64
```

``` r
hist(evol[20, ], breaks = "FD", freq = FALSE,
     main = "", xlab = "Número de cromos distintos (en 20 sobres)")
abline(v = limits[, 20], lty = 2)
```



\begin{center}\includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/unnamed-chunk-13-1} \end{center}

Representar las realizaciones del proceso y los intervalos de predicción puntuales:


``` r
matplot(1:max_len, evol, type = "l", col = "lightgray", lty = 1,
    xlab="Número de sobres", ylab="Número de cromos distintos")
matlines(1:max_len, t(limits), lty = c(2, 1, 2), col = 1)
```



\begin{center}\includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/unnamed-chunk-14-1} \end{center}


## Capítulo 2 [Generación de números pseudoaleatorios](gen-pseudo.html)


### Ejercicio [2.1](ejercicios.html#exr:RANDVN) {#sol-RANDVN}

Enunciado \@ref(exr:RANDVN):

Uno de los primeros generadores utilizados fue el denominado método de los
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

Este algoritmo está implementado en la función `simres::rvng()` (ver también `simres::rng()`; fichero [*rng.R*](R/rng.R)):


``` r
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
 ## <bytecode: 0x000001bcdf0c7388>
 ## <environment: namespace:simres>
```

Estudiar las características del generador de cuadrados medios a partir de una secuencia de 500 valores. 
Emplear únicamente métodos gráficos.

---


## Capítulo 5 [Simulación de variables discretas](discretas.html)

### Ejercicio [5.1](ejercicios-discretas.html#exr:mixta-cuantil) {#sol-mixta-cuantil}

Enunciado \@ref(exr:mixta-cuantil) (Simulación de una distribución mixta mediante el método de inversión generalizado):

Consideramos la variable aleatoria con función de distribución dada por: 
$$F(x)=\left\{
\begin{array}
[c]{cl}0 & \mbox{si $x<0$}\\
\frac{x}{2}+\frac{1}{10} & \mbox{si $x\in[0,\frac{1}{5})$}\\
x+\frac{1}{10} & \mbox{si $x\in[\frac{1}{5},\frac{9}{10}]$}\\
1 & \mbox{en otro caso}
\end{array}
\right.$$

Esta función está implementada en el siguiente código:

``` r
fdistr <- function(x) {
ifelse(x < 0, 0,
    ifelse(x < 1/5, x/2 + 1/10,
        ifelse(x <= 9/10, x + 1/10, 1) ) )
}
```

Como es una función vectorial podemos emplear `curve()` para representarla:

``` r
curve(fdistr, from = -0.1, to = 1.1, main = '')
abline(h = c(1/10, 2/10, 3/10), lty = 2) # Discontinuidades
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/mixta-cuantil-plot-sol-1} 

}

\caption{Función de distribución mixta (con discontinuidades en 0 y 1/5).}(\#fig:mixta-cuantil-plot-sol)
\end{figure}

**Nota**: Esta variable toma los valores 0 y 1/5 con probabilidad 1/10.

a)  Diseñar un algoritmo basado en el método de inversión generalizado 
    para generar observaciones de esta variable.
     
b)  Implementar el algoritmo en una función que permita generar $nsim$ 
    valores de esta variable.

---

a)  El algoritmo general es siempre el mismo. Empleando la función cuantil:
    $$Q\left( u\right) = \inf \left\{ x\in \mathbb{R}:F\left( x\right) 
    \geq u\right\},$$
    el algoritmo sería:
    
    1. Generar $U\sim \mathcal{U}\left( 0,1\right)$
    
    2. Devolver $X=Q\left( U\right)$
    
    En este caso concreto:
    
    1. Generar $U\sim \mathcal{U}\left( 0,1\right)$
    
    2. Si $U < \frac{1}{10}$ devolver $X = 0$
    
    3. En caso contrario, si $U < \frac{2}{10}$ devolver $X = 2(U - \frac{1}{10})$
    
    4. En caso contrario, si $U < \frac{3}{10}$ devolver $X = \frac{2}{10}$
    
    5. En caso contrario devolver $X = U - \frac{1}{10}$

       
b)  El algoritmo de simulación se puede implementar a partir de la función cuantil
    (vectorial):   

    
    ``` r
    # Función cuantil:
    fquant <- function(u) 
      ifelse(u < 1/10, 0,
             ifelse(u < 2/10, 2*(u - 1/10),
                    ifelse(u < 3/10, 1/5, u - 1/10) ) )
    # Función para generar nsim valores:
    rx <- function(nsim) fquant(runif(nsim))
    ```

    Ejemplo:
    
    
    ``` r
    set.seed(1)
    nsim <- 10^4
    system.time(simx <- rx(nsim))
    ```
    
    ```
     ##    user  system elapsed 
     ##       0       0       0
    ```
    
    ``` r
    hist(simx, breaks = "FD", freq = FALSE)
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/unnamed-chunk-18-1} \end{center}
    
    En este caso como no es una variable absolutamente continua mejor emplear 
    la función de distribución para compararla con la teórica:
    
    
    ``` r
    curve(ecdf(simx)(x), from= -0.1, to = 1.1, type = "s")
    curve(fdistr(x), type = "s", lty = 2, add = TRUE)
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{25-Soluciones_files/figure-latex/unnamed-chunk-19-1} \end{center}

