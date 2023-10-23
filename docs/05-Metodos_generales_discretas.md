# Simulación de variables discretas {#discretas}

<!-- Capítulo \@ref(discretas) -->




<!-- 
---
title: "Simulación de variables discretas"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    pandoc_args: ["--number-offset", "4,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("05-Metodos_generales_discretas.Rmd")
knitr::purl("05-Metodos_generales_discretas.Rmd", documentation = 2)
knitr::spin("05-Metodos_generales_discretas.R",knit = FALSE)
-->


Se trata de simular una variable aleatoria discreta $X$ con función de masa de
probabilidad:
$$\begin{array}{l|ccccc}
 x_{i}  &  x_{1}  &  x_{2}  &  \cdots   &  x_{n}  &  \cdots   \\ \hline
 P\left( X=x_{i}\right)   &  p_{1}  &  p_{2}  &  \cdots   &  p_{n}  &  
\cdots  
\end{array}$$

Considerando como partida una $\mathcal{U}\left( 0,1\right)$, la
idea general consiste en asociar a cada posible valor $x_{i}$ de $X$
un subintervalo de $\left( 0,1\right)$ de longitud igual a la correspondiente
probabilidad.
Por ejemplo, como ya se mostró en capítulos anteriores, es habitual emplear
código de la forma:

```r
x <- runif(nsim) < p
```
para simular una distribución $Bernoulli(p)$.

Otro ejemplo es la distribución uniforme discreta 
en $\{1,2,\ldots,n\}$, con función de masa de probabilidad dada por
$$p_{j}=\frac{1}{n}\text{, para }j=1,2,\ldots n.$$ 
Pueden generarse valores de esta distribución de forma muy eficiente truncando la distribución uniforme:

::: {.conjecture #unif-discr name="distribución uniforme discreta"}
<br>

1. Generar $U\sim \mathcal{U}\left( 0,1\right)$.

2. Devolver $X=\left\lfloor nU\right\rfloor + 1$.

:::

Sin embargo, para generar variables discretas con dominio finito en R, si no se dispone de un algoritmo específico más eficiente, es recomendable emplear la función [`sample()`](https://rdrr.io/r/base/sample.html) del paquete base.
En el caso general[^discretas-1]:

```r
sample(x = valores, size = nsim, replace = TRUE, prob)
```
implementa eficientemente el método "alias" que describiremos más adelante, en la Sección \@ref(alias).

[^discretas-1]: En el caso de que `x = n` sea un entero (o un escalar) generará valores de una distribución uniforme discreta (equivalente a `sample(x = 1:n, ...)` o `sample.int(n, ...)`).


## Método de la transformación cuantil {#transcuant}

Este método es una adaptación del método de inversión (válido para el caso continuo) a distribuciones discretas, por lo que también se denomina *método de inversión generalizado*. 
En este caso, la función de distribución es:
$$F\left( x\right)  =\sum_{x_{j}\leq x}p_{j},$$
y la distribución de la variable aleatoria $F\left( X\right)$ no es uniforme (es una variable aleatoria discreta que toma los valores $F\left( x_{i} \right)$ con probabilidad $p_{i}$, $i=1,2,\ldots$).
Sin embargo, se puede generalizar el método de inversión a situaciones en las que $F$ no es invertible considerando la función cuantil.

Se define la función cuantil o inversa generalizada de una función de distribución $F$ como:
$$Q\left( u\right) =\inf \left\{ x\in \mathbb{R}:F\left( x\right) \geq
u\right\} ,\ \forall u\in \left( 0,1\right).$$
Si $F$ es invertible $Q=F^{-1}$.

::: {.theorem #invgen name="de inversión generalizada"}
<br>
  
Si $U\sim \mathcal{U}\left( 0,1\right)$, la variable aleatoria $Q\left( U\right)$ tiene función de distribución $F$.
:::

::: {.proof}

Bastaría ver que: 
$$Q\left( u\right) \leq x \Longleftrightarrow u\leq F(x).$$

Como $F$ es monótona y por la definición de $Q$: 
$$Q\left( u\right) \leq x \Rightarrow u \leq F(Q\left( u\right)) \leq F(x).$$
Por otro lado como $Q$ también es monótona: 
$$u \leq F(x) \Rightarrow Q\left( u\right) \leq Q(F(x)) \leq x$$
:::


A partir de este resultado se deduce el siguiente algoritmo general para simular una distribución de probabilidad discreta.

::: {.conjecture name="de transformación cuantil" }
<br>
  
1.  Generar $U\sim \mathcal{U}\left( 0,1\right)$.

2.  Devolver $X=Q\left( U\right)$.
:::

El principal problema es el cálculo de $Q\left( U\right)$. 
En este caso, suponiendo por comodidad que los valores que toma la variable están ordenados ($x_{1}<x_{2}<\cdots$), la función cuantil será:
$$\begin{array}{ll}
Q\left( U\right) &=\inf \left\{ x_{j}:\sum_{i=1}^{j}p_{i}\geq U\right\} \\
&=x_{k}\text{, tal que }\sum_{i=1}^{k-1}p_{i}<U\leq \sum_{i=1}^{k}p_{i}
\end{array}$$

Para encontrar este valor se puede emplear el siguiente algoritmo:

::: {.conjecture name="de transformación cuantil con búsqueda secuencial" }
<br>
  
1.  Generar $U\sim \mathcal{U}\left( 0,1\right)$.

2.  Hacer $I=1$ y $S=p_{1}$.

3.  Mientras $U>S$ hacer $I=I+1$ y $S=S+p_{I}$

4.  Devolver $X=x_{I}$.
:::

Este algoritmo no es muy eficiente, especialmente si el número de posibles valores de la variable es grande.


\begin{figure}[!htb]

{\centering \includegraphics[width=0.8\linewidth]{images/cuantil-sim} 

}

\caption{Ilustración de la simulación de una distribución discreta mediante transformación cuantil (con búsqueda secuencial).}(\#fig:cuantil-movie)
\end{figure}

::: {.remark}
El algoritmo anterior es válido independientemente de que los valores que toma la variable estén ordenados.
:::

Si la variable toma un número finito de valores, se podría implementar en R de la siguiente forma:


```r
rpmf0 <- function(x, prob = 1/length(x), n = 1000) {
  X <- numeric(n)
  U <- runif(n)
  for(j in 1:n) {
    i <- 1
    Fx <- prob[1]
    while (Fx < U[j]) {
      i <- i + 1
      Fx <- Fx + prob[i]
    }
    X[j] <- x[i]
  }
  return(X)
}
```

Adicionalmente, para disminuir ligeramente el tiempo de computación, se pueden almacenar las probabilidades acumuladas en una tabla. 
Este algoritmo está implementado en la función  [`rpmf()`](https://rubenfcasal.github.io/simres/reference/rpmf.html) del paquete [`simres`](https://rubenfcasal.github.io/simres) (fichero [*rpmf.R*](R/rpmf.R)), que ademas devuelve el número de comparaciones en un atributo `ncomp`:


```r
library(simres)
rpmf
```

```
 ## function(x, prob = 1/length(x), n = 1000, as.factor = FALSE) {
 ##   # Numero de comparaciones
 ##   ncomp <- 0
 ##   # Inicializar FD
 ##   Fx <- cumsum(prob)
 ##   # Simular
 ##   X <- numeric(n)
 ##   U <- runif(n)
 ##   for(j in 1:n) {
 ##     i <- 1
 ##     while (Fx[i] < U[j]) i <- i + 1
 ##     X[j] <- x[i]
 ##     ncomp <- ncomp + i
 ##   }
 ##   if(as.factor) X <- factor(X, levels = x)
 ##   attr(X, "ncomp") <- ncomp
 ##   return(X)
 ## }
 ## <bytecode: 0x000001b7dca82ec0>
 ## <environment: namespace:simres>
```

::: {.example #binom-cuant name="Simulación de una binomial mediante el método de la transformación cuantil"}
<br>

Empleamos la rutina anterior para para generar una muestra de $nsim=10^{5}$ observaciones de una variable $\mathcal{B}(10,0.5)$ y obtenemos el tiempo de CPU empleado: 


Calcular también la
media muestral (compararla con la teórica $np$) y el número
medio de comparaciones para generar cada observación.

Empleamos la rutina anterior para generar las simulaciones:

```r
set.seed(1)
n <- 10
p <- 0.5
nsim <- 10^5
x <- 0:n
pmf <- dbinom(x, n, p)
system.time( rx <- rpmf(x, pmf, nsim) )
```

```
 ##    user  system elapsed 
 ##    0.06    0.00    0.06
```

A partir de ellas podríamos aproximar el valor esperado:

```r
mean(rx)
```

```
 ## [1] 4.99697
```
aunque en este caso el valor teórico es conocido $np$ = 5.

Calculamos el número medio de comparaciones para generar cada observación:

```r
ncomp <- attr(rx, "ncomp")
ncomp/nsim
```

```
 ## [1] 5.99697
```

```r
# Se verá más adelante que el valor teórico es sum((1:length(x))*pmf)
```

Representamos la aproximación por simulación de la función de masa de probabilidad y la comparamos con la teórica:

```r
res <- table(rx)/nsim
# res <- table(factor(rx, levels = x))/nsim
plot(res, ylab = "frecuencia relativa", xlab = "valor")
points(x, pmf, pch = 4, col = "blue")  # Comparación teórica
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{05-Metodos_generales_discretas_files/figure-latex/comprfmp-1} 

}

\caption{Comparación de las frecuencias relativas de los valores generados, mediante el método de la transformación cuantil, con las probabilidades teóricas.}(\#fig:comprfmp)
\end{figure}

También podríamos realizar comparaciones numéricas:

```r
res <- as.data.frame(res)
names(res) <- c("x", "psim")
res$pteor <- pmf
print(res, digits = 2)
```

```
 ##     x    psim   pteor
 ## 1   0 0.00100 0.00098
 ## 2   1 0.00981 0.00977
 ## 3   2 0.04457 0.04395
 ## 4   3 0.11865 0.11719
 ## 5   4 0.20377 0.20508
 ## 6   5 0.24477 0.24609
 ## 7   6 0.20593 0.20508
 ## 8   7 0.11631 0.11719
 ## 9   8 0.04415 0.04395
 ## 10  9 0.01021 0.00977
 ## 11 10 0.00083 0.00098
```

```r
# Máximo error absoluto 
max(abs(res$psim - res$pteor))
```

```
 ## [1] 0.0014625
```

```r
# Máximo error absoluto porcentual 
100*max(abs(res$psim - res$pteor) / res$pteor)
```

```
 ## [1] 15.008
```
:::

<!-- 
PENDIENTE: Error estándar de las aproximaciones de las probabilidades
-->

::: {.remark}
Puede ocurrir que no todos los valores sean generados en la simulación.
En el código anterior si `length(x) > length(psim)`, la sentencia `res$pteor <- pmf` generará un error.
Una posible solución sería trabajar con factores (llamar a la función `rpmf()` con `as.factor = TRUE` o emplear `res <- table(factor(rx, levels = x))/nsim`). 
:::

<!-- 
# Alternativamente se podría emplear `match()`: 
res <- data.frame(x = x, pteor = pmf, psim = 0)
res.sim <- table(rx)/nsim
index <- match(names(res.sim), x)
res$psim[index] <- res.sim 
-->


### Eficiencia del algoritmo

Si consideramos la variable aleatoria $\mathcal{I}$ correspondiente a las etiquetas, su función de masa de probabilidad sería:
$$\begin{array}{l|ccccc}
i & 1 & 2 & \cdots & n & \cdots \\ \hline
P\left( \mathcal{I}=i\right) & p_{1} & p_{2} & \cdots & p_{n} & \cdots 
\end{array}$$
y el número de comparaciones en el paso 3 sería un valor aleatorio de esta variable.
Una medida de la eficiencia del algoritmo de la transformación cuantil es el número medio de comparaciones[^transcuant-1]:
$$E\left( \mathcal{I}\right) =\sum_i ip_{i}.$$

[^transcuant-1]: Realmente, cuando la variable toma un número finito de valores: $x_{1}$, $x_{2}$, $\ldots$, $x_{n}$, no sería necesario hacer
la última comprobación $U>\sum_{i=1}^{n}p_{i}=1$ y se generaría directamente $x_{n}$, por lo que el número medio de comparaciones sería: $$\sum_{i=1}^{n-1}ip_{i}+\left( n-1\right)  p_{n}.$$

Para disminuir el número esperado de comparaciones podemos
reordenar los valores $x_{i}$ de forma que las probabilidades
correspondientes sean decrecientes. Esto equivale a considerar
un etiquetado $l$ de forma que:
$$p_{l\left( 1\right) }\geq p_{l\left( 2\right) }\geq \cdots \geq p_{l\left(
n\right) }\geq \cdots$$


::: {.example #binom-cuantb name="Simulación de una binomial, continuación"}
<br>

Podemos repetir la simulación del Ejemplo \@ref(exm:binom-cuant) anterior ordenando previamente las probabilidades en orden decreciente y también empleando la función `sample()` de R.


```r
set.seed(1)
tini <- proc.time()
# Ordenar
ind <- order(pmf, decreasing = TRUE)
# Generar
rx <- rpmf(x[ind], pmf[ind], nsim)
# Tiempo de CPU
tiempo <- proc.time() - tini
tiempo
```

```
 ##    user  system elapsed 
 ##    0.03    0.00    0.03
```

```r
# Número de comparaciones
ncomp <- attr(rx, "ncomp")
ncomp/nsim
```

```
 ## [1] 3.08369
```

```r
sum((1:length(x))*pmf[ind]) # Valor teórico
```

```
 ## [1] 3.083984
```

Como ya se comentó, en R se recomienda emplear la función `sample()` 
(implementa eficientemente el método de Alias descrito en la Sección \@ref(alias)):


```r
system.time( rx <- sample(x, nsim, replace = TRUE, prob = pmf) )
```

```
 ##    user  system elapsed 
 ##    0.01    0.00    0.01
```

:::

## Método de la tabla guía

También conocido como método de búsqueda indexada (*Indexed Search*), la idea consiste en construir $m$ subintervalos equiespaciados en $[0,1]$ de la forma:
$$I_{j}=\left[ u_{j},u_{j+1}\right) =\left[ \frac{j-1}{m},\frac{j}{m}\right) 
\text{ para }j=1,2,\ldots ,m$$
y utilizarlos como punto de partida para la búsqueda.
En una tabla guía se almacenan los indices de los cuantiles correspondientes a los extremos inferiores de los intervalos:
$$g_{j}=Q_{\mathcal{I}}(u_{j})=\inf \left\{ i:F_{i}\geq u_{j}=\frac{j-1}{m}\right\}$$

El punto de partida para un valor $U$ será $g_{j_{0}}$ con:
$$j_{0}=\left\lfloor mU\right\rfloor +1$$


\begin{figure}[!htb]

{\centering \includegraphics[width=0.8\linewidth]{images/tabla-sim} 

}

\caption{Ilustración de la simulación de una distribución discreta mediante tabla guía.}(\#fig:tabla-movie)
\end{figure}

En este caso, puede verse que una cota del número medio de comparaciones es:
$$E\left( N\right) \leq 1+\frac{n}{m}$$

::: {.conjecture name="de simulación mediante tabla guía; Chen y Asau, 1974" }
<br>

Inicialización:

1.  Hacer $F_{1}=p_{1}$.

2.  Desde $i=2$ hasta $n$ hacer $F_{i}=F_{i-1}+p_{i}$.

Cálculo de la tabla guía:

1.  Hacer $g_{1}=1$ e $i=1$.

2.  Desde $j=2$ hasta $m$ hacer

    2.a  Mientras $(j-1)/m>F_{i}$ hacer $i=i+1$.

    2.b  $g_{j}=i$

Simulación mediante tabla guía:

1.  Generar $U\sim \mathcal{U}\left( 0,1\right)$.

2.  Hacer $j=\left\lfloor mU\right\rfloor +1$.

3.  Hacer $i=g_{j}$.

4.  Mientras $U>F_{i}$ hacer $i=i+1$.

5.  Devolver $X=x_{i}$.

:::

Este algoritmo está implementado en la función  [`simres::rpmf.table()`](https://rubenfcasal.github.io/simres/reference/rpmf.table.html) (fichero [*rpmf.R*](R/rpmf.R)) y devuelve también el número de comparaciones en un atributo `ncomp`:


```r
rpmf.table
```

```
 ## function(x, prob = 1/length(x), m, n = 1000, as.factor = FALSE) {
 ##   # Inicializar tabla y FD
 ##   Fx <- cumsum(prob)
 ##   g <- rep(1,m)
 ##   i <- 1
 ##   for(j in 2:m) {
 ##     while (Fx[i] < (j-1)/m) i <- i + 1
 ##     g[j] <- i
 ##   }
 ##   ncomp <- i - 1
 ##   # Generar valores
 ##   X <- numeric(n)
 ##   U <- runif(n)
 ##   for(j in 1:n) {
 ##     i <- i0 <- g[floor(U[j] * m) + 1]
 ##     while (Fx[i] < U[j]) i <- i + 1
 ##     ncomp <- ncomp + i - i0
 ##     X[j] <- x[i]
 ##   }
 ##   if(as.factor) X <- factor(X, levels = x)
 ##   attr(X, "ncomp") <- ncomp
 ##   return(X)
 ## }
 ## <bytecode: 0x000001b7e1982108>
 ## <environment: namespace:simres>
```

::: {.example #binom-tabla name="Simulación de una binomial mediante tabla guía"}
<br>

Repetimos la simulación del Ejemplo \@ref(exm:binom-cuant) anterior empleando esta rutina con $m=n-1$.


```r
set.seed(1)
system.time( rx <- rpmf.table(x, pmf, n-1, nsim) )
```

```
 ##    user  system elapsed 
 ##    0.06    0.00    0.06
```

Número medio de comparaciones:

```r
ncomp <- attr(rx, "ncomp")
ncomp/nsim
```

```
 ## [1] 0.55951
```

```r
sum((1:length(x))*pmf) # Numero esperado con búsqueda secuencial
```

```
 ## [1] 6
```

Análisis de los resultados:

```r
res <- table(rx)/nsim
plot(res, ylab = "frecuencia relativa", xlab = "valores")
points(x, pmf, pch = 4, col = "blue")  # Comparación teórica
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{05-Metodos_generales_discretas_files/figure-latex/comptabla-1} 

}

\caption{Comparación de las frecuencias relativas de los valores generados, mediante el método de la tabla guía, con las probabilidades teóricas.}(\#fig:comptabla)
\end{figure}

:::

## Método de Alias {#alias}

Se basa en representar la distribución de $X$ como una mixtura
(uniforme) de variables dicotómicas (Walker, 1977):
$$Q^{(i)}=\left\{ 
\begin{array}{ll}
x_{i} & \text{con prob. } q_{i} \\ 
x_{a_{i}} & \text{con prob. } 1-q_{i}
\end{array}
\ \right.$$

Hay varias formas de construir las tablas de probabilidades $q_i$ y de alias $a_i$.
Se suele emplear el denominado algoritmo “Robin Hood” de inicialización (Kronmal y Peterson, 1979).
La idea es "tomar prestada" parte de la probabilidad de los valores más probables (ricos) para asignársela a los valores menos probables (pobres), recordando el valor de donde procede (almacenando el índice en $a_i$). 


::: {.conjecture name='"Robin Hood" de inicialización; Kronmal y Peterson, 1979' #robin-hood }
<br>

1.  Desde $i=1$ hasta $n$ hacer $q_{i}=np_{i}$.

2.  Establecer $L=\left\{ l:q_{l}<1\right\}$ y 
    $H=\left\{ h:q_{h}\geq 1\right\}$.

3.  Si $L$ ó $H$ vacios terminar.

4.  Seleccionar $l\in L$ y $h\in H$.

5.  Hacer $a_{l}=h$.

6.  Eliminar $l$ de $L$.

7.  Hacer $q_{h}=q_{h}-\left( 1-q_{l}\right)$.

8.  Si $q_{h}<1$ mover $h$ de $H$ a $L$.

9.  Ir al paso 3.

:::

\begin{figure}[!htb]

{\centering \includegraphics[width=0.9\linewidth]{05-Metodos_generales_discretas_files/figure-latex/unnamed-chunk-14-1} 

}

\caption{Pasos del algoritmo de inicialización del método Alias.}(\#fig:unnamed-chunk-14)
\end{figure}

El algoritmo para generar las simulaciones es el estándar del método de composición:

::: {.conjecture name="método alias de simulación; Walker, 1977" #walker}
<br>
  
1.  Generar $U,V\sim \mathcal{U}\left( 0,1\right)$.

2.  Hacer $i=\left\lfloor nU\right\rfloor +1$.

3.  Si $V<q_{i}$ devolver $X=x_{i}$.

4.  En caso contrario devolver $X=x_{a_{i}}$.

:::

Este algoritmo es muy eficiente y es el empleado en la función `sample()` de R^[R implementa este algoritmo en el fichero fuente [random.c](https://svn.r-project.org/R/trunk/src/main/random.c) (para muestreo probabilístico con reemplazamiento, función C `walker_ProbSampleReplace()`), aunque el paso 2 del algoritmo de simulación empleado por defecto cambió ligeramente a partir de la versión 3.6.0 para evitar posibles problemas de redondeo (ver Sección \@ref(oprrng)).].

Este método también está implementado en la función  [`simres::rpmf.alias()`](https://rubenfcasal.github.io/simres/reference/rpmf.alias.html) (fichero [*rpmf.R*](R/rpmf.R)), empleando código R menos eficiente:


```r
rpmf.alias
```

```
 ## function(x, prob = 1/length(x), n = 1000, as.factor = FALSE) {
 ##   # Inicializar tablas
 ##   a <- numeric(length(x))
 ##   q <- prob*length(x)
 ##   low <- q < 1
 ##   high <- which(!low)
 ##   low <- which(low)
 ##   while (length(high) && length(low)) {
 ##     l <- low[1]
 ##     h <- high[1]
 ##     a[l] <- h
 ##     q[h] <- q[h] - (1 - q[l])
 ##     if (q[h] < 1) {
 ##       high <- high[-1]
 ##       low[1] <- h
 ##     } else low <- low[-1]
 ##   } # while
 ##   # Generar valores
 ##   V <- runif(n)
 ##   i <- floor(runif(n)*length(x)) + 1
 ##   X <- x[ ifelse( V < q[i], i, a[i]) ]
 ##   if(as.factor) X <- factor(X, levels = x)
 ##   return(X)
 ## }
 ## <bytecode: 0x000001b7dc810628>
 ## <environment: namespace:simres>
```


::: {.example #binom-alias name="Simulación de una binomial mediante en método de Alias"}
<br>

Repetimos la simulación del Ejemplo \@ref(exm:binom-cuant) anterior empleando esta rutina.



```r
set.seed(1)
system.time( rx <- rpmf.alias(x, pmf, nsim) )
```

```
 ##    user  system elapsed 
 ##    0.02    0.00    0.02
```

Análisis de los resultados:

```r
res <- table(rx)/nsim
plot(res, ylab = "frecuencia relativa", xlab = "valores")
points(x, pmf, pch = 4, col = "blue")  # Comparación teórica
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{05-Metodos_generales_discretas_files/figure-latex/compalias-1} 

}

\caption{Comparación de las frecuencias relativas de los valores generados, mediante el método de alias, con las probabilidades teóricas.}(\#fig:compalias)
\end{figure}

:::


## Simulación de una variable discreta con dominio infinito

Los métodos anteriores están pensados para variables que toman un número finito de valores.
Si la variable discreta tiene dominio infinito no se podrían almacenar las probabilidades acumuladas, aunque en algunos casos podrían calcularse de forma recursiva.

::: {.example name="distribución de Poisson"}
<br>

Una variable $X$ con distribución de Poisson de parámetro $\lambda$, 
toma los valores $x_{1}=0$, $x_{2}=1$, $\ldots$ con probabilidades:
$$p_{j}=P\left( X=x_{j}\right)  =P\left( X=j-1\right)  =\frac{e^{-\lambda
}\lambda^{j-1}}{\left( j-1\right)  !}\text{, }j=1,2,\ldots$$
En este caso, como:
$$P\left( X=j\right)  =\frac{e^{-\lambda}\lambda^{j}}{j!}
=\frac{\lambda}{j}\frac{e^{-\lambda}\lambda^{j-1}}{\left( j-1\right)  !}
=\frac{\lambda}{j}P\left( X=j-1\right),$$
el algoritmo de inversión con búsqueda secuencial sería:

1. Generar $U\sim \mathcal{U}\left( 0,1\right)$.

2. Hacer $I=0$, $p=e^{-\lambda}$ y $S=p$.

3. Mientras $U>S$ hacer $I=I+1$, $p=\frac{\lambda}{I}p$ y $S=S+p$.

4. Devolver $X=I$.

:::

Hay modificaciones de los algoritmos anteriores, por ejemplo el de tabla guía con búsqueda secuencial en la cola de la distribución, para variables con dominio infinito.

Como alternativa, siempre se puede pensar en truncar la distribución,
eliminando los valores muy poco probables (teniendo en cuenta el número de generaciones que se pretenden realizar), 
aunque la distribución de las simulaciones será aproximada.


## Otros métodos

Muchos de los métodos descritos en el Capítulo \@ref(continuas) para variables continuas son directamente aplicables al caso discreto:

* Aceptación-Rechazo (Sección \@ref(AR)): En principio habría que considerar una variable auxiliar discreta con el mismo soporte, pero también hay modificaciones para variables auxiliares continuas.

* Método de composición (Sección \@ref(composicion)): es uno de los más empleados, por ejemplo en el método de Alias (Sección \@ref(alias)) y para simular la distribución binomial negativa (Sección \@ref(notables-disc)).

Hay otros métodos que tratan de reducir el número medio de comparaciones de la búsqueda secuencial, por ejemplo los árboles (binarios) de Huffman (e.g. Cao, 2002, Sección 4.2).
Estos métodos son muy poco eficientes para simular variables discretas pero pueden resultar de utilidad para diseñar experimentos de simulación más complejos (la idea es la misma, preocuparse principalmente por los sucesos más probables).

En ocasiones el método de la transformación cuantil puede acelerarse computacionalmente porque, mediante cálculos directos, es posible encontrar el valor de la función cuantil en cualquier $U$, evitando el bucle de búsqueda. 
Normalmente se realiza mediante truncamiento de una distribución continua.

::: {.example #transcuant-directo name="Cálculo directo de la función cuantil"}
<br>

La función de masa de probabilidad de una distribución geométrica es:
$$P\left( X=j\right)  =P\left( I=j+1\right)  =p\left( 1-p\right)^{j}\text{,
}j=0,1,\ldots$$

Si se considera como variable aleatoria continua auxiliar una exponencial, con función de distribución 
$G\left( x\right) = 1-e^{-\lambda x}$ si $x\geq0$,
se tiene que:
$$\begin{aligned}
G\left( i\right) - G\left( i-1\right)   
& = 1-e^{-\lambda i}-\left(1-e^{-\lambda\left( i-1\right) }\right)  
= e^{-\lambda\left( i-1\right)}-e^{-\lambda i}\\
& = e^{-\lambda\left( i-1\right)  }\left( 1-e^{-\lambda}\right)  
= \left( 1-e^{-\lambda}\right)  \left( e^{-\lambda}\right)^{i-1} \\
& = p\left(1-p\right)^{i-1},
\end{aligned}$$ 
tomando $p=1-e^{-\lambda}$.
De donde se obtendría el algoritmo:

:::

::: {.conjecture #geometrica name="distribución geométrica"}
<br>

0. Hacer $\lambda=-\ln\left( 1-p\right)$.

1. Generar $U\sim \mathcal{U}\left( 0,1\right)$.

2. Hacer $T=-\frac{\ln U}{\lambda}$.

3. Devolver $X=\left\lfloor T\right\rfloor$.

:::


## Métodos específicos para generación de distribuciones notables {#notables-disc}

Los comentarios al principio de la Sección \@ref(notables-cont) para el caso de variables continuas serían válidos también para distribuciones notables discretas.

Entre los distintos métodos disponibles para la generación de las distribuciones discretas más conocidas podríamos destacar el de la distribución binomial negativa mediante el método de composición (Sección \@ref(composicion)).

La distribución binomial negativa, $BN(r, p)$, puede interpretarse como el número de fracasos antes del $r$-ésimo éxito^[La distribución binomial negativa es una generalización de la geométrica y, debido a su reproductividad en el parámetro $r$, podría simularse como suma de $r$ variables geométricas. Sin embargo, este algoritmo puede ser muy costoso en tiempo de computación si $r$ es elevado.] y su función de masa de probabilidad es
$$P(X = i) = \binom{i+r-1}i p^r (1-p)^i \text{, para }i=0,1,\ldots$$

A partir de la propiedad
$$X|_{Y} \sim \text{Pois}\left(  Y\right)  \text{, }Y \sim \operatorname{Gamma} \left( r, \frac{p}{1-p}\right)  \Rightarrow X \sim BN(r, p)$$
se puede deducir el siguiente método específico de simulación.

::: {.conjecture #bin-neg name="distribución binomial negativa"}
<br>

1. Simular $L \sim \operatorname{Gamma}\left( r, \frac{p}{1-p} \right)$.

2. Simular $X \sim Pois \left(  L\right)$.

3. Devolver $X$.

:::

Por este motivo se denominada también a esta distribución *Gamma-Poisson*. Empleando una aproximación similar podríamos generar otras distribuciones, como la *Beta-Binomial*, empleadas habitualmente en inferencia bayesiana.


## Ejercicios {#ejercicios-discretas}

::: {.exercise #mixta-cuantil name="Simulación de una distribución mixta mediante el método de inversión generalizado"}
<br>

Consideramos la variable aleatoria con función de distribución dada por: 
$$F(x)=\left\{
\begin{array}
[c]{cl}0 & \mbox{si $x<0$}\\
\frac{x}{2}+\frac{1}{10} & \mbox{si $x\in[0,\frac{1}{5})$}\\
x+\frac{1}{10} & \mbox{si $x\in[\frac{1}{5},\frac{9}{10}]$}\\
1 & \mbox{en otro caso}
\end{array}
\right.$$

**Nota**: Esta variable toma los valores 0 y 1/5 con probabilidad 1/10.

Esta función está implementada en el siguiente código:

```r
fdistr <- function(x) {
ifelse(x < 0, 0,
    ifelse(x < 1/5, x/2 + 1/10,
        ifelse(x <= 9/10, x + 1/10, 1) ) )
}
```

Como es una función vectorial podemos emplear `curve()` para representarla:

```r
curve(fdistr, from = -0.1, to = 1.1, type = 's', main = '')
abline(h = c(1/10, 2/10, 3/10), lty = 2) # Discontinuidades
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{05-Metodos_generales_discretas_files/figure-latex/mixta-cuantil-plot-1} 

}

\caption{Función de distribución mixta (con discontinuidades en 0 y 1/5).}(\#fig:mixta-cuantil-plot)
\end{figure}

a)  Diseñar un algoritmo basado en el método de inversión generalizado 
    para generar observaciones de esta variable.
     
b)  Implementar el algoritmo en una función que permita generar $nsim$ 
    valores de esta variable.   
    
:::



Ver solución en Sección \@ref(sol-mixta-cuantil).


::: {.exercise #hipergeom name="distribución hipergeométrica"}
<br>

Se pretende simular $nsim=10^{4}$ observaciones de una variable
hipergeométrica (`dhyper(x, m, n, k)`) de parámetros $m$ los dos últimos dígitos del DNI, $n=100-m$ y $k=20$.

a)  Comprobar que el rango de posibles valores de esta variable es
    `max(0, k-n):min(m, k)`. Generar los valores empleando el método
    de la transformación cuantil usando búsqueda secuencial. Obtener
    el tiempo de CPU empleado. Aproximar por simulación la función
    de masa de probabilidad, representarla gráficamente y compararla
    con la teórica. Calcular también la media muestral (compararla
    con la teórica $km/(m+n)$) y el número medio de comparaciones
    para generar cada observación.

b)  Repetir el apartado anterior: ordenando previamente las
    probabilidades en orden decreciente, empleando la
    función `sample` de R, mediante una tabla guía (con
    $k-1$ subintervalos) y usando el método de Alias.

:::


