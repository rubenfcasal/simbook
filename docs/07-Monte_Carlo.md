# Métodos Monte Carlo {#monte-carlo}

<!-- Capítulo \@ref(monte-carlo) -->



<!-- 
---
title: "Métodos Monte Carlo"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    pandoc_args: ["--number-offset", "7,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("07-Monte_Carlo.Rmd")
knitr::purl("07-Monte_Carlo.Rmd", documentation = 2)
knitr::spin("07-Monte_Carlo.R",knit = FALSE)

Pendiente:
Métodos Monte Carlo en inferencia estadística
  Sección Estudios de simulación
Temple simulado
  Tomar como referencia el código de optim para SANN
  https://svn.r-project.org/R/trunk/src/appl/optim.c
  Cambiar gráfica
-->


Uno de los objetivos habituales en inferencia es la aproximación de una esperanza (o el caso particular de una probabilidad), es decir, se trataría de evaluar una integral, que en ocasiones puede ser compleja. 
Esto puede ser de interés en otros campos, aunque la integral no esté relacionada con procesos estocásticos. 
Si el número de dimensiones es pequeño puede ser recomendable emplear métodos numéricos (ver Apéndice \@ref(int-num)), pero si la dimensión del dominio de la integral es grande, puede ser mucho más eficiente emplear métodos basados en simulación o incluso ser la única aproximación realmente viable.
En las secciones \@ref(int-MC) y \@ref(muestreo-importancia) se describen este tipo de procedimientos, denominados métodos de integración Monte Carlo.

Los métodos Monte Carlo emplean simulación para resolver problemas matemáticos, como la evaluación de expresiones, la aproximación de integrales o la resolución de sistemas, entre muchos otros.
Estos problemas son de interés en muchos campos (Física, Economía, Informática...) y pueden ser estrictamente deterministas.
Otro ejemplo es la resolución de problemas de optimización. 
Para evitar problemas de mínimos locales se puede recurrir a herramientas que emplean búsquedas aleatorias de los valores óptimos. 
En la Sección \@ref(opt-MC) se describen brevemente los principales métodos de optimización Monte Carlo.

En muchos casos, especialmente en Inferencia Estadística, existe una aleatoriedad inherente al modelo empleado para resolver el problema y los métodos Monte Carlo surgen de forma natural.
Como se comentó en la Sección \@ref(conceptos), en ocasiones no se pueden obtener soluciones analíticas a problemas de inferencia, o solo se disponen de resultados asintóticos que pueden no ser suficientemente buenos para muestras finitas, y se puede recurrir a la simulación.

Los estudios Monte Carlo son una herramienta habitual para investigar las ventajas e inconvenientes de un método de inferencia, y para entender mejor su funcionamiento. 
Por este motivo suelen ser el primer paso en el desarrollo de métodos Monte Carlo (que posteriormente pueden ser objeto de estudios analíticos si producen buenos resultados).
También se pueden emplear para tratar de determinar, entre los métodos disponibles, el más adecuado para resolver el problema de interés.

En Estadística Computacional (que incluiría los métodos de Aprendizaje Estadístico/Automático) se emplean métodos de inferencia computacionalmente intensivos y muchos de ellos emplean simulación. 
Entre las técnicas empleadas destacan los métodos de remuestreo, como el jackknife o el bootstrap, que trataremos en capítulos siguientes.
Como ya se comentó, la simulación empleando un modelo estimado también se denomina bootstrap paramétrico y la mayoría de los métodos Monte Carlo de inferencia estadística los podríamos clasificar como métodos de remuestreo.

<!-- 
Bootstrap paramétrico 
Eficiencia experimento simulación -> técnicas de reducción de la varianza
-->

## Integración Monte Carlo {#int-MC}

La integración Monte Carlo se emplea principalmente para aproximar integrales multidimensionales:
$$I = \int \cdots \int _D s\left( x_1,\ldots ,x_n\right) dx_1 \cdots dx_n$$ 
donde puede presentar ventajas respecto a los métodos tradicionales de integración numérica (ver Apéndice \@ref(int-num)), ya que la velocidad de convergencia no depende del número de dimensiones.

La idea es reescibir la expresión de la integral, encontrando una función de densidad $f$ definida sobre $D$, de forma que:
$$I = \int _D s(\mathbf{x}) d \mathbf{x} = \int h(\mathbf{x})f(\mathbf{x}) d \mathbf{x} = E\left( h(\mathbf{X}) \right)$$ 
donde $\mathbf{X} \sim f$ (y preferiblemente fácil de simular).

### Integración Monte Carlo clásica

En el caso de que el dominio $D$ sea acotado, la aproximación más simple consiste en considerar una distribución uniforme en $D$ (i.e. $f(\mathbf{x})=1_D(\mathbf{x})/|D|$ y $h(\mathbf{x}) = |D|s(\mathbf{x})$).

Por simplicidad nos centraremos en el caso unidimensional (el orden de convergencia es independiente del número de dimensiones).
Supongamos que nos interesa aproximar:
$$I = \int_0^1 s(x) dx$$
Si $x_1,x_2,\ldots ,x_n$ *i.i.d.* $\mathcal{U}(0, 1)$
entonces:
$$I = E\left( s\left( \mathcal{U}(0, 1) \right) \right)
\approx \frac{1}{n}\sum\limits_{i=1}^n s\left( x_i\right)$$

Si el intervalo de integración es genérico:
$$I = \int_a^b s(x) dx = 
(b-a)\int_a^b s(x) \frac1{(b-a)}dx = 
(b-a)E\left( s\left( \mathcal{U}(a, b) \right) \right).$$
Si $x_1,x_2,\ldots ,x_n$ *i.i.d.* $\mathcal{U}(a, b)$:
$$I \approx \frac{b-a}{n}\sum\limits_{i=1}^n s\left( x_i\right)$$


::: {.example #mc-integral-clas name="integración Monte Carlo clásica"}
<br>

Como primera aproximación para implementar la integración Monte Carlo clásica
para aproximar integrales del tipo:
$$I = \int_a^b s(x)  dx,$$ 
podríamos considerar la siguiente función:


``` r
mc.integral0 <- function(fun, a, b, n) {
  # Integración Monte Carlo de `fun()` entre `a` y `b` utilizando una muestra 
  # (pseudo) aleatoria de tamaño `n`. Se asume que `fun()` es una función de 
  # una sola variable (no vectorial), `a < b` y `n` entero positivo.
  # -----------------------  
  x <- runif(n, a, b)
  fx <- sapply(x, fun) # Si fun fuese vectorial bastaría con: fx <- fun(x)
  return(mean(fx) * (b - a))
}
```

Como ejemplo la empleamos para aproximar: 
$$\int_0^1 4x^4 dx = \frac{4}{5},$$


``` r
fun <- function(x) ifelse((x > 0) & (x < 1), 4 * x^4, 0)
# return(4 * x^4)
curve(fun, 0, 1)
abline(h = 0, lty = 2)
abline(v = c(0, 1), lty = 2)

set.seed(1)
mc.integral0(fun, 0, 1, 20)
```

```
 ## [1] 0.97766
```

``` r
mc.integral0(fun, 0, 1, 100)
```

```
 ## [1] 0.73112
```

``` r
mc.integral0(fun, 0, 1, 100)
```

```
 ## [1] 0.83043
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/int-mc-clas-1} 

}

\caption{Ejemplo de integral en dominio acotado.}(\#fig:int-mc-clas)
\end{figure}

La función `mc.integral0()` no es adecuada para analizar la convergencia de la aproximación por simulación.
Una alternativa más eficiente para representar gráficamente la convergencia está implementada en la función [`mc.integral()`](https://rubenfcasal.github.io/simres/reference/mc.integral.html) del paquete [`simres`](https://rubenfcasal.github.io/simres) (fichero [*mc.plot.R*](R/mc.plot.R)):


``` r
library(simres)
mc.integral
```

```
 ## function(fun, a, b, n, level = 0.95, plot = TRUE, ...) {
 ##   fx <- sapply(runif(n, a, b), fun) * (b - a)
 ##   result <- if (plot) conv.plot(fx, level = level, ...) else {
 ##     q <- qnorm((1 + level)/2)
 ##     list(approx = mean(fx), max.error = q * sd(fx)/sqrt(n))
 ##   }
 ##   return(result)
 ## }
 ## <bytecode: 0x000002015b488830>
 ## <environment: namespace:simres>
```

``` r
set.seed(1)
mc.integral(fun, 0, 1, 5000, ylim = c(0.2, 1.4))
```

```
 ## $approx
 ## [1] 0.81422
 ## 
 ## $max.error
 ## [1] 0.030282
```

``` r
abline(h = 4/5, lty = 2, col = "blue")
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/mc-integral-1} 

}

\caption{Convergencia de la aproximación de la integral mediante simulación.}(\#fig:mc-integral)
\end{figure}

Si sólo interesa la aproximación:


``` r
set.seed(1)
mc.integral(fun, 0, 1, 5000, plot = FALSE)
```

```
 ## $approx
 ## [1] 0.81422
 ## 
 ## $max.error
 ## [1] 0.030282
```

**Nota**: Es importante tener en cuenta que la función [`mc.integral()`](https://rubenfcasal.github.io/simres/reference/mc.integral.html) solo es válida para dominio finito.

:::

### Caso general

En lo que resta de esta sección (y en las siguientes) asumiremos que nos interesa aproximar una esperanza:
$$\theta = E\left( h\left( X\right) \right) = \int h\left( x\right) f(x)dx$$
siendo $X\sim f$. 
Entonces, si $x_1,x_2,\ldots ,x_n$ *i.i.d.* $X$:
$$\theta \approx \frac{1}{n}\sum\limits_{i=1}^nh\left( x_i\right)$$

Por ejemplo, como en el ejercicio anterior se considera de una función de densidad, 
se correspondería con el caso general de $h(x) = x$ y $f(x) = 4x^3$ para $0<x<1$.
La idea es que, en lugar de considerar una distribución uniforme, 
es preferible generar más valores donde hay mayor "área" (ver Figura \@ref(fig:int-mc-clas)).

Los pasos serían simular `x` con densidad $f$ y aproximar la integral por `mean(h(x))`.
En este caso podemos generar valores de la densidad objetivo fácilmente mediante el método de inversión, ya que $F(x) = x^4$ si $0<x<1$.


``` r
rfun <- function(nsim) runif(nsim)^(1/4) # Método de inversión
nsim <- 5000
set.seed(1)
x <- rfun(nsim)
# h <- function(x) x
# res <- mean(h(x)) # Aproximación por simulación 
res <- mean(x)
res
```

```
 ## [1] 0.79678
```

``` r
# error <- 2*sd(h(x))/sqrt(nsim)
error <- 2*sd(x)/sqrt(nsim)
error
```

```
 ## [1] 0.0047282
```

Esta forma de proceder permite aproximar integrales impropias en las que el dominio de integración no es acotado.

::: {.example #mc-intinf name="integración Monte Carlo con dominio no acotado"}
<br>

Aproximar:
$$\phi(t)=\int_{t}^{\infty}\frac1{\sqrt{2\pi}}e^{-\frac{x^2}2}dx,$$
para $t=4.5$, empleando integración Monte Carlo (aproximación tradicional con dominio infinito).


``` r
# h <- function(x) x > 4.5
# f <- function(x) dnorm(x)
set.seed(1)
nsim <- 10^3
x <- rnorm(nsim)
mean(x > 4.5) # mean(h(x))
```

```
 ## [1] 0
```

``` r
pnorm(-4.5)  # valor teórico P(X > 4.5) 
```

```
 ## [1] 3.3977e-06
```

De esta forma es difícil que se generen valores (en este caso ninguno) en la región que interesaría para la aproximación de la integral:


``` r
any(x > 4.5)
```

```
 ## [1] FALSE
```

Como ya se comentó anteriormente, sería preferible generar más valores donde hay mayor "área", pero en este caso $f$ concentra la densidad en una región que no resulta de utilidad.
Por ese motivo puede ser preferible recurrir a una densidad auxiliar que solvente este problema.

:::


## Muestreo por importancia {#muestreo-importancia}

Para aproximar la integral:
$$\theta = E\left( h\left( X\right) \right) = \int h\left( x\right) f(x)dx,$$
puede ser preferible generar observaciones de una densidad $g$ que tenga una forma similar al producto $hf$.

Si $Y\sim g$:
$$\theta  = \int h\left( x\right) f(x)dx 
 = \int \frac{h\left( x\right) f(x)}{g(x)}g(x)dx
 = E\left( q\left( Y\right) \right).$$
siendo
$q\left( x\right)  = \frac{h\left( x\right) f(x)}{g(x)}$.

Si $y_1,y_2,\ldots ,y_n$ *i.i.d.* $Y\sim g$:
$$\theta \approx \frac{1}{n}\sum\limits_{i=1}^nq\left( y_i\right) 
= \frac{1}{n}\sum\limits_{i=1}^nw(y_i)h\left( y_i\right)  
= \hat{\theta}_{g}$$
con $w(y) = \frac{f(y)}{g(y)}$.

En este caso $Var(\hat{\theta}_{g}) = Var\left( q\left( Y\right) \right) /n$,  pudiendo reducirse significativamente respecto al método clásico si:
$$g(x)\underset{aprox.}{\propto } \left\vert h(x) \right\vert f(x),$$
ya que en ese caso $\left\vert q(x) \right\vert$ sería aproximadamente constante (puede demostrarse fácilmente que la varianza es mínima si esa relación es exacta).

Para garantizar la convergencia de la aproximación por simulación, la varianza del estimador $\hat{\theta}_{g}$ debería ser finita, i.e.:
$$E\left( q^2\left( Y\right) \right)  
= \int \frac{h^2\left( x\right)f^2(x)}{g(x)}dx 
= E\left( h^2\left( X\right) \frac{f(X)}{g(X)}\right)
< \infty.$$
La idea básica es que si la densidad $g$ tiene colas más pesadas que la densidad $f$ con mayor facilidad puede dar lugar a "simulaciones" con varianza finita (podría emplearse en casos en los que no existe $E \left( h^2 \left( X \right) \right)$; ver Sección \@ref(convergencia)).

La distribución de los pesos $w(y_i)$ debería ser homogénea para evitar datos influyentes (que introducirían inestabilidad en la aproximación).


::: {.example #mc-imp}
<br>

Podríamos aproximar la integral del Ejemplo \@ref(exm:mc-intinf) anterior empleando muestreo por importancia considerando como densidad auxiliar una exponencial de parámetro $\lambda=1$ truncada en $t$:
$$g(x) = \lambda e^{-\lambda (x - t)}\text{, }x>t,$$
(podemos emplear `dexp(y - t)` para evaluar esta densidad y `rexp(n) + t` para generar valores). 

En primer lugar comparamos $h(x)f(x)$ con la densidad auxiliar reescalada, $g(x)f(4.5)$, para comprobar si es una buena elección:

``` r
curve(dnorm(x), 4.5, 6, ylab = "dnorm(x) y dexp(x-4.5)*k")
abline(v = 4.5)
abline(h = 0)
escala <- dnorm(4.5)  # Reescalado para comparación...
curve(dexp(x - 4.5) * escala, add = TRUE, lty = 2)  
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/mc-imp-int-1} 

}

\caption{Objetivo a integrar (densidad objetivo truncada) y densidad auxiliar reescalada.}(\#fig:mc-imp-int)
\end{figure}

Generamos valores de la densidad auxiliar y calculamos los pesos:

``` r
set.seed(1)
nsim <- 10^3
y <- rexp(nsim) + 4.5    #  Y ~ g
w <- dnorm(y)/dexp(y - 4.5)
```

La aproximación por simulación sería `mean(w * h(y))`:

``` r
# h(x) <- function(x) x > 4.5  # (1 si x > 4.5 => h(y) = 1)
mean(w) # mean(w*h(y))
```

```
 ## [1] 3.1448e-06
```

``` r
pnorm(-4.5)  # valor teórico
```

```
 ## [1] 3.3977e-06
```

Representamos gráficamente la aproximación en función del número de simulaciones:

``` r
plot(cumsum(w)/1:nsim, type = "l", ylab = "Aproximación", xlab = "Iteraciones")
abline(h = pnorm(-4.5), lty = 2, col = "blue")
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/mc-imp-conv-1} 

}

\caption{Convergencia de la aproximación de la integral mediante muestreo por importancia.}(\#fig:mc-imp-conv)
\end{figure}

El error estándar de la aproximación sería `sqrt(var(w * h(y))/nsim)`:

``` r
sqrt(var(w)/nsim) # sd(w*h(y))/sqrt(nsim)   
```

```
 ## [1] 1.3712e-07
```
Mientras que empleando la aproximación tradicional:

``` r
est <- mean(rnorm(nsim) > 4.5)
est
```

```
 ## [1] 0
```

``` r
sqrt(est * (1 - est)/nsim)
```

```
 ## [1] 0
```
:::


::: {.example #mc-imp2 name="muestreo por importancia con mala densidad auxiliar"}
<br>

Supongamos que se pretende aproximar $P\left(2<X<6\right)$ siendo $X\sim Cauchy(0,1)$ empleando muestreo por importancia y considerando como densidad auxiliar la normal estándar $Y\sim N(0,1)$. Representaremos gráficamente la aproximación y estudiaremos los pesos $w(y_i)$.
    
**Nota**: En este caso van a aparecer problemas 
(la densidad auxiliar debería tener colas más pesadas que la densidad objetivo;
sería adecuado si intercambiáramos las distribuciones objetivo y auxiliar,
como en el Ejemplo \@ref(exm:mc-imp-sample) siguiente).

Se trata de aproximar `pcauchy(6) - pcauchy(2)`,
i.e. `f(y) = dcauchy(y)` y `h(y) = (y > 2) * (y < 6)`,
empleando muestreo por importancia con `g(y) = dnorm(y)`.

``` r
nsim <- 10^5
set.seed(4321)
y <- rnorm(nsim)
w <- dcauchy(y)/dnorm(y) # w <- w/sum(w) si alguna es una cuasidensidad
```

La aproximación por simulación es `mean(w(y) * h(y))`:

``` r
mean(w * (y > 2) * (y < 6)) 
```

```
 ## [1] 0.099293
```

``` r
pcauchy(6) - pcauchy(2)  # Valor teórico
```

```
 ## [1] 0.095015
```

Si se estudia la convergencia:

``` r
plot(cumsum(w * (y > 2) * (y < 6))/1:nsim, type = "l", ylab = "Aproximación", xlab = "Iteraciones")
abline(h = pcauchy(6) - pcauchy(2), lty = 2, col = "blue")
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/mc-imp2-conv-1} 

}

\caption{Gráfico de convergencia de la aproximación mediante muestreo por importancia con mala densidad auxiliar.}(\#fig:mc-imp2-conv)
\end{figure}
Lo que indica es una mala elección de la densidad auxiliar.

La distribución de los pesos debería ser homogénea.
Por ejemplo, si los reescalamos para que su suma sea el número de valores generados, deberían tomar valores en torno a uno:

``` r
boxplot(nsim * w/sum(w))  
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/mc-imp2-boxplot-1} 

}

\caption{Gráfico de cajas de los pesos del muestreo por importancia reescalados (de forma que su media es 1).}(\#fig:mc-imp2-boxplot)
\end{figure}

:::


### Remuestreo (del muestreo) por importancia

Cuando $f$ y/o $g$ son cuasi-densidades, para evitar calcular constantes normalizadoras, se emplea como aproximación:
$$\theta \approx \frac{\sum\limits_{i=1}^n w(y_i) h\left( y_i\right) }{ \sum\limits_{i=1}^n w(y_i)}.$$
De hecho este estimador es empleado muchas veces en lugar del anterior ya que, aunque en general no es insesgado, puede ser más eficiente si $w(Y)$ y $w(Y)h(Y)$ están altamente correlacionadas [e.g. @liu2004, p.35].

Adicionalmente, puede verse que con un muestreo de $\left\{y_1, y_2, \ldots, y_n \right\}$ ponderado por $w(y_i)$ (prob. $=w(y_i)\left/ \sum\nolimits_{i=1}^n w(y_i) \right.$ ) se obtiene una simulación aproximada de $f$ [*Sample importance resampling*, @rubin1987].


::: {.example #mc-imp-sample name="simulación de normal estándar a partir de Cauchy; Sampling Importance Resampling"}
<br>

Generamos 1000 simulaciones de una distribución (aprox.) $N(0,1)$ (densidad objetivo) mediante remuestreo del muestreo por importancia de $10^{5}$ valores de una $Cauchy(0,1)$ (densidad auxiliar).
    
**Nota**: En este caso `f(y) = dnorm(y)` y `g(y) = dcauchy(y)`, al revés del Ejemplo \@ref(exm:mc-imp2) anterior.


``` r
# Densidad objetivo
# f <- dnorm # f <- function(x) ....

nsim <- 10^3
# El nº de simulaciones de la densidad auxiliar debe ser mucho mayor:
nsim2 <- 10^5
set.seed(4321)
y <- rcauchy(nsim2)
w <- dnorm(y)/dcauchy(y) # w <- w/sum(w) si alguna es una cuasidensidad

# Si se pidiera aproximar una integral
# h(y) = y si es la media # h <- function(y) y
# mean(w * h(y))
```

Sampling Importance Resampling: 


``` r
rx <- sample(y, nsim, replace = TRUE, prob = w/sum(w))
hist(rx, freq = FALSE, breaks = "FD", ylim = c(0, 0.5))
lines(density(rx))
curve(dnorm, col = "blue", add = TRUE)
```

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/imp-res-1} 

}

\caption{Distribución de los valores generados mediante remuestreo por importancia y densidad objetivo.}(\#fig:imp-res)
\end{figure}

:::

**Nota**: Si f o g fuesen cuasidensidades y se pidiese aproximar la integral, habría que reescalar los pesos  `w <- f(y)/g(y)` en la aproximación por simulación, resultando `sum(w * h(y))/sum(w)` (media ponderada) y en el análisis de convergencia se emplearía `cumsum(w * h(y))/cumsum(w)`.

::: {.exercise #mc-imp-sample2 name="propuesto"}
<br>

Consideramos una variable aleatoria con densidad:
    $$f(x)\propto e^{-x}\cos^{2}(x),\text{ si }x>0.$$

a.  Aproximar mediante integración Monte Carlo la media de esta
    distribución ($h(x)=x$) empleando muestreo de importancia con
    distribución auxiliar una exponencial de parámetro $\lambda=1$ a
    partir de 10000 simulaciones (OJO: se conoce la cuasi-densidad
    de la variable aleatoria de interés, emplear la aproximación descrita
    en apuntes).

b.  Generar 500 simulaciones (aprox.) de la distribución de interés
    mediante remuestreo del muestreo por importancia.
    
**NOTA**: En el último apartado, para comprobar que los valores generados proceden de la distribución objetivo, si representamos la cuasidensidad $f^{\ast}(x) = e^{-x}\cos^{2}(x)$ junto con el histograma (en escala de densidades, `freq = FALSE`), hay que tener en cuenta que faltaría dividir la cuasidensidad por una constante normalizadora para poder compararlos directamente. 
Si no se reescala la cuasidensidad, podríamos compobar si la forma es similar (si la distribución de los valores generados es proporcional a la cuasidensidad, con mayor concentración donde la cuasidensidad se aleja de 0). 
En este caso (como $g$ es una densidad) podríamos estimar la constante normalizadora ($f(x) = \frac{1}{c}f^{\ast}(x)$) a partir de los pesos del muestreo por importancia (`c.approx <- mean(w)`; en este caso concreto $c=\frac{3}{5}$).

:::


<!-- 
---  

***LA MATERIA EVALUABLE EN EL CURSO 2022/2023 TERMINA AQUÍ***

--- 
-->

## Optimización Monte Carlo {#opt-MC}

Supongamos que estamos interesados en la minimización de una función:
$$\underset{\mathbf{x}\in D}{\min }f(\mathbf{x}).$$

Hay una gran cantidad de algoritmos numéricos para resolver problemas de optimización no lineal multidimensional, por ejemplo los basados en el método de Newton-Raphson (implementados en la función `nlm`, entre otras).

La idea original consiste en buscar los ceros de su primera derivada (o del gradiente) empleando una aproximación iterativa:
$$\mathbf{x}_{i+1} = \mathbf{x}_i-[Hf(\mathbf{x}_i)]^{-1}\nabla f(\mathbf{x}_i),$$
donde $Hf(\mathbf{x}_i)$ es el hessiano de la función (matriz de segundas derivadas) y $\nabla f(\mathbf{x}_i)$ el gradiente (vector de primeras derivadas).
Estos métodos normalmente funcionan muy bien cuando la función objetivo no tiene mínimos locales (ideal $f$ cuadrática).
Los resultados obtenidos pueden ser muy malos en caso contrario (especialmente en el caso multidimensional) y dependen en gran medida del punto inicial^[Este tipo de algoritmos se denominan *codiciosos* o *voraces*, porque buscan la mejor opción a "corto plazo".] 
Un ejemplo donde es habitual que aparezcan este tipo de problemas es en la estimación por máxima verosimilitud (la función objetivo puede ser multimodal).


::: {.example #mv-nlm name="Estimación por máxima verosimilitud mediante un algoritmo de Newton"}
<br>

La mixtura de distribuciones normales:
$$\frac1{4}N(\mu_1,1)+\frac{3}{4}N(\mu_2,1),$$ 
tiene una función de verosimilitud asociada bimodal.
Generaremos una muestra de 200 valores de esta distribución con $\mu_1=0$ y $\mu_2=2.5$, construiremos la correspondiente función de verosimilitud y la representaremos gráficamente.

Obtención de la muestra (simulación mixtura dos normales):


``` r
nsim <- 200
mu1 <- 0 
mu2 <- 2.5
sd1 <- sd2 <- 1

set.seed(12345)
p.sim <- rbinom(nsim, 1, 0.25)
data <- rnorm(nsim, mu1*p.sim + mu2*(1-p.sim), sd1*p.sim + sd2*(1-p.sim))

hist(data, freq = FALSE, breaks = "FD", ylim = c(0, 0.3))
curve(0.25 * dnorm(x, mu1, sd1) + 0.75 * dnorm(x, mu2, sd2), add = TRUE)
```



\begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-13-1} \end{center}

Podemos obtener la estimación por máxima verosimilitud de los parámetros empleando la rutina `nlm` para minimizar el logaritmo (negativo) de la función de verosimilitud:

``` r
like <- function(mu)
  -sum(log((0.25 * dnorm(data, mu[1], sd1) + 0.75 * dnorm(data, mu[2], sd2))))
  # NOTA: Pueden aparecer NA/Inf por log(0)
```

Si queremos capturar los valores en los que se evalúa esta función, podemos proceder de forma similar a como se describe en el capítulo 
[Function operators](http://adv-r.had.co.nz/Function-operators.html#behavioural-fos) 
de la primera edición del libro "Advanced R" de Hadley Wickham:
"Behavioural FOs leave the inputs and outputs of a function unchanged, 
but add some extra behaviour". 

``` r
tee <- function(f) {
  function(...) {
    input <- if(nargs() == 1) c(...) else list(...)
    output <- f(...)
    # Hacer algo ...
    # ... con output e input
    return(output)
  }
}
```

En este caso queremos representar los puntos en los que el algoritmo de optimización evalúa la función objetivo (especialmente como evoluciona el valor óptimo) 


``` r
tee.optim2d <- function(f) {
  best.f <- Inf   # Suponemos que se va a minimizar (opción por defecto)
  best.par <- c(NA, NA)   
  function(...) {
    input <- c(...)
    output <- f(...)
    ## Hacer algo ...
    points(input[1], input[2], col = "lightgray")
    if(best.f > output) {
      lines(rbind(best.par, input), lwd = 2, col = "blue")
      best.f <<- output
      best.par <<- input
      # points(best.par[1], best.par[2], col = "blue", pch = 20)
      # cat("par = ", best.par, "value = ", best.f, "\n")
    } 
    ## ... con output e input
    return(output)
  }
}
```

Representar la superficie del logaritmo de la verosimilitud, 
los puntos iniciales y las iteraciones en la optimización numérica con `nlm`:


``` r
mmu1 <- mmu2 <- seq(-2, 5, length = 128)
lli <- outer(mmu1, mmu2, function(x,y) apply(cbind(x,y), 1, like))

par(mar = c(4, 4, 1, 1))
image(mmu1, mmu2, -lli, xlab = expression(mu[1]), ylab = expression(mu[2]))
contour(mmu1, mmu2, -lli, nlevels = 50, add = TRUE)

# Valores iniciales aleatorios
nstarts <- 5
set.seed(1)
starts <- matrix(runif(2*nstarts, -2, 5), nrow = nstarts)
points(starts, col = "blue", pch = 19)

# Minimización numérica con nlm
for (j in 1:nstarts){
  # Normalmente llamaríamos a nlm(like, start)
  res <- nlm(tee.optim2d(like), starts[j, ]) # nlm(like, starts[j, ])
  points(res$estimate[1],res$estimate[2], pch = 19)
  cat("par = ", res$estimate, ", value =", res$minimum, "\n")
}
```

```
 ## par =  -0.038925 2.4946 , value = 361.57 
 ## par =  -0.038925 2.4946 , value = 361.57
```

```
 ## par =  -0.038925 2.4946 , value = 361.57 
 ## par =  3.1322 0.96285 , value = 379.37
```

```
 ## par =  20.51 1.712 , value = 474.14
```



\begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-17-1} \end{center}

:::


### Algoritmos de optimización Monte Carlo

Una alternativa sería tratar de generar valores aleatoriamente de forma que las regiones donde la función objetivo es menor tuviesen mayor probabilidad y menor probabilidad las regiones donde la función objetivo es mayor.
Por ejemplo, se podría pensar en generar valores de acuerdo a una densidad (tranformación Boltzman-Gibbs):
$$g(x)\propto \exp \left( -f(x)/T\right) ,$$
donde $T>0$ es un parámetro (denominado temperatura) seleccionado de forma que se garantice la integrabilidad.

Entre los métodos de optimización Monte Carlo podríamos destacar:

-   Métodos con gradiente aleatorio.
    <!-- 
    Ver paquete https://cran.r-project.org/package=maxLik
    Henningsen A, Toomet O (2011). maxLik: A package for maximum likelihood estimation in R. Computational Statistics, 26(3), 443-458. https://doi.org/10.1007/s00180-010-0217-1. 
    -->

-   Temple simulado.

-   Algoritmos genéticos.

-   Monte Carlo EM.

-   ...

### Temple simulado

Método inspirado en el templado de un metal (se calienta el metal a alta temperatura y se va enfriando lentamente).
En cada paso se reemplaza la aproximación actual por un valor aleatorio “cercano”, elegido con una probabilidad que depende de la mejora en la función objetivo y de un parámetro $T$ (denominado temperatura) que disminuye gradualmente durante el proceso.

-   Cuando la temperatura es grande los cambios son bastante
    probables en cualquier dirección.

-   Al ir disminuyendo la temperatura los cambios tienden a ser
    siempre “cuesta abajo”.

Al tener una probabilidad no nula de aceptar una modificación “cuesta arriba” se trata de evitar quedar atrapado en un óptimo local (ver Figura \@ref(fig:templesimulado)).

\begin{figure}[!htbp]

{\centering \includegraphics[width=0.75\linewidth]{images/templesimulado} 

}

\caption{Fuente: Premchand Akella ([ppt](https://www.presentica.com/doc/11473134/simulated-annealing-pdf-document)).}(\#fig:templesimulado)
\end{figure}

Este procedimiento se puede ver como una adaptación del método de Metropolis-Hastings que se tratará en el Capítulo XX (Introducción a los métodos de cadenas de Markov Monte Carlo).


***Algoritmo***


``` r
temp <- temp.ini
par <- par.ini
fun.par <- FUN(par)
iter <- 1
while(temp > temp.final && iter < iter.max) {
  iter.temp <- 1
  while(iter.temp < iter.temp.max)) { # iteraciones con temperatura constante
    par.new <- PERTURB(par, temp)
    fun.new <- FUN(par.new)
    fun.inc <- fun.new - fun.par
    if ((fun.inc < 0) || (runif(1) > exp(-(fun.inc/temp)))) break
    iter.temp <- iter.temp + 1
  }
  iter <- iter + iter.temp
  par <- par.new
  fun.par <- fun.new
  temp <- SCHEDULE(temp)
}

FUN <- function(par, ...) {...}
SCHEDULE <- function(temp, temp.ini, iter) 
  temp.ini / log(iter + exp(1) - 1)
  # temp.ini / log(((temp - 1) %/% tmax)*tmax + exp(1))
PERTURB <- function(par, temp, scale = 1/temp.ini) 
  rnorm(length(par), par, 1/(scale*temp))
```

Una versión de este método está implementado^[En el fichero fuente [optim.c](https://svn.r-project.org/R/trunk/src/appl/optim.c).] en la función `optim()`:

``` r
optim(par, fn, gr = NULL, ..., method = "SANN", control = list(maxit = 10000, temp = 10, tmax = 10)
```

El argumento `gr` permite especificar la función para generar posiciones candidatas (por defecto núcleo gausiano con escala proporcional a la temperatura actual) y permitiría resolver problemas de optimización combinatoria.
El argumento `control` permite establecer algunas opciones adicionales:

- `maxit`: número total de evaluaciones de la función (único criterio de parada), por defecto 10000.  

- `temp`: temperatura inicial, por defecto 10.

- `tmax`: número de evaluaciones de la función para cada temperatura, por defecto 10.

::: {.example #mv-sann name="Estimación máximo-verosimil empleando temple simulado"}
<br>

Repetimos el Ejemplo \@ref(exm:mv-nlm) anterior empleando el método "SANN" de la función `optim()`:


``` r
# Representar la superficie del logaritmo de la verosimilitud
image(mmu1, mmu2, -lli, xlab = expression(mu[1]), ylab = expression(mu[2]))
contour(mmu1, mmu2, -lli, nlevels = 50, add = TRUE)
points(starts, col = "blue", pch = 19)

set.seed(1)
for (j in 1:nstarts){
  # Normalmente llamaríamos a optim(start, like, method = "SANN")
  # Note that the "SANN" method depends critically on the settings of the control parameters.
  # For "SANN" maxit gives the total number of function evaluations: there is no other stopping criterion. 
  # Defaults to 10000.
  res <- optim(starts[j, ], tee.optim2d(like), method = "SANN", control = list(temp = 100, maxit = 2000))
  points(res$par[1],res$par[2], pch = 19)
  cat("par = ", res$par, ", value =", res$value, "\n")
}
```

```
 ## par =  0.00020235 2.4734 , value = 361.64 
 ## par =  -0.18274 2.4559 , value = 362.03 
 ## par =  -0.028134 2.4845 , value = 361.58 
 ## par =  -0.036429 2.4886 , value = 361.57 
 ## par =  0.68142 2.37 , value = 374.84
```



\begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-20-1} \end{center}


Como alternativa podríamos emplear la siguiente función basada en el algoritmo del Ejemplo 5.9 de Robert y Casella (2010):


``` r
SA <- function(fun, pini, lower = -Inf, upper = Inf, tolerance = 1e-04, factor = 1) {
  temp <- scale <- iter <- dif <- 1
  npar <- length(pini) 
  par <- matrix(pini, ncol = npar)
  curfun <- hval <- fun(pini)
  while (dif > tolerance) {
    prop <- par[iter, ] + rnorm(npar) * scale[iter]
    # Se decide si se acepta la propuesta
    if (any(prop < lower) || any(prop > upper) || 
        (temp[iter] * log(runif(1)) > curfun - fun(prop)))  prop <- par[iter, ]
    curfun <- fun(prop)
    hval <- c(hval, curfun)
    par <- rbind(par, prop)
    iter <- iter + 1
    temp <- c(temp, 1/log(iter + 1))  # Actualizar la temperatura
    # Se controla el número de perturbaciones aceptadas
    ace <- length(unique(par[(iter/2):iter, 1]))
    if (ace == 1) 
      # si es muy pequeño se disminuye la escala de la perturbación
      factor <- factor/10
    if (2 * ace > iter) 
      # si es muy grande se aumenta
      factor <- factor * 10
    scale <- c(scale, max(2, factor * sqrt(temp[iter])))  # Actualizar la escala de la perturbación
    dif <- (iter < 100) + (ace < 2) + (max(hval) - max(hval[1:(iter/2)]))
  }
  list(par = par, value = hval, iter = iter)
}

# Representar la superficie del logaritmo de la verosimilitud
image(mmu1, mmu2, -lli, xlab = expression(mu[1]), ylab = expression(mu[2]))
contour(mmu1, mmu2, -lli, nlevels = 50, add = TRUE)
points(starts, col = "blue", pch = 19)

set.seed(1)
for (j in 1:nstarts) {
  sar <- SA(like, starts[j, ])
  with(sar, lines(par[, 1], par[, 2], lwd = 2, col = "blue"))
  with(sar, points(par[iter, 1], par[iter, 2], pch = 19))
  with(sar, cat("par = ", par[iter, ], ", value =", value[iter], "\n"))
}
```

```
 ## par =  -0.20913 2.3415 , value = 363 
 ## par =  -0.29867 2.5733 , value = 363.66 
 ## par =  -0.47085 2.426 , value = 365.33 
 ## par =  -0.34544 2.4463 , value = 363.51 
 ## par =  -0.12363 2.4648 , value = 361.74
```



\begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-21-1} \end{center}

:::


### Algoritmos genéticos

Los algoritmos genéticos tratan de encontrar la mejor solución
(entre un conjunto de soluciones posibles) imitando los procesos de
evolución biológica:

-   **Población**: formada por $n$ individuos $\mathbf{x}_i$
    codificados en **cromosomas**.

-   $f(\mathbf{x}_i)$ ajuste/capacidad/**adaptación** del
    individuo $\mathbf{x}_i$.

-   **Selección**: los individuos mejor adaptados tienen mayor
    probabilidad de ser **padres**.

-   **Cruzamiento**: los cromosomas de dos padres se combinan para
    generar hijos.

-   **Mutación**: modificación al azar del cromosoma de los
    hijos (variabilidad).

-   **Elitismo**: el mejor individuo pasa a la siguiente generación.

Los paquetes de R `DEOptim` y `gafit` implementan algunos de estos
tipos de algoritmos.


::: {.example #mv-DEoptim name="Estimación máximo-verosimil empleando un algoritmo genético"}
<br>

Repetimos el ejemplo anterior empleando el algoritmo genético implementado en la función `DEoptim::DEOptim()`:


``` r
require(DEoptim)

# Representar la superficie del logaritmo de la verosimilitud
image(mmu1, mmu2, -lli, xlab = expression(mu[1]), ylab = expression(mu[2]))
contour(mmu1, mmu2, -lli, nlevels = 50, add = TRUE)
# Estos algoritmos no requieren valores iniciales (los generan al azar en el rango)

lower <- c(-2, -2)
upper <- c(5, 5)
set.seed(1)
# DEoptim(like, lower, upper)
der <- DEoptim(tee.optim2d(like), lower, upper, DEoptim.control(itermax = 10))
```

```
 ## Iteration: 1 bestvalit: 373.132461 bestmemit:   -0.764103    2.196961
 ## Iteration: 2 bestvalit: 367.580379 bestmemit:   -0.430095    2.196961
 ## Iteration: 3 bestvalit: 367.580379 bestmemit:   -0.430095    2.196961
 ## Iteration: 4 bestvalit: 367.580379 bestmemit:   -0.430095    2.196961
 ## Iteration: 5 bestvalit: 361.906887 bestmemit:    0.058951    2.455186
 ## Iteration: 6 bestvalit: 361.906887 bestmemit:    0.058951    2.455186
 ## Iteration: 7 bestvalit: 361.906887 bestmemit:    0.058951    2.455186
 ## Iteration: 8 bestvalit: 361.657986 bestmemit:   -0.064005    2.452184
 ## Iteration: 9 bestvalit: 361.657986 bestmemit:   -0.064005    2.452184
 ## Iteration: 10 bestvalit: 361.657986 bestmemit:   -0.064005    2.452184
```

``` r
# Por defecto fija el tamaño de la población a NP = 10*npar = 20
# Puede ser mejor dejar el valor por defecto itermax = 200
points(der$optim$bestmem[1], der$optim$bestmem[2], pch = 19)
```



\begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-22-1} \end{center}

:::


## Métodos Monte Carlo en Inferencia Estadística {#aplic-inf}

***Work in progress***: Esta sección es muy preliminar y variará en siguientes ediciones.

<!-- 
Pendiente: 
  Mejorar redacción.
  Incluir leyendas de las figuras.
  Añadir distribución muestreo multidimensional, 
    media y varianza muestrales.
    regresión simple: estimaciones coeficientes y varianza.
-->


Como ya se comentó en la introducción muchas de las aplicaciones de la simulación serían de utilidad en Estadística:

-   Distribución de estimadores puntuales/estadísticos:

    -   Aproximación de la distribución.

        -   Aproximación de características de la distribución.

        -   Valided de la distribución asintótica.

    -   Comparación de estimadores.

-   Estimación por intervalo de confianza:

    -   Obtención de intervalos/bandas de confianza (probabilidad).

    -   Análisis de un estimador por intervalo de confianza.

-   Contrastes de hipótesis:

    -   Aproximación del $p$-valor.

    -   Análisis de un contraste de hipótesis.

    -   Validación teoría.

-   Métodos de remuestro bootstrap.

-   Inferencia Bayesiana

-   ...


En esta sección nos centraremos en estudios de simulación Monte Carlo en algunas de estas aplicaciones y daremos algún ejemplo de métodos Monte Carlo para inferencia estadística.
La mayoría de los métodos Monte Carlo los podríamos clasificar como métodos de remuestreo y se tratarán con mayor profundidad en capítulos siguientes.


Observación: 
En esta sección se obtendrán simulaciones de estadísticos a partir de muestras (podemos pensar que se parte de generaciones de una variable multivariante).
En la mayoría de los ejemplos se generan todas las muestras de una vez, se guardan y se procesan vectorialmente (normalmente empleando la función `apply`).
Como ya se comentó en el Capítulo \@ref(rrng), en problemas mas complejos, en los que no es necesario almacenar todas las muestras, puede ser preferible emplear un bucle para generar y procesar las muestras iterativamente. 


### Distribución en el muestreo

::: {.exercise #distr-media name="Distribución de la media muestral"}
<br> 

Si $X_{1},\ldots,X_{n}$ es una muestra aleatoria simple de una
variable aleatoria $X \sim \mathcal{N}\left( \mu, \sigma^2 \right)$, la
distribución en el muestreo de:
$$\hat{\mu}=\overline{X}=\dfrac{1}{n}\sum_{i=1}^{n}X_{i}$$
es:
$$\overline{X} \sim \mathcal{N}\left(  \mu,\dfrac{\sigma^2}{n}\right)$$ 
Confirmar este resultado mediante simulación, para ello:

a)  Crear un conjunto de datos `muestras` con 500 muestras de tamaño
    $n=10$ de una $N(1,2)$. Añadir al conjunto de datos las
    estimaciones de la media y desviación típica obtenidas con cada
    una de las muestras.
    
    Valores iniciales:
    
    
    ``` r
    set.seed(54321) # Fijar semilla para reproducibilidad
    nsim <- 500
    nx <- 10
    ```
    
    Valores teóricos:
    
    
    ``` r
    mux <- 1
    sdx <- 2
    ```
    
    Simulación de las muestras (al estilo `Rcmdr`):
    
    
    ``` r
    muestras <- as.data.frame(matrix(rnorm(nsim*nx, mean=mux, sd=sdx), ncol=nx))
    rownames(muestras) <- paste("muestra", 1:nsim, sep="")
    colnames(muestras) <- paste("obs", 1:nx, sep="")
    str(muestras)
    ```
    
    ```
     ## 'data.frame':	500 obs. of  10 variables:
     ##  $ obs1 : num  0.642 -0.856 -0.568 -2.301 0.184 ...
     ##  $ obs2 : num  3.483 2.216 1.1 4.305 0.677 ...
     ##  $ obs3 : num  1.24 -1.51 -3.98 2.29 2.46 ...
     ##  $ obs4 : num  3.286 0.947 0.953 -1.663 2.623 ...
     ##  $ obs5 : num  3.77 -1.34 1.61 -2.46 1.11 ...
     ##  $ obs6 : num  -2.044 0.32 3.046 0.136 3.555 ...
     ##  $ obs7 : num  0.6186 -1.8614 4.3386 0.0996 0.8334 ...
     ##  $ obs8 : num  -0.829 2.202 -1.688 1.534 -0.114 ...
     ##  $ obs9 : num  0.4904 -0.6713 0.5451 -0.6517 0.0168 ...
     ##  $ obs10: num  2.79 2.84 1.27 3.93 2.17 ...
    ```
    
    Estimaciones:
    
    
    ``` r
    muestras$mean <- rowMeans(muestras[,1:nx])
    muestras$sd <- apply(muestras[,1:nx], 1, sd)
    ```
    
    
    La fila `muestras[i,]` contiene las observaciones de la i-ésima muestra y 
    la correspondiente media y desviación típica.
    
    
    ``` r
    muestras[1,]
    ```
    
    ```
     ##            obs1   obs2   obs3   obs4   obs5    obs6    obs7     obs8
     ## muestra1 0.6422 3.4827 1.2425 3.2856 3.7669 -2.0444 0.61863 -0.82936
     ##             obs9  obs10   mean     sd
     ## muestra1 0.49038 2.7901 1.3445 1.9513
    ```
    
Normalmente emplearemos sin embargo una ordenación por columnas (cada fila se corresponderá con una generación).

b)  Generar el histograma (en escala de densidades) de las medias
    muestrales y compararlo con la densidad teórica.
    
    Distribución de la media muestral:
    
    
    ``` r
    hist(muestras$mean, freq = FALSE, breaks = "FD", 
         xlab = "Medias", ylab = "Densidad")
    # Densidad observada (estimación)
    lines(density(muestras$mean)) 
    # Densidad teórica (bajo normalidad)
    curve(dnorm(x, mux, sdx/sqrt(nx)), lwd = 2, col = "blue", add = TRUE) 
    # Aproximación del valor esperado de la media muestral mediante simulación
    abline(v = mean(muestras$mean), lty = 2)  
    # Valor esperado de la media muestral (teórico)
    abline(v = mux, col = "blue")
    ```
    
    \begin{figure}[!htbp]
    
    {\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/mednorm-1} 
    
    }
    
    \caption{Distribución de la media muestral de una distribución normal.}(\#fig:mednorm)
    \end{figure}

:::

<br> \vspace{0.5cm}

::: {.exercise #distr-mediab name="Distribución de la media muestral continuación"}
<br> 

Si $X_{1},\ldots,X_{n}$ es una m.a.s. de una variable aleatoria
$X$ (cualquiera) con $E\left( X \right) = \mu$ y
$Var\left( X \right) = \sigma^{2}$, por el Teorema Central del Límite, 
la distribución en el muestreo de $\hat{\mu}=\overline{X}$ se aproxima a la
normalidad:
$$\overline{X}\underset{n\rightarrow\infty}{\longrightarrow}
N\left( \mu, \dfrac{\sigma}{\sqrt{n}}\right)$$ 
Típicamente se suele considerar que esta aproximación es buena
para tamaños muestrales $n>30$,
aunque dependerá de las características de la distribución de $X$.

a)  Repetir el Ejercicio \@ref(exr:distr-media) anterior considerando muestras de una $Exp(1)$ (tener en cuenta que $X\sim Exp(\lambda)\Rightarrow\mu_{X}=\sigma_{X}=1/\lambda$).
    ¿Qué ocurre con la distribución de la media muestral?
    
    
    ``` r
    set.seed(54321) # Fijar semilla para reproducibilidad
    nsim <- 500
    nx <- 10    
    # nx <- 50
    ```
    
    Valores teóricos:
    
    
    ``` r
    lambda <- 1
    muexp <- 1/lambda
    sdexp <- muexp
    ```
    
    Simulación de las muestras:
    
    
    ``` r
    muestras2 <- as.data.frame(matrix(rexp(nsim*nx, rate=lambda), ncol=nx))
    rownames(muestras2) <- paste("muestra", 1:nsim, sep="")
    colnames(muestras2) <- paste("obs", 1:nx, sep="")
    ```
    
    Estimaciones:
    
    
    ``` r
    muestras2$mean <- rowMeans(muestras2[,1:nx])
    muestras2$sd <- apply(muestras2[,1:nx], 1, sd)
    ```
    
    Distribución de la media muestral:
    
    
    ``` r
    hist(muestras2$mean, xlim = c(-0.1, 2.5), freq = FALSE, breaks = "FD", 
         xlab = "Medias", ylab = "Densidad")
    # Densidad observada (estimación)
    lines(density(muestras2$mean)) 
    # Distribución asintótica (TCL)
    curve(dnorm(x,muexp,sdexp/sqrt(nx)), lwd=2, col="blue", add=TRUE) 
    # Aproximación del valor esperado de la media muestral mediante simulación
    abline(v=mean(muestras2$mean),lty=2)  
    # Valor esperado de la media muestral (teórico)
    abline(v=muexp, col="blue")
    ```
    
    \begin{figure}[!htbp]
    
    {\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/medexp-1} 
    
    }
    
    \caption{Distribución de la media muestral de una distribución exponencial y distribución asintótica.}(\#fig:medexp)
    \end{figure}


b)  Aumentar el tamaño muestral a 50. ¿Se aproxima más la
    distribución de las medias muestrales a la teórica bajo
    normalidad?
        
    Ejecutar el código del apartado anterior fijando `nx <- 50`.
       
:::

### Intervalos de confianza


::: {.exercise #ic-media name="Intervalo de confianza para la media"}
<br>

A partir del enunciado del Ejercicio \@ref(exr:distr-media), se deduce que el intervalo de confianza (de nivel $1-\alpha$) para la media $\mu$ de una población normal con varianza conocida es:
$$IC_{1-\alpha}\left(  \mu\right)  = 
\left(  \overline{X}-z_{1-\alpha/2}\dfrac{\sigma}{\sqrt{n}},\ \overline{X} 
+ z_{1-\alpha/2}\dfrac{\sigma}{\sqrt{n}} \right).$$
La idea es que el $100(1-\alpha)\%$ de los intervalos así construidos contentrán el verdadero valor del parámetro.

a)  Utilizando el conjunto de datos `muestras` del ejercicio 1 (500
    muestras de tamaño $n=10$ de una $N(1,2)$), añadir en dos nuevas
    variables los extremos del intervalo de confianza para la media
    con varianza conocida al conjunto de datos. Analizar la
    cobertura de estas estimaciones por IC.
        
    IC para la media con varianza conocida (bajo normalidad):
    
    
    ``` r
    alfa <- 0.05
    z <- qnorm(1 - alfa/2)
    muestras$ici <- muestras$mean - z*sdx/sqrt(nx)
    muestras$ics <- muestras$mean + z*sdx/sqrt(nx)
    ```
    
    Cobertura de las estimaciones por IC:
    
    
    ``` r
    muestras$cob <- (muestras$ici < mux) & (mux < muestras$ics) 
    ncob <- sum(muestras$cob) # Nº de intervalos que contienen la verdadera media
    ncob
    ```
    
    ```
     ## [1] 480
    ```
    
    ``` r
    100*ncob/nsim     # Proporción de intervalos
    ```
    
    ```
     ## [1] 96
    ```
    
    ``` r
    100*(1 - alfa)    # Proporción teórica bajo normalidad
    ```
    
    ```
     ## [1] 95
    ```
    
    Como ejemplo ilustrativo, generamos el gráfico de los primeros 50 intervalos:
    
    
    ``` r
    m <- 50
    tmp <- muestras[1:m,]
    attach(tmp)
    color <- ifelse(cob,"blue","red")
    plot(1:m, mean, col = color, ylim = c(min(ici),max(ics)), 
         xlab = "Muestra", ylab = "IC")
    arrows(1:m, ici, 1:m, ics, angle = 90, length = 0.05, code = 3, col = color)
    abline(h = mux, lty = 3)
    detach(tmp)
    ```
    
    \begin{figure}[!htbp]
    
    {\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/cobicnorm-1} 
    
    }
    
    \caption{Cobertura de las estimaciones por IC.}(\#fig:cobicnorm)
    \end{figure}

b)  Repetir el apartado anterior considerando muestras de una
    $Exp(1)$. ¿Qué ocurre con la cobertura del intervalo de
    confianza obtenido bajo normalidad?

    Ejecutar el código del apartado a) del ejercicio 2.
    
    IC para la media con varianza conocida (bajo normalidad)
    
    
    ``` r
    alfa <- 0.05
    z <- qnorm(1 - alfa/2)
    muestras2$ici <- muestras2$mean - z*sdexp/sqrt(nx)
    muestras2$ics <- muestras2$mean + z*sdexp/sqrt(nx)
    ```
    
    Cobertura de las estimaciones por IC:
    
    
    ``` r
    muestras2$cob <- (muestras2$ici < muexp) & (muexp < muestras2$ics) 
    ncob <- sum(muestras2$cob) # Nº de intervalos que contienen la verdadera media
    ncob
    ```
    
    ```
     ## [1] 469
    ```
    
    ``` r
    100*ncob/nsim     # Proporción de intervalos
    ```
    
    ```
     ## [1] 93.8
    ```
    
    ``` r
    100*(1 - alfa)    # Proporción teórica bajo normalidad
    ```
    
    ```
     ## [1] 95
    ```
    
    Como ejemplo ilustrativo, generamos el gráfico de los primeros 100 intervalos:
    
    
    ``` r
    m <- 100
    tmp <- muestras2[1:m,]
    attach(tmp)
    color <- ifelse(cob,"blue","red")
    plot(1:m, mean, col = color, ylim = c(min(ici),max(ics)), 
         xlab = "Muestra", ylab = "IC")
    arrows(1:m, ici, 1:m, ics, angle = 90, length = 0.05, code = 3, col = color)
    abline(h = muexp, lty = 3)
    detach(tmp)
    ```
    
    \begin{figure}[!htbp]
    
    {\centering \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/cobicexp-1} 
    
    }
    
    \caption{Cobertura de las estimaciones por IC (bajo normalidad).}(\#fig:cobicexp)
    \end{figure}

c)  ¿Qué ocurre si aumentamos el tamaño muestral a 50?
    
    Ejecutar el código del ejercicio anterior fijando `nx <- 50` 
    y el del apartado anterior.

:::

En los apartados b) y c) podíamos considerar bootstrap descrito en siguientes capítulos.

Podemos aproximar por simulación los intervalos de probabilidad de la media muestral (tendríamos una idea del valor esperado de lo que obtendríamos
con el bootstrap percentil; en este caso el estimador es insesgado...):


``` r
# Distribución de la media muestral
hist(muestras2$mean, freq=FALSE, breaks="FD", 
     main="Distribución de la media muestral", xlab="Medias", ylab="Densidad")
# Densidad observada (estimación)
lines(density(muestras2$mean), lwd=2, col='red') 
# Densidad teórica (bajo normalidad)
curve(dnorm(x,muexp,sdexp/sqrt(nx)), col="blue", add=TRUE) 
# Aproximación por simulación del valor esperado de la media muestral
abline(v=mean(muestras2$mean), lty=2)
# Valor esperado de la media muestral (teórico)
abline(v=muexp, col="blue") 
# IP bajo normalidad
ic.aprox <- apply(muestras2[ ,c('ici','ics')], 2, mean)
## ic.aprox 
##       ici       ics 
## 0.3865199 1.6261099
# Intervalo de probabilidad para la media muestral aproximado bajo normalidad
abline(v = ic.aprox, col='blue')

# Intervalo de probabilidad para la media muestral (aproximado por simulación)
ic.sim <- quantile(muestras2$mean, c(alfa/2, 1 - alfa/2))
## ic.sim
##      2.5%     97.5% 
## 0.4714233 1.8059094
# IP (aprox.) 
abline(v=ic.sim, lty=2, col='red') 
```



\begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-34-1} \end{center}


::: {.remark}
Estimaciones puntuales, por intervalo de confianza y contrastes de hipótesis
para la media con varianza desconocida bajo normalidad
se pueden obtener con la función `t.test`.
:::

<br> \vspace{0.5cm}

::: {.exercise #ic-agresti-coull name="Intervalo de confianza Agresti-Coull para una proporción"}
<br>

El Intervalo de confianza para una proporción construido usando la
aproximación normal tiene un mal comportamiento cuando el tamaño de
la muestra es pequeño. Una simple y efectiva mejora consiste en
añadir a la muestra $2a$ elementos, $a$ exitos y $a$ fracasos. Así
el intervalo de confianza al $\left(  1-\alpha\right)  100\%$ para
una proporción mejorado es:
$$\begin{aligned}
IC_{1-\alpha}^{a}\left(  p\right)   
& =\left(  \tilde{p}-z_{1-\alpha/2}\sqrt{\frac{\tilde{p}(1-\tilde{p})}{\tilde{n}}} \text{ , }
\tilde{p}+z_{1-\alpha/2}\sqrt{\frac{\tilde{p}(1-\tilde{p})}{\tilde{n}}}\right)  ,\\
\text{siendo }\tilde{n} & = n+2a \text{, } \tilde{p} = \frac{np+a}{\tilde{n}}.
\end{aligned}$$
En el caso de $a=2$ se denomina IC Agresti-Coull.

(Los apartados a) y b) están basados en los ejemplos 1.5 y 1.6 de [Suess y Trumbo](http://www.springer.com/gp/book/9780387402734), 2010)


a)  Teniendo en cuenta que la variable aleatoria $X=n\hat{p}\sim\mathcal{B}(n,p)$,
    obtener y representar gráficamente la cobertura teórica del
    intervalo de confianza estándar ($a=0$) de una proporción para
    una muestra de tamaño $n=30$, $\alpha=0.05$ y distintos valores
    de $p$ (`p.teor <- seq(1/n, 1 - 1/n, length = 1000)`).
        
    Parámetros:
    
    
    ``` r
    n <- 30
    alpha <- 0.05
    adj <- 0  # (adj <- 2 para Agresti-Coull)
    ```
    
    Probabilidades teóricas:
    
    
    ``` r
    m <- 1000
    p.teor <- seq(1/n, 1 - 1/n, length = m) 
    ```
    
    Posibles resultados:
    
    
    ``` r
    x <- 0:n
    p.est <- (x + adj)/(n + 2 * adj) 
    ic.err <- qnorm(1 - alpha/2) * sqrt(p.est * (1 - p.est)/(n + 2 * adj))  
    lcl <- p.est - ic.err 
    ucl <- p.est + ic.err 
    ```
    
    Recorrer prob. teóricas:
    
    
    ``` r
    p.cov <- numeric(m)
    for (i in 1:m) {
      # cobertura de los posibles intervalos
      cover <- (p.teor[i] >= lcl) & (p.teor[i] <= ucl)  
      # prob. de los posibles intervalos
      p.rel <- dbinom(x[cover], n, p.teor[i])           
      # prob. total de cobertura
      p.cov[i] <- sum(p.rel)                            
    }
    ```
    
    Gráfico coberturas:
    
    
    ``` r
    plot(p.teor, p.cov, type = "l", ylim = c(1 - 4 * alpha, 1))
    abline(h = 1 - alpha, lty = 2) 
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-39-1} \end{center}
    
    Fuente [Suess y Trumbo (2010)](http://www.springer.com/gp/book/9780387402734).

b)  Repetir el apartado anterior considerando intervalos de
    confianza Agresti-Coull ($a=2$).
        
    Parámetros:
    
    
    ``` r
    n <- 30
    alpha <- 0.05
    adj <- 2  # Agresti-Coull
    
    # Probabilidades teóricas:
    m <- 1000
    p.teor <- seq(1/n, 1 - 1/n, length = m) 
    # Posibles resultados:
    x <- 0:n
    p.est <- (x + adj)/(n + 2 * adj) 
    ic.err <- qnorm(1 - alpha/2) * sqrt(p.est * (1 - p.est)/(n + 2 * adj))  
    lcl <- p.est - ic.err 
    ucl <- p.est + ic.err 
    # Recorrer prob. teóricas:
    p.cov <- numeric(m)
    for (i in 1:m) {
      # cobertura de los posibles intervalos
      cover <- (p.teor[i] >= lcl) & (p.teor[i] <= ucl)  
      # prob. de los posibles intervalos
      p.rel <- dbinom(x[cover], n, p.teor[i])           
      # prob. total de cobertura
      p.cov[i] <- sum(p.rel)                            
    }
    # Gráfico coberturas:
    plot(p.teor, p.cov, type = "l", ylim = c(1 - 4 * alpha, 1))
    abline(h = 1 - alpha, lty = 2) 
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-40-1} \end{center}

c)  Repetir el apartado anterior empleando simulación para aproximar
    la cobertura.
        
    Parámetros:
    
    
    ``` r
    n <- 30
    alpha <- 0.05
    adj <- 2  #' (2 para Agresti-Coull)
    
    set.seed(54321)
    nsim <- 500 
    # Probabilidades teóricas:
    m <- 1000
    p.teor <- seq(1/n, 1 - 1/n, length = m) 
    ```
    
    Recorrer prob. teóricas:
    
    
    ``` r
    # m <- length(p.teor)
    p.cov <- numeric(m)
    for (i in 1:m) {
      # Equivalente a simular nsim muestras de tamaño n
      # ry <- matrix(rbinom(n*nsim, 1, p.teor[i]), ncol=n)
      # rx <- apply(ry, 1, sum)
      rx <- rbinom(nsim, n, p.teor[i])
      p.est <- (rx + adj)/(n + 2 * adj)  
      ic.err <- qnorm(1 - alpha/2) * sqrt(p.est * (1 - p.est)/(n + 2 * adj))
      p.cov[i] <- mean( abs(p.est - p.teor[i]) < ic.err )
    }
    ```
    
    Representar:
    
    
    ``` r
    plot(p.teor, p.cov, type = "l", ylim = c(1 - 4 * alpha, 1))
    abline(h = 1 - alpha, lty = 2) 
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-43-1} \end{center}

:::

Como ya se comentó, el caso de ajustar un modelo a los datos y realizar simulaciones a partir de ese modelo ajustado para aproximar las características de interés de un estadístico, se denomina también bootstrap paramétrico.
Para más detalles ver por ejemplo la [Sección 3.1](https://rubenfcasal.github.io/book_remuestreo/modunif-boot-par.html) de Cao y Fernández-Casal (2020).
En este libro, en las secciones [4.6.2](https://rubenfcasal.github.io/book_remuestreo/icboot-ejem.html#estudio-sim-exp) y [B.3.2](https://rubenfcasal.github.io/book_remuestreo/ejemplos-3.html#estudio-sim-boot), se incluyen ejemplos adicionales de estudios de simulación.


### Contrastes de hipótesis {#contrastes}

Ver [Capítulo 5](https://rubenfcasal.github.io/book_remuestreo/contrastes.html) de Cao y Fernández-Casal (2020).

::: {.exercise #ks-test-sim name="Test de Kolmogorov-Smirnov"}
<br>

En la Sección \@ref(calgen) del Tema \@ref(gen-pseudo) se propuso el análisis de la bondad de ajuste de un generador de números pseudo-aleatorios mediante el test de Kolmogorov-Smirnov (ver Sección \@ref(ks-test)). 
Sin embargo, si $H_{0}$ es compuesta (los parámetros desconocidos se estiman por máxima verosimilitud y se trabaja con $\hat{F}_{0}$) los cuantiles de la distribución (asintótica) de $D_{n}$ pueden ser demasiado conservativos y sería preferible utilizar la distribución exacta.

a)  Analizar el comportamiento del contraste de Kolmogorov-Smirnov
    para contrastar normalidad empleando repetidamente este test,
    considerando 1000 pruebas con muestras de tamaño 30 de 
    una $\mathcal{N}(0,1)$. Comparar gráficamente el ajuste 
    de la distribución del $p$-valor a la de referencia 
    (estudiar el tamaño del contraste).
        
    Valores iniciales:
    
    
    ``` r
    set.seed(54321)
    nx <- 30
    mx <- 0
    sx <- 1
    nsim <- 1000
    estadistico <- numeric(nsim)
    pvalor <- numeric(nsim)
    ```
    
    Realizar contrastes
    
    
    ``` r
    for(isim in 1:nsim) {
      rx <- rnorm(nx, mx, sx)
      tmp <- ks.test(rx, "pnorm", mean(rx), sd(rx))
      estadistico[isim] <- tmp$statistic
      pvalor[isim] <- tmp$p.value
    }
    ```
    
    Proporción de rechazos:
    
    
    ``` r
    {
      cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
      cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
      cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")
    }
    ```
    
    ```
     ## 
     ## Proporción de rechazos al 1% = 0 
     ## Proporción de rechazos al 5% = 0 
     ## Proporción de rechazos al 10% = 0.001
    ```
    
    Análisis de los p-valores:
    
    
    ``` r
    hist(pvalor, freq=FALSE)
    abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
    # Distribución empírica
    curve(ecdf(pvalor)(x), type = "s", lwd = 2, 
          main = 'Tamaño del contraste', ylab = 'Proporción de rechazos', 
          xlab = 'Nivel de significación')
    abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE)
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-46-1} \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-46-2} \end{center}

b)  Repetir el apartado anterior considerando el test de Lilliefors
    (rutina `lillie.test` del paquete `nortest`).

    
    
    ``` r
    library(nortest, quietly = TRUE)
    ```
    
    Valores iniciales:
    
    
    ``` r
    set.seed(54321)
    nx <- 30
    mx <- 0
    sx <- 1
    nsim <- 1000
    estadistico <- numeric(nsim)
    pvalor <- numeric(nsim)
    ```
    
    Realizar contrastes
    
    
    ``` r
    for(isim in 1:nsim) {
      rx <- rnorm(nx, mx, sx)
      # tmp <- ks.test(rx, "pnorm", mean(rx), sd(rx))
      tmp <- lillie.test(rx)
      estadistico[isim] <- tmp$statistic
      pvalor[isim] <- tmp$p.value
    }
    ```
    
    Proporción de rechazos:
    
    
    ``` r
    {
      cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
      cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
      cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")
    }
    ```
    
    ```
     ## 
     ## Proporción de rechazos al 1% = 0.01 
     ## Proporción de rechazos al 5% = 0.044 
     ## Proporción de rechazos al 10% = 0.089
    ```
    
    Análisis de los p-valores:
    
    
    ``` r
    hist(pvalor, freq=FALSE)
    abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
    
    # Distribución empírica
    curve(ecdf(pvalor)(x), type = "s", lwd = 2, main = 'Tamaño del contraste', 
          ylab = 'Proporción de rechazos', xlab = 'Nivel de significación')
    abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE)
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-51-1} \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-51-2} \end{center}

c)  Repetir el apartado a) contrastando una distribución exponencial
    y considerando 500 pruebas con muestras de tamaño 30 de una $Exp(1)$.

    Valores iniciales:
    
    
    ``` r
    set.seed(54321)
    nx <- 30
    ratex <- 1
    nsim <- 500
    estadistico <- numeric(nsim)
    pvalor <- numeric(nsim)
    ```
    
    Realizar contrastes
    
    
    ``` r
    for(isim in 1:nsim) {
      rx <- rexp(nx, ratex)
      tmp <- ks.test(rx, "pexp", 1/mean(rx))
      estadistico[isim] <- tmp$statistic
      pvalor[isim] <- tmp$p.value
    }
    ```
    
    Proporción de rechazos:
    
    
    ``` r
    {
      cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
      cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
      cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")
    }
    ```
    
    ```
     ## 
     ## Proporción de rechazos al 1% = 0 
     ## Proporción de rechazos al 5% = 0.004 
     ## Proporción de rechazos al 10% = 0.008
    ```
    
    Análisis de los p-valores:
    
    
    ``` r
    hist(pvalor, freq=FALSE)
    abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
    
    # Distribución empírica
    curve(ecdf(pvalor)(x), type = "s", lwd = 2, 
          main = 'Tamaño del contraste', ylab = 'Proporción de rechazos', 
          xlab = 'Nivel de significación')
    abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE) 
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-55-1} \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-55-2} \end{center}

d)  Diseñar una rutina que permita realizar el contraste KS de
    bondad de ajuste de una variable exponencial aproximando el
    $p$-valor por simulación y repetir el apartado anterior
    empleando esta rutina.


    
    ``` r
    ks.exp.sim <- function(x, nsim = 10^3) {
      DNAME <- deparse(substitute(x))
      METHOD <- "Kolmogorov-Smirnov Test of pexp by simulation" 
      n <- length(x)
      RATE <- 1/mean(x)
      ks.exp.stat <- function(x, rate=1/mean(x)) {
        DMinus <- pexp(sort(x), rate=rate) - (0:(n - 1))/n
        DPlus <- 1/n - DMinus
        Dn = max(c(DMinus, DPlus))
      }  
      STATISTIC <- ks.exp.stat(x, rate = RATE) 
      names(STATISTIC) <- "Dn"
      # PVAL <- 0
      # for(i in 1:nsim) {
      #   rx <- rexp(n, rate = RATE)
      #   if (STATISTIC <= ks.exp.stat(rx)) PVAL <- PVAL+1
      # }
      # PVAL <- PVAL/nsim
      # PVAL <- PVAL/(nsim + 1)
      # PVAL <- (PVAL + 1)/(nsim + 2)
      rx <- matrix(rexp(n*nsim, rate = RATE), ncol=n)
      PVAL <- mean(STATISTIC <= apply(rx, 1, ks.exp.stat))
      return(structure(list(statistic = STATISTIC, alternative = "two.sided", 
                       p.value = PVAL, method = METHOD, data.name = DNAME), 
                       class = "htest"))
    }
    ```
    
    Simulación:
    
    
    ``` r
    set.seed(54321)
    nx <- 30
    ratex <- 1
    nsim <- 500
    estadistico <- numeric(nsim)
    pvalor <- numeric(nsim)
    ```
    
    Realizar contrastes
    
    
    ``` r
    for(isim in 1:nsim) {
      rx <- rexp(nx, ratex)
      # tmp <- ks.test(rx, "pexp", 1/mean(rx))
      tmp <- ks.exp.sim(rx, nsim = 200)
      estadistico[isim] <- tmp$statistic
      pvalor[isim] <- tmp$p.value
    }
    ```
    
    Proporción de rechazos:
    
    
    ``` r
    {
      cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
      cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
      cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")
    }
    ```
    
    ```
     ## 
     ## Proporción de rechazos al 1% = 0.008 
     ## Proporción de rechazos al 5% = 0.058 
     ## Proporción de rechazos al 10% = 0.106
    ```
    
    Análisis de los p-valores:
    
    
    ``` r
    hist(pvalor, freq=FALSE)
    abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
    
    # Distribución empírica
    curve(ecdf(pvalor)(x), type = "s", lwd = 2, 
          main = 'Tamaño del contraste', ylab = 'Proporción de rechazos', 
          xlab = 'Nivel de significación')
    abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE) 
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-59-1} \includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-59-2} \end{center}

e)  Estudiar la potencia de los contrastes de los apartados c) y d),
    considerando como alternativa una distribución Weibull.

    La distribución exponencial es un caso particular de la Weibull:
    `dexp(x, ratex) == dweibull(x, 1, 1/ratex)`.
    Estudiamos lo que ocurre al desplazar  `dweibull(x, shape, 1/ratex)` con `0 < shape < 2`.
    
    CUIDADO: las simulaciones pueden requerir de mucho tiempo de computación
    (consideramos valores pequeños de `nx` y `nsim` en datos y en `ks.exp.sim`).
    
    
    ``` r
    set.seed(54321)
    nx <- 20
    ratex <- 1    # Puede ser interesante representarlo variando rate
    nsim <- 200
    alfa <- 0.1   # Puede ser interesante representarlo variando alfa
    
    
    shapex <- seq(0.25, 1.75, len=21)
    preject <- numeric(length(shapex)) # Porporciones de rechazos con ks.test
    ks.test.p <- function(x) ks.test(x, "pexp", 1/mean(x))$p.value
    preject2 <- preject # Porporciones de rechazos con ks.exp.sim
    ks.exp.sim.p <- function(x) ks.exp.sim(x, 200)$p.value
    
    for (i in seq_along(shapex)) { 
      rx <- matrix(rweibull(nx*nsim, shape = shapex[i], scale = 1/ratex), ncol=nx)
      preject[i] <- mean( apply(rx, 1, ks.test.p) <= alfa )
      preject2[i] <- mean( apply(rx, 1, ks.exp.sim.p) <= alfa )
    }
    
    plot(shapex, preject, type="l", main = paste("Potencia del contraste ( alfa =", alfa, ")"), 
         xlab = "shape", ylab = "Proporción de rechazos")
    lines(shapex, preject2, lty = 2)
    abline(h = alfa, v = 1, lty = 3)
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/potencia-1} \end{center}

:::

El estadístico de Kolmogorov-Smirnov `Dn = max(c(DMinus, DPlus))` tiene ventajas desde el
punto de vista teórico, pero puede no ser muy potente para detectar diferencias entre la
distribución bajo la hipótesis nula y la distribución de los datos. 
La ventaja de la aproximación por simulación es que no estamos atados a resultados teóricos
y podemos emplear el estadístico que se considere oportuno 
(la principal desventaja es el tiempo de computación). 
Por ejemplo, podríamos pensar en utilizar como estadístico la suma de los errores en
valor absoluto del correspondiente gráfico PP, y solo habría que cambiar el estadístico 
`Dn` en la función `ks.exp.sim` por `Dn = sum(abs( (1:n - 0.5)/n -  pexp(sort(x), rate=rate) ))`.


### Comparación de estimadores

::: {.exercise #media-median name="Comparación de la eficiencia de la media muestral y de la mediana bajo contaminación"}
<br>

Supongamos que estamos interesados en estudiar el efecto de datos atípicos en la estimación de la media teórica mediante la media y la mediana muestrales. 
Consideramos una variable aleatoria con distribución normal contaminada, en la que una observación procede de una $N(0,1)$ con probabilidad 0.95 y de una $N(3,3^2)$ con probabilidad 0.05 (mixtura). 
Se puede generar una muestra de esta variable (mixtura) mediante el método de composición descrito en la Sección \@ref(composicion), por ejemplo empleando el siguiente código:

```
p.sim <- rbinom(n, 1, 0.05)
dat.sim <- rnorm(n, 3*p.sim, 1+2*p.sim)
```

Podemos comparar la densidad objetivo con la de los valores contaminados:


``` r
curve(dnorm(x, 0, 1), -3, 12, ylab = 'densidad', lty = 3)
curve(0.95*dnorm(x, 0, 1) + 0.05*dnorm(x, 3, 3), add = TRUE)
```



\begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/contaminada-1} \end{center}
:::

::: {.remark}
Como se comentó en la Sección \@ref(composicion), también es habitual simular este tipo de datos generando un porcentaje alto de valores (en este caso un 95%) de la distribución base ($N(0,1)$) y el resto (5%) de la distibución "contaminadora" ($N(3,3^2)$), aunque se suele considerar un porcentaje de contaminación del 1% o inferior (además, como en este caso concreto no va importar el orden, no sería necesario combinar aleatoriamente los valores).
:::

<br> \vspace{0.5cm}

a)  Aproximar mediante simulación (500 generaciones) el sesgo y
    error estándar de la media y la mediana en el caso de una
    muestra de tamaño $n=100$ (suponiendo que se pretende estimar la
    media no contaminada 0).

    
    ``` r
    # media y mediana
    xsd <- 1
    xmed <- 0
    ndat <- 100
    nsim <- 500
    
    # for (isim in 1:nsim) # evitar matrix y apply
    set.seed(1)
    ntsim <- ndat*nsim
    p.sim <- rbinom(ntsim, 1, 0.05)
    dat.sim <- rnorm(ntsim, 3*p.sim, 1+2*p.sim)
    dat.sim <- matrix(dat.sim, ncol=nsim)
    ```
    
    Cada columna es una muestra
    
    
    ``` r
    str(dat.sim[,1])
    ```
    
    ```
     ##  num [1:100] 0.197 -0.42 1.163 -0.406 0.744 ...
    ```
    
    ``` r
    hist(dat.sim[,1])
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-61-1} \end{center}
    
    Calculamos los estimadores:
    
    
    ``` r
    mean.sim <- apply(dat.sim, 2, mean)
    median.sim <- apply(dat.sim, 2, median)
    ```
    
    Estimamos sus características:
    
    
    ``` r
    mean(mean.sim) # Coincide con el sesgo (media teórica es 0)
    ```
    
    ```
     ## [1] 0.146
    ```
    
    ``` r
    sd(mean.sim)
    ```
    
    ```
     ## [1] 0.13495
    ```
    
    ``` r
    mean(median.sim) # Coincide con el sesgo (media teórica es 0)
    ```
    
    ```
     ## [1] 0.044535
    ```
    
    ``` r
    sd(median.sim)
    ```
    
    ```
     ## [1] 0.13006
    ```
    
    Sesgo:
    
    
    ``` r
    boxplot(mean.sim-xmed, median.sim-xmed, 
          names=c("Media","Mediana"), ylab="Sesgo")
    abline(h = 0, lty = 2)
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-64-1} \end{center}
    
    Error cuadrático:
    
    
    ``` r
    boxplot((mean.sim-xmed)^2, (median.sim-xmed)^2, 
          names=c("Media","Mediana"), ylab="Error cuadrático")
    ```
    
    
    
    \begin{center}\includegraphics[width=0.75\linewidth]{07-Monte_Carlo_files/figure-latex/unnamed-chunk-65-1} \end{center}
    
    Estadísticos error cuadrático:
    
    
    ``` r
    # SE media
    summary((mean.sim-xmed)^2) 
    ```
    
    ```
     ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     ## 0.00000 0.00451 0.02063 0.03949 0.05915 0.36196
    ```
    
    ``` r
    # SE mediana
    summary((median.sim-xmed)^2) 
    ```
    
    ```
     ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     ## 0.00000 0.00165 0.00706 0.01886 0.02439 0.26184
    ```



