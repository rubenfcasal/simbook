Análisis de resultados de simulación {#cap4}
====================================




Work in progress...

En este capítulo nos centraremos en la aproximación mediante simulación de la media teórica de un estadístico a partir de la media muestral de una secuencia de simulaciones de dicho estadístico.
La aproximación de una probabilidad sería un caso particular considerando una variable de Bernouilli.

En primer lugar se tratará el análisis de la convergencia y la precisión de la aproximación por simulación. 
Al final del capítulo se incluye una breve introducción a los problemas de estabilización y dependencia (con los que nos solemos encontrar en simulación dinámica y MCMC).


Convergencia {#convergencia}
------------

Supongamos que estamos interesados en aproximar la media teórica 
$\mu = E\left( X\right)$ a partir de una secuencia i.i.d. $X_{1}$,
$X_{2}$, $\cdots$, $X_{n}$ mediante la media muestral $\bar{X}_{n}$.
Una justificación teórica de la validez de la aproximación obtenida
mediante simulación es *la ley (débil) de los grandes números*:

-   Si $X_{1}$, $X_{2}$, $\cdots$ es una secuencia de v.a.'s
    independientes con:
    $$E\left( X_{i}\right) =\mu \text{ y }Var\left( X_{i}\right) 
    =\sigma^{2}<\infty,$$
    entonces $\overline{X}_{n}=\left( X_{1}+\cdots +X_{n}\right) /n$ 
    converge en probabilidad a $\mu$. i.e. para cualquier $\varepsilon >0$:
    $$\lim\limits_{n\rightarrow \infty }P\left( \left\vert \overline{X}_{n}-\mu
    \right\vert <\varepsilon \right) = 1.$$

-   *La ley fuerte* establece la convergencia casi segura.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-2"><strong>(\#exm:unnamed-chunk-2) </strong></span>Aproximación de una probabilidad</div>\EndKnitrBlock{example}
Simulamos una distribución de Bernoulli de parámetro $p=0.5$:

```r
p <- 0.5
set.seed(1)
nsim <- 10000
# nsim <- 100
rx <- runif(nsim) <= p
```
La aproximación por simulación de $p$ será:

```r
mean(rx) 
```

```
## [1] 0.5047
```
Podemos generar un gráfico con la evolución de la aproximación con el siguiente código:

```r
plot(cumsum(rx)/1:nsim, type="l", lwd=2, xlab="Número de generaciones", 
     ylab="Proporción muestral", ylim=c(0,1))
abline(h = mean(rx), lty = 2)
# valor teórico
abline(h = p) 
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/proporcion-1} 

}

\caption{Aproximación de la proporción en función del número de generaciones.}(\#fig:proporcion)
\end{figure}

### Detección de problemas de convergencia

Una suposición crucial es que las variables $X_{i}$ deben tener
varianza finita (realmente esta suposición puede relajarse:
$E\left( \left\vert X_{i} \right\vert \right) < \infty$).
En caso contrario la media muestral puede no converger a
una constante. Un ejemplo conocido es la distribución de Cauchy:


```r
set.seed(1)
nsim <- 10000
rx <- rcauchy(nsim)
plot(cumsum(rx)/1:nsim, type="l", lwd=2, 
     xlab="Número de generaciones", ylab="Media muestral")
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/cauchy-1} 

}

\caption{Evolución de la media muestral de una distribución de Cauchy en función del número de generaciones.}(\#fig:cauchy)
\end{figure}


Para detectar problemas de convergencia es recomendable representar la evolución de la aproximación de la característica de interés (sobre el número de generaciones), 
además de realizar otros análisis descriptivos de las simulaciones. 
Por ejemplo, en este caso podemos observar los valores que producen estos saltos mediante un gráfico de cajas: 


```r
boxplot(rx)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/cauchy-box-1} 

}

\caption{Gráfico de cajas de 10000 generaciones de una distribución de Cauchy.}(\#fig:cauchy-box)
\end{figure}


Estimación de la precisión
--------------------------

En el caso de la media muestral $\overline{X}_{n}$, un estimador
insesgado de $Var\left( \overline{X}_{n}\right) =\sigma ^{2}/n$ es la cuasi-varianza muestral:
$$\widehat{Var}\left( \overline{X}_{n}\right) =\frac{\widehat{S}^{2}}{n}$$
con:
$$\widehat{S}_{n}^{2}=\dfrac{1}{n-1}\sum\limits_{i=1}^{n}\left( X_{i}-
\overline{X}\right) ^{2}.$$

En el caso de una proporción $\hat{p}_{n}$:
$$\widehat{Var}\left( \hat{p}_{n}\right) 
=\frac{\hat{p}_{n}(1-\hat{p}_{n})}{n-1},$$
aunque se suele emplear la varianza muestral.

Los valores obtenidos servirían como medidas básicas de la precisión
de la aproximación, aunque su principal aplicación es la
construcción de intervalos de confianza.

Teorema central del límite
--------------------------

Si $X_{1}$, $X_{2}$, $\cdots$ es una secuencia de v.a.'s
independientes con $E\left( X_{i}\right) =\mu$ y
$Var\left( X_{i}\right) = \sigma ^{2}<\infty$, entonces:
$$Z_{n}=\frac{\overline{X}_{n}-\mu }{\frac{\sigma }{\sqrt{n}}}
\overset{d}{ \rightarrow } N(0,1)$$
i.e. $\lim\limits_{n\rightarrow \infty }F_{Z_{n}}(z)=\Phi (z)$.
Por tanto, un intervalo de confianza asintótico para $\mu$ es:
$$IC_{1-\alpha }(\mu ) = \left( \overline{X}_{n}
- z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}},\ 
\overline{X}_n+z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}} \right).$$

Podemos considerar que
$z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}}$ 
es la precisión obtenida (con nivel de confianza $1-\alpha$).

La convergencia de la aproximación, además de ser aleatoria, se podría considerar lenta.
La idea es que para doblar la precisión (disminuir el error a la mitad), necesitaríamos un número de generaciones cuatro veces mayor. Pero una ventaja, es que este error no depende del número de dimensiones (en el caso multidimensional puede ser mucho más rápida que otras alternativas numéricas).

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-5"><strong>(\#exm:unnamed-chunk-5) </strong></span>Aproximación de la media de una distribución normal</div>\EndKnitrBlock{example}


```r
xsd <- 1
xmed <- 0
set.seed(1)
nsim <- 1000
rx <- rnorm(nsim, xmed, xsd)
```
La aproximación por simulación de la media será:

```r
mean(rx)
```

```
## [1] -0.01164814
```
Como medida de la precisión de la aproximación podemos considerar (se suele denominar error de la aproximación):

```r
2*sd(rx)/sqrt(nsim)
```

```
## [1] 0.06545382
```
(es habitual emplear 2 en lugar de 1.96, 
lo que se correspondería con $1 - \alpha = 0.9545$ en el caso de normalidad).
Podemos añadir también los correspondientes intervalos de confianza al gráfico de convergencia:

```r
n <- 1:nsim
est <- cumsum(rx)/n
# Errores estándar calculados con la varianza muestral por comodidad:
esterr <- sqrt(cumsum((rx-est)^2))/n  
plot(est, type = "l", lwd = 2, xlab = "Número de generaciones", 
     ylab = "Media y rango de error", ylim = c(-1, 1))
abline(h = est[nsim], lty=2)
lines(est + 2*esterr, lty=3)
lines(est - 2*esterr, lty=3)
abline(h = xmed)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/conv-esterr-1} 

}

\caption{Gráfico de convergencia incluyendo el error de la aproximación.}(\#fig:conv-esterr)
\end{figure}

 
Determinación del número de generaciones
----------------------------------------

Normalmente el valor de $n$ se toma del orden de varias centenas o millares. 
En los casos en los que la simulación se utiliza para aproximar una característica central de la distribución (como una media) puede bastar un número de simulaciones del orden de $n = 100, 200, 500$. 
Sin embargo, en otros casos pueden ser necesarios valores del tipo $B = 1000, 2000, 5000, 10000$.

En muchas ocasiones puede interesar obtener una aproximación con un nivel de precisión fijado.
Para una precisión absoluta $\varepsilon$, se trata de determinar
$n$ de forma que:
$$z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}}<\varepsilon$$

Un algoritmo podría ser el siguiente:

1.  Hacer $j=0$
    y fijar un tamaño inicial $n_{0}$ (e.g. 30 ó 60).

2.  Generar $\left\{ X_{i}\right\} _{i=1}^{n_{0}}$ 
    y calcular $\overline{X}_{n_0}$ y $\widehat{S}_{n_{0}}$.

3.  Mientras $\left. z_{1-\alpha /2}\widehat{S}_{n_j}\right/ \sqrt{n_{j}}>\varepsilon$ hacer:

    3.1. $j=j+1$.
    
    3.2. $n_{j}=\left\lceil \left( \left. z_{1-\alpha /2}\widehat{S}
         _{n_{j-1}}\right/ \varepsilon \right)^{2}\right\rceil$.
    
    3.3. Generar $\left\{ X_{i}\right\}_{i=n_{j-1}+1}^{n_j}$ 
         y calcular $\overline{X}_{n_j}$ y $\widehat{S}_{n_j}$.

3.  Devolver $\overline{X}_{n_j}$ y $\left. z_{1-\alpha /2}\widehat{S}_{n_j}\right/ \sqrt{n_{j}}$.

Para una precisión relativa $\varepsilon \left\vert \mu \right\vert$ se procede análogamente de forma que:
$$z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}}<\varepsilon \left\vert 
\overline{X}_{n}\right\vert .$$


El problema de la dependencia
-----------------------------

En el caso de dependencia, la estimación de la precisión se complica:
$$Var\left( \overline{X}\right) =\frac{1}{n^{2}}\left( 
\sum_{i=1}^{n}Var\left( X_{i} \right) + 2\sum_{i<j}Cov\left( X_{i},X_{j}\right) \right).$$

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:mmc"><strong>(\#exm:mmc) </strong></span>Aproximación de una proporción bajo dependencia (cadena de Markov)</div>\EndKnitrBlock{example}
Supongamos que en A Coruña llueve de media 1/3 días al año,
y que la probabilidad de que un día llueva solo depende de lo que ocurrió el día anterior, 
siendo 0.94 si el día anterior llovió y 0.03 si no.
Podemos generar valores de la variable indicadora de día lluvioso con el siguiente código:

```r
# Variable dicotómica 0/1 (FALSE/TRUE)  
set.seed(1)
nsim <- 10000
alpha <- 0.03 # prob de cambio si seco
beta <- 0.06  # prob de cambio si lluvia
rx <- logical(nsim) # x == "llueve"
rx[1] <- FALSE # El primer día no llueve
for (i in 2:nsim)
  rx[i] <- if (rx[i-1]) runif(1) > beta else runif(1) < alpha
```
Se podría pensar en emplear las expresiones anteriores:

```r
n <- 1:nsim
est <- cumsum(rx)/n
esterr <- sqrt(est*(1-est)/(n-1)) # OJO! Supone independencia
plot(est, type="l", lwd=2, ylab="Probabilidad", 
     xlab="Número de simulaciones", ylim=c(0,0.6))
abline(h = est[nsim], lty=2)
lines(est + 2*esterr, lty=2) 
lines(est - 2*esterr, lty=2)
abline(h = 1/3, col="darkgray")     # Prob. teor. cadenas Markov
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/conv-dep-1} 

}

\caption{Gráfico de convergencia incluyendo el error de la aproximación (calculado asumiendo independencia).}(\#fig:conv-dep)
\end{figure}

La aproximación de la proporción sería correcta (es consistente):

```r
est[nsim]
```

```
## [1] 0.3038
```

Sin embargo, al ser datos dependientes esta aproximación del error estandar no es adecuada:

```r
esterr[nsim]
```

```
## [1] 0.004599203
```

En este caso al haber dependencia positiva se produce una 
subestimación del verdadero error estandar.


```r
acf(as.numeric(rx))
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/acf-depsec-1} 

}

\caption{Correlograma de la secuencia indicadora de días de lluvia.}(\#fig:acf-depsec)
\end{figure}

El gráfico de autocorrelaciones sugiere que si tomamos 1 de cada 25 
podemos suponer independencia.


```r
lag <- 24
xlag <- c(rep(FALSE, lag), TRUE)
rxi <- rx[xlag]
acf(as.numeric(rxi))
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/acf-depsec2-1} 

}

\caption{Correlograma de la subsecuencia de días de lluvia obtenida seleccionando uno de cada 25 valores.}(\#fig:acf-depsec2)
\end{figure}


```r
nrxi <- length(rxi)
n <- 1:nrxi
est <- cumsum(rxi)/n
esterr <- sqrt(est*(1-est)/(n-1))
plot(est, type="l", lwd=2, ylab="Probabilidad", 
     xlab=paste("Número de simulaciones /", lag + 1), ylim=c(0,0.6))
abline(h = est[length(rxi)], lty=2)
lines(est + 2*esterr, lty=2) # Supone independencia
lines(est - 2*esterr, lty=2)
abline(h = 1/3, col="darkgray")     # Prob. teor. cadenas Markov
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/conv-dep2-1} 

}

\caption{Gráfico de convergencia de la aproximación de la probabilidad a partir de la subsecuencia de días de lluvia (calculando el error de aproximación asumiendo independencia).}(\#fig:conv-dep2)
\end{figure}


Esta forma de proceder podría ser adecuada para tratar de aproximar la precisión:

```r
esterr[nrxi]
```

```
## [1] 0.02277402
```
pero no sería eficiente para aproximar la media. Siempre será preferible emplear
todas las observaciones. 

Por ejemplo, se podría pensar en considerar las medias de grupos de 24 valores 
consecutivos y suponer que hay independencia entre ellas:


```r
rxm <- rowMeans(matrix(rx, ncol = lag, byrow = TRUE))
nrxm <- length(rxm)
n <- 1:nrxm
est <- cumsum(rxm)/n
esterr <- sqrt(cumsum((rxm-est)^2))/n  # Error estándar
plot(est, type="l", lwd=2, ylab="Probabilidad", 
     xlab=paste("Número de simulaciones /", lag + 1), ylim=c(0,0.6))
abline(h = est[length(rxm)], lty=2)
lines(est + 2*esterr, lty=2) # OJO! Supone independencia
lines(est - 2*esterr, lty=2)
abline(h = 1/3, col="darkgray")     # Prob. teor. cadenas Markov
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/conv-dep-lotes-1} 

}

\caption{Gráfico de convergencia de las medias por lotes.}(\#fig:conv-dep-lotes)
\end{figure}

Esta es la idea del método de medias por lotes 
(*batch means*; *macro-micro replicaciones*) para la estimación de la varianza.
En el ejemplo anterior se calcula el error estándar de la aproximación por simulación de la proporción:

```r
esterr[nrxm]
```

```
## [1] 0.01569017
```
pero si el objetivo es la aproximación de la varianza (de la variable y no de las medias por lotes), habrá que reescalarlo adecuadamente. 
Supongamos que la correlación entre $X_i$ y $X_{i+k}$ es aproximadamente nula,
y consideramos las subsecuencias (lotes) $(X_{t+1},X_{t+2},\ldots,X_{t+k})$ con $t=(j-1)k$, $j=1,\ldots,m$ y $n = mk$. 
Entonces:

$$\begin{aligned}
Var \left(\bar X \right) &= Var \left(\frac{1}{n} \sum_{i=1}^n X_i\right) 
= Var \left( \frac{1}{m}\sum_{j=1}^m \left(\frac{1}{k} \sum_{t=(i-1)k + 1}^{ik} X_t\right) \right) \\
&\approx \frac{1}{m^2} \sum_{j=1}^m Var \left(\frac{1}{k} \sum_{t=(i-1)k + 1}^{ik} X_t\right)
\approx \frac{1}{m} Var \left(\bar{X}_k \right)
\end{aligned}$$
donde $\bar{X}_k$ es la media de una subsecuencia de longitud $k$.

<!-- Ecuaciones p.237 Gentle, Random numbers & MC methods -->


```r
var.aprox <- nsim * esterr[length(rxm)]^2
var.aprox
```

```
## [1] 2.461814
```

Obtenida asumiendo independencia entre las medias por lotes, y que será
una mejor aproximación que asumir independencia entre las generaciones
de la variable:


```r
var(rx)
```

```
## [1] 0.2115267
```

Alternativamente se podría recurrir a la generación de múltiples secuencias
independientes entre sí: 

```r
# Variable dicotómica 0/1 (FALSE/TRUE)  
set.seed(1)
nsim <- 1000
nsec <- 10
alpha <- 0.03 # prob de cambio si seco
beta <- 0.06  # prob de cambio si lluvia
rxm <- matrix(FALSE, nrow = nsec, ncol= nsim)
for (i in 1:nsec) {
  # rxm[i, 1] <- FALSE # El primer día no llueve
  # rxm[i, 1] <- runif(1) < 1/2 # El primer día llueve con probabilidad 1/2
  rxm[i, 1] <- runif(1) < 1/3 # El primer día llueve con probabilidad 1/3 (ideal)
  for (j in 2:nsim)
    rxm[i, j] <- if (rxm[i, j-1]) runif(1) > beta else runif(1) < alpha
}
```
La idea sería considerar las medias de las series como una muestra independiente
de una nueva variable y estimar su varianza de la forma habitual:

```r
# Media de cada secuencia
n <- 1:nsim
est <- apply(rxm, 1, function(x) cumsum(x)/n)
matplot(n, est, type = 'l', lty = 3, col = "lightgray",
     ylab="Probabilidad", xlab="Número de simulaciones")
# Aproximación
mest <- apply(est, 1, mean)
lines(mest, lwd = 2)
abline(h = mest[nsim], lty = 2)
# Precisión
mesterr <- apply(est, 1, sd)/sqrt(nsec)
lines(mest + 2*mesterr, lty = 2)
lines(mest - 2*mesterr, lty = 2)
# Prob. teor. cadenas Markov
abline(h = 1/3, col="darkgray")     
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/conv-dep-nsec-1} 

}

\caption{Gráfico de convergencia de la media de 10 secuencias generadas de forma independiente.}(\#fig:conv-dep-nsec)
\end{figure}

```r
# Aproximación final
mest[nsim] # mean(rxm)
```

```
## [1] 0.3089
```

```r
# Error estándar
mesterr[nsim]
```

```
## [1] 0.02403491
```
Trataremos este tipo de problemas en la diagnosis de algoritmos de
simulación Monte Carlo de Cadenas de Markov (MCMC). 
Aparecen también en la simulación dinámica (por eventos o cuantos). 


### Periodo de calentamiento

En el caso de simulación de datos dependientes (simulación dinámica) 
pueden aparecer problemas de estabilización. Puede ocurrir que el sistema 
evolucione lentamente en el tiempo hasta alcanzar su distribución estacionaria,
siendo muy sensible a las condiciones iniciales con las que se comienzó la
simulación. En tal caso resulta conveniente ignorar los resultados obtenidos
durante un cierto período inicial de tiempo (denominado período de calentamiento 
o estabilización), cuyo único objeto es conseguir que se estabilice la distribución de
probabilidad.

Como ejemplo comparamos la simulación del Ejemplo \@ref(exm:mmc) con la obtenida considerando como punto de partida un día lluvioso (con una semilla distinta para evitar dependencia).


```r
set.seed(2)
nsim <- 10000
rx2 <- logical(nsim)
rx2[1] <- TRUE # El primer día llueve
for (i in 2:nsim)
  rx2[i] <- if (rx2[i-1]) runif(1) > beta else runif(1) < alpha
n <- 1:nsim
est <- cumsum(rx)/n
est2 <- cumsum(rx2)/n
plot(est, type="l", ylab="Probabilidad", 
     xlab="Número de simulaciones", ylim=c(0,0.6))
lines(est2, lty = 2)
# Ejemplo periodo calentamiento nburn = 2000
abline(v = 2000, lty = 3)
# Prob. teor. cadenas Markov
abline(h = 1/3, col="darkgray")     
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-17-1} \end{center}


En estos casos puede ser recomendable ignorar los primeros valores generados (por ejemplo los primeros 2000) y recalcular los 
estadísticos deseados.

También trataremos este tipo de problemas en la diagnosis de algoritmos MCMC. 

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-18"><strong>(\#exm:unnamed-chunk-18) </strong></span>Simulación de un proceso autorregresivo (serie de tiempo)</div>\EndKnitrBlock{example}

$$X_t = \mu + \rho * (X_{t-1} - \mu) + \varepsilon_t$$
Podemos tener en cuenta que en este caso la varianza es:
$$\textrm{var}(X_t)=\operatorname{E}(X_t^2)-\mu^2=\frac{\sigma_\varepsilon^2}{1-\rho^2}.$$

Establecemos los parámetros:

```r
nsim <- 200   # Numero de simulaciones
xmed <- 0     # Media
rho <- 0.5    # Coeficiente AR
nburn <- 10   # Periodo de calentamiento (burn-in)
```
Se podría fijar la varianza del error:

```r
evar <- 1
# Varianza de la respuesta
xvar <- evar / (1 - rho^2)
```
pero la recomendación sería fijar la varianza de la respuesta:

```r
xvar <- 1     
# Varianza del error
evar <- xvar*(1 - rho^2)
```

Para simular la serie, al ser un $AR(1)$, normalmente simularíamos el primer valor

```r
rx[1] <- rnorm(1, mean = xmed, sd = sqrt(xvar))
```
o lo fijamos a la media (en este caso nos alejamos un poco de la distribución estacionaria, para que el "periodo de calentamiento" sea mayor). 
Después generamos los siguientes valores de forma recursiva:

```r
set.seed(1)
x <- numeric(nsim + nburn)
# Establecer el primer valor 
x[1] <- -10
# Simular el resto de la secuencia
for (i in 2:length(x))
  x[i] <- xmed + rho*(x[i-1] - xmed) + rnorm(1, sd=sqrt(evar))
x <- as.ts(x)
plot(x)
abline(v = nburn, lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/sim-ar1-1} 

}

\caption{Ejemplo de una simulación de una serie de tiempo autorregresiva.}(\#fig:sim-ar1)
\end{figure}
y eliminamos el periodo de calentamiento:

```r
rx <- x[-seq_len(nburn)]
```



Para simular una serie de tiempo en `R`
se puede emplear la función `arima.sim` del paquete base `stats`.
En este caso el periodo de calentamiento se establece mediante el
parámetro `n.start` (que se fija automáticamente a un valor adecuado).

Por ejemplo, podemos generar este serie autoregressiva con:

```r
rx2 <- arima.sim(list(order = c(1,0,0), ar = rho), n = nsim, n.start = nburn, sd = sqrt(evar))
```
La recomendación es fijar la varianza de las series simuladas si se quieren comparar
resultados considerando distintos parámetros de dependencia.



Observaciones
-------------

-   En el caso de que la característica de interés de la
    distribución de $X$ no sea la media, los resultados anteriores
    no serían en principio aplicables.

-   Incluso en el caso de la media, las "bandas de confianza"
    obtenidas con el TCL son puntuales (si generamos nuevas
    secuencias de simulación es muy probable que no
    estén contenidas).

-   En muchos casos (p.e. la generación de múltiples secuencias de
    simulación puede suponer un coste computacional importante),
    puede ser preferible emplear un método de remuestreo.

