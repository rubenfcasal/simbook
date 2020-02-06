Análisis de resultados de simulación
====================================




Convergencia
------------

Supongamos que estamos interesados en aproximar la media teórica 
$E\left( X\right)$ a partir de una secuencia i.i.d. $X_{1}$,
$X_{2}$, $\cdots$, $X_{n}$ mediante la media muestral $\bar{X}_{n}$
(aproximar una probabilidad sería un caso particular considerando una distribución de Bernouilli).
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


Para detectar problemas de convergencia es recomendable representar la evolución de la aproximación de la
característica de interés (sobre el número de generaciones). 
Además del análisis descriptivo de las simulaciones. 
Por ejemplo, en el siguiente gráfico de cajas observamos los valores que producen estos saltos:


```r
boxplot(rx)
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-5-1} \end{center}


Estimación de la precisión
--------------------------

En el caso de la media muestral $\overline{X}_{n}$, un estimador
insesgado de $Var\left( \overline{X}_{n}\right) =\sigma ^{2}/n$ es:
$$\widehat{Var}\left( \overline{X}_{n}\right) =\frac{\widehat{S}^{2}}{n}$$
con:
$$\widehat{S}_{n}^{2}=\dfrac{1}{n-1}\sum\limits_{i=1}^{n}\left( X_{i}-
\overline{X}\right) ^{2}.$$

En el caso de una proporción $\hat{p}_{n}$:
$$\widehat{Var}\left( \hat{p}_{n}\right) 
=\frac{\hat{p}_{n}(1-\hat{p}_{n})}{n-1}.$$

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
\overline{X}+z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}} \right).$$

Podemos considerar que
$z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}}$ 
es la precisión obtenida (con nivel de confianza $1-\alpha$).


\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-6"><strong>(\#exm:unnamed-chunk-6) </strong></span>Aproximación de la media de una distribución normal</div>\EndKnitrBlock{example}


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
Como medida de la precisión de la aproximación podemos considerar:

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
esterr <- sqrt(cumsum((rx-est)^2))/(n - 1)  # Error estandar
plot(est, type = "l", lwd = 2, xlab = "Número de generaciones", 
     ylab = "Media y rango de error", ylim = c(-1, 1))
abline(h = est[nsim], lty=2)
lines(est + 2*esterr, lty=3)
lines(est - 2*esterr, lty=3)
abline(h = xmed)
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-10-1} \end{center}

 
Determinación del número de generaciones
----------------------------------------

En muchas ocasiones puede interesar obtener una aproximación con un
nivel de precisión fijado.

Para una precisión absoluta $\varepsilon$, se trata de determinar
$n$ de forma que:
$$z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}}<\varepsilon$$

Un algoritmo podría ser el siguiente:

1.  Hacer $j=0$
    y fijar un tamaño inicial $n_{0}$ (e.g. 30 ó 60).

2.  Generar $\left\{ X_{i}\right\} _{i=1}^{n_{0}}$ y calcular $\widehat{S}_{n_{0}}$.

3.  Mientras $\left. z_{1-\alpha /2}\widehat{S}_{n_{j-1}}\right/ 
    \sqrt{n_{j}}>\varepsilon$ hacer:

4.  $\qquad j=j+1$.

5.  $\qquad n_{j}=\left\lceil \left( \left. z_{1-\alpha /2}\widehat{S}
      _{n_{j-1}}\right/ \varepsilon \right)^{2}\right\rceil$.

6.  $\qquad$ Generar $\left\{ X_{i}\right\}
      _{i=n_{j-1}+1}^{n_{j}}$ y calcular $\widehat{S}_{n_{j}}$.



Para una precisión relativa $\varepsilon \left\vert \mu
\right\vert$ se procede análogamente de forma
que:$$z_{1-\alpha /2}\dfrac{\widehat{S}_{n}}{\sqrt{n}}<\varepsilon \left\vert 
\overline{X}_{n}\right\vert .$$


El problema de la dependencia
-----------------------------

En el caso de dependencia, la estimación de la precisión se complica:
$$Var\left( \overline{X}\right) =\frac{1}{n^{2}}\left( 
\sum_{i=1}^{n}Var\left( X_{i} \right) + 2\sum_{i<j}Cov\left( X_{i},X_{j}\right) \right).$$

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-11"><strong>(\#exm:unnamed-chunk-11) </strong></span>Aproximación de una proporción bajo dependencia (cadena de Markov)</div>\EndKnitrBlock{example}
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



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-13-1} \end{center}

La aproximación de la proporción sería correcta (es consistente):

```r
est[nsim]
```

```
## [1] 0.3038
```

Sin embargo, al ser datos dependientes **esta aproximación del error estandar no es adecuada**:

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



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-16-1} \end{center}

El gráfico de autocorrelaciones sugiere que si tomamos 1 de cada 25 
podemos suponer independencia.



```r
lag <- 24
xlag <- c(rep(FALSE, lag), TRUE)
rxi <- rx[xlag]
acf(as.numeric(rxi))
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-17-1} \end{center}

```r
n <- 1:length(rxi)
est <- cumsum(rxi)/n
esterr <- sqrt(est*(1-est)/(n-1))
plot(est, type="l", lwd=2, ylab="Probabilidad", 
     xlab=paste("Número de simulaciones /", lag + 1), ylim=c(0,0.6))
abline(h = est[length(rxi)], lty=2)
lines(est + 2*esterr, lty=2) # Supone independencia
lines(est - 2*esterr, lty=2)
abline(h = 1/3, col="darkgray")     # Prob. teor. cadenas Markov
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-17-2} \end{center}


Esta forma de proceder podría ser adecuada para tratar de aproximar la precisión 
pero no sería eficiente para aproximar la media. Siempre será preferible emplear
todas las observaciones. 

Por ejemplo, se podría pensar en considerar las medias de grupos de 24 valores 
consecutivos y suponer que hay independencia entre ellas:


```r
rxm <- rowMeans(matrix(rx, ncol = lag, byrow = TRUE))
n <- 1:length(rxm)
est <- cumsum(rxm)/n
esterr <- sqrt(cumsum((rxm-est)^2))/(n - 1)  # Error estándar
plot(est, type="l", lwd=2, ylab="Probabilidad", 
     xlab=paste("Número de simulaciones /", lag + 1), ylim=c(0,0.6))
abline(h = est[length(rxm)], lty=2)
lines(est + 2*esterr, lty=2) # OJO! Supone independencia
lines(est - 2*esterr, lty=2)
abline(h = 1/3, col="darkgray")     # Prob. teor. cadenas Markov
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-18-1} \end{center}

Esta es la idea del método de medias por lotes 
(*batch means*; *macro-micro replicaciones*) para la estimación de la varianza.
Alternativamente se podría recurrir a la generación de múltiples secuencias
independientes.

Trataremos este tipo de problemas en la diagnosis de algoritmos de
simulación Monte Carlo de Cadenas de Markov (MCMC). 
Aparecen también en la simulación dinámica (por eventos o cuantos). 

**Nota**: En el caso anterior se calcula el error estándar de la aproximación por simulación de la proporción,
pero si el objetivo es la aproximación de la varianza (de la variable y no de las medias por lotes), 
habrá que reescalarlo adecuadamente:


```r
var.aprox <- nsim * esterr[length(rxm)]^2
var.aprox
```

```
## [1] 2.473664
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


### Periodo de calentamiento

En el caso de simulación de datos dependientes (simulación dinámica) 
pueden aparecer problemas de estabilización. Puede ocurrir que el sistema 
evolucione lentamente en el tiempo hasta alcanzar su distribución estacionaria,
siendo muy sensible a las condiciones iniciales con las que se comienzó la
simulación. En tal caso resulta conveniente ignorar los resultados obtenidos
durante un cierto período inicial de tiempo (denominado período de calentamiento 
o estabilización), cuyo único objeto es conseguir que se estabilice la distribución de
probabilidad.

Como ejemplo comparamos la simulación anterior con la obtenida 
considerando como punto de partida un día lluvioso.


```r
set.seed(1)
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
abline(v = 2000, lty = 3)
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-21-1} \end{center}


En estos casos puede ser recomendable ignorar los primeros valores generados (por ejemplo los primeros 2000) y recalcular los 
estadísticos deseados.

También trataremos este tipo de problemas en la diagnosis de algoritmos MCMC. 

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-22"><strong>(\#exm:unnamed-chunk-22) </strong></span>Simulación de un proceso autorregresivo (serie de tiempo)</div>\EndKnitrBlock{example}

$$X_t = \mu + \rho * (X_{t-1} - \mu) + \varepsilon_t$$
Podemos tener en cuenta que en este caso la varianza es:
$$\textrm{var}(X_t)=\operatorname{E}(X_t^2)-\mu^2=\frac{\sigma_\varepsilon^2}{1-\rho^2}.$$

Establecer parámetros

```r
nsim <- 200   # Numero de simulaciones
xvar <- 1     # Varianza
xmed <- 0     # Media
rho <- 0.5    # Coeficiente AR
nburn <- 10   # Periodo de calentamiento (burn-in)
# Varianza del error
evar <- xvar*(1 - rho^2)
```
Alternativamente se podría fijar la varianza del error:

```r
#   evar <- 1     
#   xvar <- evar / (1 - rho^2)
```
Simular


```r
set.seed(1)
x <- numeric(nsim + nburn)
# Normalmente simularíamos el primer valor
#   rx[1] <- rnorm(1, mean = xmed, sd = sqrt(xvar))
# Aunque en este caso lo fijamos para alejarnos un poco de la distribución estacionaria
# (lo que requerirá de un mayor periodo de calentamiento)
x[1] <- -10
# Simular el resto de la secuencia
for (i in 2:length(x))
  x[i] <- xmed + rho*(x[i-1] - xmed) + rnorm(1, sd=sqrt(evar))
x <- as.ts(x)
plot(x)
abline(v = nburn, lty = 2)
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-Analisis_resultados_files/figure-latex/unnamed-chunk-25-1} \end{center}

```r
# Eliminar periodo de calentamiento
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

-   Incluso en el caso de la media, las bandas de confianza
    obtenidas con el TCL son puntuales (si generamos nuevas
    secuencias de simulación es muy probable que no
    estén contenidas).

-   En muchos casos (p.e. la generación de múltiples secuencias de
    simulación puede suponer un coste computacional importante),
    puede ser preferible emplear un método de remuestreo.

