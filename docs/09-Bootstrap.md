# Métodos de remuestreo {#bootstrap}

<!-- Capítulo \@ref(bootstrap) -->



<!-- 
---
title: "Métodos de remuestreo"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    pandoc_args: ["--number-offset", "8,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("09-Bootstrap.Rmd")
knitr::purl("09-Bootstrap.Rmd", documentation = 2)
knitr::spin("09-Bootstrap.R",knit = FALSE)

Pendiente: 
Selección del estadístico
  Interesaría un estadístico pivotal
  Bootstrap percentil: invariante frente a transformaciones monótonas, aunque debería ser insesgada y con varianza independiente del parámetro
  Bootstrap básico o natural: varianza independiente del parámetro
  Boostrap estudentizado
  
Esto da pie a una de las consideraciones más importantes a la hora
de diseñar un buen método de remuestreo bootstrap: ha de procurarse que el bootstrap imite todas las condiciones que cumple la población original.  

Aplicaciones bootstrap
Ejercicios
  mediana
  coeficiente de correlación
  
Incluir expresiones teóricas en ejemplo en introducción
Incluir aquí el caso multidimensional?

-->

Un par de notas:

- Etimología: se denomina *bootstrap* a la cinta de la bota (oreja lateral o trasera para ayudar a calzarse las botas). 
 
- Modismo anglosajón: "to pull oneself up by one’s bootstraps", que podríamos traducir como resolver un problema con medios propios, sin la ayuda de otros.
    Se cree que esta frase está basada en un libro del siglo XVIII:
    
    > "The Baron had fallen to the bottom of a deep lake. 
    >  Just when it looked like all was lost, 
    >  he thought to pick himself up by his own bootstraps".
    >
    > --- Raspe, Rudolph Erich (The Surprising Adventures of Baron Munchausen, 1785)


El bootstrap es un procedimiento estadístico que sirve para aproximar características de la distribución en el muestreo de un estadístico.
Para ello se emplea (normalmente) simulación, generando un gran número de muestras mediante algún tipo de remuestreo de la muestra original.

Su ventaja principal es que no requiere hipótesis sobre el mecanismo generador de los datos (aunque los resultados asintóticos requieren de hipótesis generales). 
Por lo que son de especial utilidad cuando no se dispone la distribución exacta del estadístico y no es
posible o adecuado emplear la distribución asintótica.

En este capítulo se incluye una breve introducción al bootstrap desde un punto de vista aplicado.
Para información adicional, inluyendo resultados teóricos, ver por ejemplo Davison y Hinkley (1997) o [Cao y Fernández-Casal (2020)](https://rubenfcasal.github.io/book_remuestreo).

En este libro nos centraremos principalmente en los métodos de remuestreo bootstrap, aunque hay otros tipos de remuestreo como el *jackknife*, para la aproximación del sesgo y varianza de un estimador (ver [Sección 2.2](https://rubenfcasal.github.io/book_remuestreo/jackknife.html) de Cao y Fernández-Casal, 2020), o los empleados en contrastes de permutaciones (ver [Sección 5.3](https://rubenfcasal.github.io/book_remuestreo/contrastes-de-permutaciones.html) de Cao y Fernández-Casal, 2020).

En los siguientes capítulos se tratarán algunas de las principales aplicaciones de los métodos bootstrap.
Entre ellas podríamos destacar:

-   Aproximación del sesgo y de la varianza de un estimador.

-   Construcción de intervalos de confianza.

-   Contraste de hipótesis.

También la simulación (condicional) de nuevas observaciones o la estimación de la probabilidad de superar un determinado umbral (probabilidad de riesgo).

## Introducción {#boot-intro}

En primer lugar nos centraremos en la idea original del bootstrap uniforme (Efron, 1979; también denominado bootstrap no paramétrico o naïve), que se expondrá de manera más formal en la Sección \@ref(boot-unif).

Supongamos que $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$ es una una muestra aleatoria simple (m.a.s.) 
de una población con distribución $F_{\theta}$ y que estamos interesados en hacer inferencia sobre $\theta$ empleando un estimador $\hat{\theta} = T\left( \mathbf{X} \right)$.

La idea es aproximar características poblacionales por las correspondientes de la distribución empírica de los
datos observados. 
Se trata de imitar el experimento de muestreo en la población real, pero empleando la distribución empírica en lugar de la distribución teórica $F_{\theta}$ desconocida.
Al conocer el mecanismo que genera los datos en el universo bootstrap, se puede emplear Monte Carlo para simularlo.
En el caso i.i.d. esto puede ser implementado mediante remuestreo, realizando repetidamente **muestreo aleatorio con reemplazamiento del conjunto de datos original** (manteniendo el tamaño muestral).

Para aproximar la distribución en el muestreo por Monte Carlo, se genera un número grande $B$ de réplicas bootstrap:

-   $\mathbf{X}^{\ast (b)}=\left( X_1^{\ast (b)},\ldots ,X_n^{\ast (b)} \right)$ muestra bootstrap (remuestra),
    obtenida mediante muestreo con reemplazamiento de $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$.

-   $\hat{\theta}^{\ast (b)} = T\left( \mathbf{X}^{\ast (b)} \right)$
    valor del estadístico en la muestra bootstrap (réplica bootstrap del estadístico).

para $b = 1,\ldots ,B$.

La idea original (bootstrap natural, Efron) es que la variabilidad de $\hat{\theta}_{b}^{\ast }$ en torno a $\hat{\theta}$ aproxima la variabilidad de $\hat{\theta}$ en torno a $\theta$: **la distribución de** $\hat{\theta}_{b}^{\ast }-\hat{\theta}$ (en el universo bootstrap) **aproxima la distribución de ** $\hat{\theta}-\theta$ (en la población).

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{images/bootstrap} 

}

\caption{Esquema de la idea del boostrap (de Efron).}(\#fig:unnamed-chunk-1)
\end{figure}

En general podríamos decir que **la muestra es a la población** **lo que la muestra bootstrap es a la muestra**.

::: {.example #mean-boot name="aproximación bootstrap de la distribución de la media muestral"}
<br>

Como ejemplo ilustrativo consideramos una muestra simulada de tamaño $n=100$ de una normal estándar y la media muestral como estimador de la media teórica:

```r
set.seed(1)
n <- 100
mean_teor <- 0
sd_teor <- 1  
muestra <- rnorm(n, mean = mean_teor, sd = sd_teor)
```

El valor del estadístico en la muestra es: 

```r
estadistico <- mean(muestra)
```

Representamos la distribución de la muestra [Figura \@ref(fig:muestra-sim)]:

```r
hist(muestra, freq = FALSE, xlim = c(-3, 3),
     main = '', xlab = 'x', ylab = 'densidad')
abline(v = estadistico, lty = 2)
curve(dnorm, col = "blue", add = TRUE)
abline(v = mean_teor, col = "blue", lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{09-Bootstrap_files/figure-latex/muestra-sim-1} 

}

\caption{Distribución de la muestra simulada (y distribución teórica en azul).}(\#fig:muestra-sim)
\end{figure}

Como aproximación de la distribución poblacional, desconocida en la práctica, siempre podemos considerar la distribución empírica (o una versión suavizada: bootstrap suavizado; Sección \@ref(modunif-boot-suav)). 
Alternativamente podríamos asumir un modelo paramétrico y estimar los parámetros a partir de la muestra (bootstrap paramétrico; Sección \@ref(modunif-boot-par) [Figura \@ref(fig:muestra-sim-aprox)].

```r
# Distribución bootstrap uniforme
curve(ecdf(muestra)(x), xlim = c(-3, 3), ylab = "F(x)", type = "s")
# Distribución bootstrap paramétrico (asumiendo normalidad)
curve(pnorm(x, mean(muestra), sd(muestra)), lty = 2, add = TRUE)
# Distribución teórica
curve(pnorm, col = "blue", add = TRUE)
legend("bottomright", legend = c("Empírica", "Aprox. paramétrica", "Teórica"),
       lty = c(1, 2, 1), col = c("black","black", "blue"))
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{09-Bootstrap_files/figure-latex/muestra-sim-aprox-1} 

}

\caption{Distribución teórica de la muestra simulada y distintas aproximaciones.}(\#fig:muestra-sim-aprox)
\end{figure}

En este caso (bootstrap uniforme) generamos las réplicas bootstrap empleando la distribución empírica:

```r
set.seed(1)
B <- 1000
estadistico_boot <- numeric(B)
for (k in 1:B) {
    remuestra <- sample(muestra, n, replace = TRUE)
    estadistico_boot[k] <- mean(remuestra)
}
```

Podríamos emplear directamente las réplicas bootstrap del estimador para aproximar la distribución en el muestreo de la media muestral (esto es lo que se conoce como bootstrap percentil directo, o simplemente bootstrap percentil):


```r
hist(estadistico_boot, freq = FALSE, xlim = c(-0.2, 0.5),
     ylab = "Densidad", main = "")
# Valor esperado bootstrap del estadístico
mean_boot <- mean(estadistico_boot)  
abline(v = mean_boot, lwd = 2)
# abline(v = estadistico, col = "blue")

# Distribución poblacional
curve(dnorm(x, mean_teor, sd_teor/sqrt(n)), col = "blue", add = TRUE)
abline(v = 0, col = "blue", lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{09-Bootstrap_files/figure-latex/mean-boot-perc-1} 

}

\caption{Aproximación de la distribución de la media muestral centrada mediante bootstrap percentil (uniforme).}(\#fig:mean-boot-perc)
\end{figure}

Sin embargo, especialmente si el estimador es sesgado, puede ser preferible emplear la distribución de $\hat{\theta}_{b}^{\ast }-\hat{\theta}$ como aproximación de la distribución de $\hat{\theta}-\theta$ (bootstrap natural, básico o percentil básico):


```r
hist(estadistico_boot - estadistico, freq = FALSE, 
     ylab = "Densidad", main = "")
# Distribución teórica
curve(dnorm(x, 0, sd_teor/sqrt(n)), col = "blue", add = TRUE)
abline(v = 0, col = "blue", lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{09-Bootstrap_files/figure-latex/mean-boot-basico-1} 

}

\caption{Aproximación de la distribución de la media muestral mediante bootstrap natural (uniforme).}(\#fig:mean-boot-basico)
\end{figure}

Sin embargo, asintóticamente ambos procedimientos son equivalentes^[Por este motivo en algunas referencias más teóricas no se diferencia entre ambos métodos y se denominan simplemente bootstrap percentil.] y pueden dar lugar a los mismos resultados en determinados problemas de inferencia.
Por ejemplo en la aproximación del sesgo y de la varianza de un estimador (Sección \@ref(sesgo-prec)):


```r
# Sesgo (teor=0)
mean_boot - estadistico # mean(estadistico_boot - estadistico)
```

```
 ## [1] 0.004714973
```

```r
# Error estándar
sd(estadistico_boot) # sd(estadistico_boot - estadistico)
```

```
 ## [1] 0.08610306
```

```r
# Error estándar teórico
sd_teor/sqrt(n) 
```

```
 ## [1] 0.1
```

:::


## El Bootstrap uniforme {#boot-unif}

Suponemos que $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$ es una una muestra aleatoria simple (m.a.s.) 
de una población con distribución $F$ y que estamos interesados en hacer inferencia sobre $\theta =\theta \left(F \right)$ empleando un estimador $\hat{\theta} = T\left( \mathbf{X} \right)$.
Para ello nos gustaría conocer la distribución en el muestreo de un estadístico $R\left( \mathbf{X},F \right)$, función del estimador (y por tanto de la muestra) y de la distribución poblacional.
Por ejemplo el estimador studentizado: 
$$R=R\left( \mathbf{X},F \right) = \frac{\hat \theta - \theta}{\sqrt{\widehat{Var}(\hat \theta)}}.$$

A veces podemos calcular directamente la distribución de $R\left( \mathbf{X},F \right)$, aunque suele depender de cantidades poblacionales, no conocidas en la práctica. 
Otras veces sólo podemos llegar a aproximar la distribución de $R\left( \mathbf{X},F \right)$ cuando $n \rightarrow \infty$.
Por ejemplo, bajo normalidad $X_i \overset{i.i.d.}{\sim} \mathcal{N}\left( \mu ,\sigma^2 \right)$, si el parámetro de interés es la media 
$$\theta \left( F \right) =\mu =\int x~dF\left( x \right) =\int xf\left( x \right) ~dx$$
y consideramos como estimador la media muestral
$$\hat{\theta} = T\left( \mathbf{X} \right) = \theta \left( F_n \right) =\int x~dF_n\left( x \right) =\bar{X}.$$
Como normalmente en la práctica la varianza es desconocida, podríamos considerar el estadístico:
$$R=R\left( \mathbf{X},F \right) =\sqrt{n}\frac{\bar{X}-\mu }{S_{n-1}} \sim t_{n-1},$$
donde $S_{n-1}^2$ es la cuasivarianza muestral:
$$S_{n-1}^2=\frac{1}{n-1}\sum_{j=1}^{n}\left( X_j-\bar{X} \right)^2.$$
Si $F$ no es normal entonces la distribución de $R$ ya no es una $t_{n-1}$, 
aunque, bajo condiciones muy generales, $R\overset{d}{\rightarrow}\mathcal{N}\left(0,1 \right)$.

En el universo bootstrap se reemplaza la distribución poblacional (desconocida) $F$ por una estimación, $\hat{F}$, de la misma. 
A partir de la aproximación $\hat{F}$ podríamos generar, condicionalmente a la muestra observada, remuestras 
$$\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots ,X_n^{\ast} \right)$$ 
con distribución $X_i^{\ast} \sim \hat{F}$, que demoninaremos remuestras bootstrap.
Por lo que podemos hablar de la distribución en el remuestreo de 
$$R^{\ast}=R\left( \mathbf{X}^{\ast},\hat{F} \right),$$ 
denominada *distribución bootstrap*.

**Una de las consideraciones más importantes** al diseñar un buen método de remuestreo bootstrap es **imitar por completo el procedimiento de muestreo en la población original** (incluyendo el estadístico y las características de su distribución muestral).  

Como ya se comentó anteriormente, en el bootstrap uniforme se emplea como aproximación la distribución empírica (ver Sección \@ref(empdistr)):
$$F_n\left( x \right) =\frac{1}{n}\sum_{i=1}^{n}\mathbf{1}\left\{ X_i\leq x\right\}.$$
Es decir, $\hat{F}=F_n$, y por tanto $R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right)$. 
En raras ocasiones (e.g. [Sección 1.3](https://rubenfcasal.github.io/book_remuestreo/c%C3%A1lculo-de-la-distribuci%C3%B3n-bootstrap-exacta-y-aproximada.html) de Cao y Fernández-Casal, 2020) es posible calcular exactamente la distribución bootstrap de $R^{\ast}$. 
Normalmente se aproxima esa distribución mediante Monte Carlo, generando una gran cantidad, $B$, de réplicas de $R^{\ast}$. 
En el caso del bootstrap uniforme, el algoritmo es:

1. Para cada $i=1,\ldots ,n$ generar $X_i^{\ast}$ a partir de $F_n$, es decir
$P^{\ast}\left( X_i^{\ast}=X_j \right) =\frac{1}{n}$, $j=1,\ldots,n$

2. Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
,X_n^{\ast} \right)$

3. Calcular $R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right)$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$

5. Utilizar esas réplicas bootstrap para aproximar la distribución en el
muestreo de $R$.

Para la elección del número de réplicas Monte Carlo $B$ se aplicarían las mismas recomendaciones de la Sección \@ref(num-gen) para el caso general de una aproximación por simulación.

Como ya se mostró anteriormente, el paso 1 se puede llevar a cabo simulando una distribución uniforme discreta mediante el método de la transformación cuantil (Algoritmo \@ref(cnj:unif-discr)):

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$

Aunque en `R` es recomendable^[De esta forma se evitan posibles problemas numéricos al emplear el método de la transformación cuantil cuando $n$ es extremadamente grande (e.g. <https://stat.ethz.ch/pipermail/r-devel/2018-September/076817.html>).] emplear la función `sample` para generar muestras aleatorias con reemplazamiento del conjunto de datos original:

```r
muestra_boot <- sample(muestra, replace = TRUE)
```

::: {.example #media-dt-desconocida name="Inferencia sobre la media con varianza desconocida"}
<br> 

Como ejemplo consideramos la muestra de tiempos de vida de microorganismos `lifetimes` del paquete `simres` [Figura \@ref(fig:microorganismos)]: 

```r
library(simres)
muestra <- lifetimes
summary(muestra)
```

```
 ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 ##  0.1430  0.2650  0.6110  0.8053  1.1200  2.0800
```

```r
hist(muestra)
rug(muestra)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{09-Bootstrap_files/figure-latex/microorganismos-1} 

}

\caption{Distribución del tiempo de vida de una muestra de microorganismos.}(\#fig:microorganismos)
\end{figure}

Supongamos que queremos obtener una estimación por intervalo de confianza de su vida media a partir de los 15 valores observados mediante bootstrap uniforme considerando el estadístico
$$R=R\left( \mathbf{X},F \right) =\sqrt{n}\frac{\bar{X}-\mu }{S_{n-1}},$$
(en la Sección \@ref(boot-ic) se tratará con más detalle la construcción de intervalos de confianza).

En el bootstrap uniforme se emplea $\hat{F}=F_n\,$, con lo cual el análogo bootstrap del estadístico $R$ será
$$R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right) =\sqrt{n}\frac{\bar{X}^{\ast}-\bar{X}}{S_{n-1}^{\ast}},$$
siendo
$$\begin{aligned}
\bar{X}^{\ast} &= \frac{1}{n}\sum_{i=1}^{n}X_i^{\ast}, \\
S_{n-1}^{\ast 2} &= \frac{1}{n-1}\sum_{i=1}^{n}\left( X_i^{\ast}-
\bar{X}^{\ast} \right)^2.
\end{aligned}$$

El algoritmo bootstrap (aproximado por Monte Carlo) procedería así:

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$

2. Obtener $\bar{X}^{\ast}$ y $S_{n-1}^{\ast 2}$

3. Calcular
$R^{\ast}=\sqrt{n}\frac{\bar{X}^{\ast}-\bar{X}}{
S_{n-1}^{\ast}}$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$

5. Aproximar la distribución en el muestreo de $R$ mediante la
distribución empírica de $R^{\ast (1)}, \ldots, R^{\ast (B)}$

Por ejemplo, podríamos emplear el siguiente código:

```r
n <- length(muestra)
alfa <- 0.05
# Estimaciones muestrales
x_barra <- mean(muestra)
cuasi_dt <- sd(muestra)
# Remuestreo
set.seed(1)
B <- 1000
estadistico_boot <- numeric(B)
for (k in 1:B) {
  remuestra <- sample(muestra, n, replace = TRUE)
  x_barra_boot <- mean(remuestra)
  cuasi_dt_boot <- sd(remuestra)
  estadistico_boot[k] <- sqrt(n) * (x_barra_boot - x_barra)/cuasi_dt_boot
}
```

Las características de interés de la distribución en el muestreo de $R$ 
se aproximan por las correspondientes de la distribución bootstrap de $R^{\ast}$.
En este caso nos interesa aproximar los puntos críticos $x_{\alpha /2}$ y
$x_{1-\alpha /2}$, tales que:
$$P\left( x_{\alpha /2} < R < x_{1-\alpha /2} \right) = 1-\alpha.$$
Para lo que podemos emplear los cuantiles muestrales^[
Se podrían considerar distintos estimadores del cuantil $x_{\alpha}$ 
(ver p.e. la ayuda de la función `quantile()`).
Si empleamos directamente la distribución empírica, el cuantil se 
correspondería con la observación ordenada en la posición $B \alpha$ 
(se suele hacer una interpolación lineal si este valor no es entero), 
lo que equivale a emplear la función `quantile()` de `R` con el parámetro 
`type = 1`. Esta función considera por defecto la posición 
$1 + (B - 1) \alpha$ (`type = 7`).
En el libro de Davison y Hinkley (1997), y en el paquete `boot`, se emplea $(B + 1) \alpha$ (equivalente a `type = 6`; lo que justifica que
consideren habitualmente 99, 199 ó 999 réplicas bootstrap).]:

```r
pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2))
pto_crit
```

```
 ##      2.5%     97.5% 
 ## -3.002197  1.877302
```

A partir de los cuales obtenemos la correspondiente estimación por IC boostrap:
$$\hat{IC}^{boot}_{1-\alpha}\left(  \mu\right)  = 
\left(  \overline{X}-x_{1-\alpha/2}\dfrac{S_{n-1}}{\sqrt{n}},\ \overline{X} 
- x_{\alpha/2}\dfrac{S_{n-1}}{\sqrt{n}} \right).$$

```r
ic_inf_boot <- x_barra - pto_crit[2] * cuasi_dt/sqrt(n)
ic_sup_boot <- x_barra - pto_crit[1] * cuasi_dt/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
 ##      2.5%     97.5% 
 ## 0.5030131 1.2888063
```

Este procedimiento para la construcción de intervalos de confianza se denomina *método percentil-t* y se tratará en la Sección \@ref(boot-ic-stud).

Como ejemplo adicional podemos comparar la aproximación de la distribución bootstrap del estadístico con la aproximación $t_{n-1}$ basada en normalidad.


```r
hist(estadistico_boot, freq = FALSE, ylim = c(0, 0.4))
abline(v = pto_crit, lty = 2)
curve(dt(x, n - 1), add = TRUE, col = "blue")
pto_crit_t <- qt(1 - alfa/2, n - 1)
abline(v = c(-pto_crit_t, pto_crit_t), col = "blue", lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{09-Bootstrap_files/figure-latex/mean-boot-stud-1} 

}

\caption{Aproximación de la distribución de la media muestral studentizada mediante bootstrap uniforme.}(\#fig:mean-boot-stud)
\end{figure}

En este caso la distribución bootstrap del estadístico es más asimétrica, por lo que el intervalo de confianza no está centrado en la media,
al contrario que el obtenido con la aproximación tradicional.
Por ejemplo, podemos obtener la estimación basada en normalidad mediante la función `t.test()`:

```r
t.test(muestra)$conf.int
```

```
 ## [1] 0.4599374 1.1507292
 ## attr(,"conf.level")
 ## [1] 0.95
```

:::

---

En el caso multidimensional, cuando trabajamos con un conjunto de datos con múltiples variables, podríamos emplear un procedimiento análogo, a partir de remuestras del vector de índices. 
Por ejemplo:

```r
data(iris)
str(iris)
```

```
 ## 'data.frame':	150 obs. of  5 variables:
 ##  $ Sepal.Length: num  5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
 ##  $ Sepal.Width : num  3.5 3 3.2 3.1 3.6 3.9 3.4 3.4 2.9 3.1 ...
 ##  $ Petal.Length: num  1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
 ##  $ Petal.Width : num  0.2 0.2 0.2 0.2 0.2 0.4 0.3 0.2 0.2 0.1 ...
 ##  $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 ..
```

```r
n <- nrow(iris)
# i_boot <- floor(n*runif(n)) + 1
# i_boot <- sample.int(n, replace = TRUE)
i_boot <- sample(n, replace = TRUE)
data_boot <- iris[i_boot, ]
str(data_boot)
```

```
 ## 'data.frame':	150 obs. of  5 variables:
 ##  $ Sepal.Length: num  5.1 5.6 6.2 4.8 5.5 6.2 5.5 5.6 5 6.5 ...
 ##  $ Sepal.Width : num  3.8 2.5 2.9 3.1 2.3 2.9 2.6 2.8 3.6 3 ...
 ##  $ Petal.Length: num  1.9 3.9 4.3 1.6 4 4.3 4.4 4.9 1.4 5.8 ...
 ##  $ Petal.Width : num  0.4 1.1 1.3 0.2 1.3 1.3 1.2 2 0.2 2.2 ...
 ##  $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 2 2 1 2 2 2 3 ..
```

Esta forma de proceder es la que emplea por defecto el paquete `boot` que describiremos más adelante (Sección \@ref(intro-pkgboot)).

::: {.exercise #unif-multi name="Bootstrap uniforme multidimensional"}
<br>
Considerando el conjunto de datos `Prestige` del paquete `carData`, supongamos que queremos realizar inferencias sobre el coeficiente de correlación entre `prestige` (puntuación de ocupaciones obtenidas a partir de una encuesta) e `income` (media de ingresos en la ocupación). 
Para ello podemos considerar el coeficiente de correlación lineal de Pearson:
$$\rho =\frac{ Cov \left( X, Y \right) }
{ \sigma \left( X \right) \sigma \left( Y \right) }$$
Su estimador es el coeficiente de correlación muestral:
$$r=\frac{\sum_{i=1}^{n}(x_i-\overline{x})(y_i-\overline{y})}
{\sqrt{ \sum_{i=1}^{n}(x_i-\overline{x})^{2}} 
\sqrt{\sum_{i=1}^{n}(y_i-\overline{y})^{2}}},$$
que podemos calcular en `R` empleando la función `cor()`:

```r
data(Prestige, package = "carData")
# with(Prestige, cor(income, prestige))
cor(Prestige$income, Prestige$prestige)
```

```
 ## [1] 0.7149057
```
Para realizar inferencias sobre el coeficiente de correlación, como aproximación más simple, se puede considerar que la distribución muestral de $r$ es aproximadamente normal de media $\rho$ y varianza
$$Var(r) \approx \frac{1 - \rho^2}{n - 2}.$$
<!-- De donde se deduciría emplear el estadístico: 
$$R = \sqrt{n - 2} \frac{r -\rho}{\sqrt{1 - r^2}}$$ -->

Aproximar mediante bootstrap uniforme (multididimensional) la distribución del estadístico $R = r -\rho$, empleando $B=1000$ réplicas, y compararla con la aproximación normal, considerando
$$\widehat{Var}(r) = \frac{1 - r^2}{n - 2}.$$



:::


## Herramientas disponibles en R sobre bootstrap  {#intro-paquetes}

En `R` hay una gran cantidad de paquetes que implementan métodos bootstrap.
Por ejemplo, al ejecutar el comando `??bootstrap` (o `help.search('bootstrap')`)
se mostrarán las funciones de los paquetes instalados que incluyen este término
en su documentación (se puede realizar la búsqueda en todos los paquetes disponibles
de `R` a través de <https://www.rdocumentation.org>).

De entre todos estas herramientas destacan dos librerías 
como las más empleadas:

* `bootstrap`: contiene las rutinas (bootstrap, cross-validation,
  jackknife) y los datos del libro "An Introduction to the Bootstrap" de B.
  Efron y R. Tibshirani, 1993, Chapman and Hall. La librería fue
  desarrollada originalmente en `S` por Rob Tibshirani y exportada a `R` por
  Friedrich Leisch. Es útil para desarrollar los ejemplos que se citan en
  ese libro.

* `boot`: incluye las funciones y conjuntos de datos utilizados en el libro 
  "Bootstrap Methods and Their Applications" de A. C. Davison y D. V. Hinkley, 1997,
  Cambridge University Press. Esta librería fue desarrollada originalmente 
  en `S` por Angelo J. Canty y posteriormente exportada a `R` (ver [Canty, 2002](http://cran.fhcrc.org/doc/Rnews/Rnews_2002-3.pdf)).
  Este paquete es mucho más completo que el paquete `bootstrap`, forma parte de la distribución estándar de `R` y es el que emplearemos como referencia en este libro (ver Sección \@ref(intro-pkgboot)).

Por otra parte existen numerosas rutinas (scripts) realizadas en `R` por
diversos autores, que están disponibles en Internet 
(por ejemplo, puede ser interesante realizar una búsqueda en 
<https://rseek.org>). 

El bootstrap uniforme se puede implementar fácilmente. Por ejemplo,
una rutina general para el caso univariante sería la siguiente:


```r
#' @param x vector que contiene la muestra.
#' @param B número de réplicas bootstrap.
#' @param statistic función que calcula el estadístico.
boot.strap0 <- function(x, B=1000, statistic=mean){
  ndat <- length(x)
  x.boot <- sample(x, ndat*B, replace=TRUE)
  x.boot <- matrix(x.boot, ncol=B, nrow=ndat)
  stat.boot <- apply(x.boot, 2, statistic)
}
```

Podríamos aplicar esta función a la muestra de tiempos de vida de
microorganismos con el siguiente código:

```r
fstatistic0 <- function(x){
  mean(x)
}

B <- 1000
set.seed(1)
stat.dat <- fstatistic0(muestra)
stat.boot <- boot.strap0(muestra, B, fstatistic0)

res.boot <- c(stat.dat, mean(stat.boot)-stat.dat, sd(stat.boot))
names(res.boot) <- c("Estadístico", "Sesgo", "Error Std.")
res.boot
```

```
 ## Estadístico       Sesgo  Error Std. 
 ## 0.805333333 0.003173267 0.154099013
```

La función `boot.strap0()` anterior no es adecuada para el caso multivariante
(por ejemplo cuando estamos interesados en regresión).
Como se mostró en la Sección \@ref(boot-intro)
sería preferible emplear remuestras del vector de índices. Por ejemplo:


```r
#' @param datos vector, matriz o data.frame que contiene los datos.
#' @param B número de réplicas bootstrap.
#' @param statistic función con al menos dos parámetros, 
#' los datos y el vector de índices de remuestreo, 
#' y que devuelve el vector de estadísticos.
#' @param ... parámetros adicionales de la función statistic.
boot.strap <- function(datos, B=1000, statistic, ...) {
  ndat <- NROW(datos)
  i.boot <- sample(ndat, ndat*B, replace=TRUE)
  i.boot <- matrix(i.boot, ncol=B, nrow=ndat)
  stat.boot <- drop(apply(i.boot, 2, function(i) statistic(datos, i, ...)))
}
```

El paquete `boot`, descrito a continuación, emplea una implementación similar.

### El paquete `boot` {#intro-pkgboot}

La función principal de este paquete es la función `boot()` que implementa distintos métodos de remuestreo para datos i.i.d..
En su forma más simple permite realizar bootstrap uniforme (que en la práctica también se denomina habitualmente *bootstrap noparamétrico*):

```r
boot(data, statistic, R)
```
donde `data` es un vector, matriz o `data.frame` que contiene los datos, 
`R` es el número de réplicas bootstrap, y `statistic` es una función 
con al menos dos parámetros (con las opciones por defecto), 
los datos y el vector de índices de remuestreo, 
y que devuelve el vector de estadísticos.

Por ejemplo, para hacer inferencia sobre la mediana del tiempo de vida de microorganismos,
podríamos emplear el siguiente código:

```r
library(boot)
muestra <- lifetimes

statistic <- function(data, i){
  # remuestra <- data[i]; mean(remuestra)
  mean(data[i])
}

set.seed(1)
res.boot <- boot(muestra, statistic, R = 1000)
```

El resultado que devuelve esta función es un objeto de clase `boot`, una lista con los siguientes componentes:

```r
names(res.boot)
```

```
 ##  [1] "t0"        "t"         "R"         "data"      "seed"      "statistic"
 ##  [7] "sim"       "call"      "stype"     "strata"    "weights"
```
Además de los parámetros de entrada (incluyendo los valores por defecto), contiene tres componentes adicionales:

* `tO`: el valor observado del estadístico 
  (su evaluación en los datos originales).
  
* `t`: la matriz de réplicas bootstrap del estadístico
  (cada fila se corresponde con una remuestra).
  
* `seed`: el valor inicial de la semilla (`.Random.seed`)
  empleada para la generación de las réplicas.

Este tipo de objetos dispone de dos métodos principales:
el método `print()` que muestra un resumen de los resultados
(incluyendo  aproximaciones bootstrap del sesgo y del error
estándar de los estadísticos; ver Sección \@ref(sesgo-prec)):

```r
res.boot
```

```
 ## 
 ## ORDINARY NONPARAMETRIC BOOTSTRAP
 ## 
 ## 
 ## Call:
 ## boot(data = muestra, statistic = statistic, R = 1000)
 ## 
 ## 
 ## Bootstrap Statistics :
 ##      original      bias    std. error
 ## t1* 0.8053333 0.003173267   0.1583306
```
y el método `plot()` que genera gráficas básicas de diagnosis
de los resultados (correspondientes al estadístico determinado por el parámetro `index`, por defecto `= 1`):  [Figura \@ref(fig:plot-res-boot)]


```r
plot(res.boot)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.9\linewidth]{09-Bootstrap_files/figure-latex/plot-res-boot-1} 

}

\caption{Gráficos de diagnóstico de los resultados bootstrap para la media muestral de los tiempos de vida de microorganismos.}(\#fig:plot-res-boot)
\end{figure}

Es recomendable examinar la distribución bootstrap del estimador (o estadístico) para detectar posibles problemas.
Como en este caso puede ocurrir que el estadístico bootstrap tome pocos valores distintos, lo que indicaría que el número de réplicas bootstrap es insuficiente o que hay algún problema con método de remuestreo empleado. 
Se darán más detalles sobre los posibles problemas del bootstrap uniforme en la Sección \@ref(deficien-unif).

::: {.exercise #boot-simul-mediana}
<br>
Repetir el ejemplo anterior considerando simultáneamente la media truncada al 10\% y la mediana como estimadores de la posición central de los datos. Estudiar si hay algún problema con su distribución bootstrap (NOTA: al igual que en el caso anterior, las distribuciones objetivo serían continuas, asumiendo que la distribución del tiempo de vida es continua).
:::

Además de estos métodos, las principales funciones de interés serían:

* `jack.after.boot()`: genera un gráfico para diagnósticar la inluencia 
  de las observaciones individuales en los resultados bootstrap 
  (se representan los cuantiles frente a las diferencias en el estadístico 
  al eliminar una observación; este gráfico también se puede obtener estableciendo
  `jack = TRUE` en `plot.boot()`).
  
* `boot.array()`: genera la matriz de índices a partir de la que se obtuvieron las remuestras (permite reconstruir las remuestras bootstrap).
  
* `boot.ci()`: construye distintos tipos de intervalos de confianza 
  (se tratarán en el Sección \@ref(boot-ic)) dependiendo del parámetro `type`:
  
    - `"norm"`: utiliza la distribución asintótica normal considerando las
        aproximaciones bootstrap del sesgo y de la varianza.
        
    - `"basic"`: emplea el estadístico $R = \hat \theta - \theta$ para la
      construcción del intervalo de confianza.
      
    - `"stud"`: calcula el intervalo a partir del estadístico studentizado 
        $R = \left( \hat \theta - \theta \right) / \sqrt{\widehat{Var}(\hat \theta)}$.
    
    - `"perc"`: utiliza directamente la distribución bootstrap del estadístico
      ($R = \hat \theta$).
    
    - `"bca"`: emplea el método $BCa$ ("bias-corrected and accelerated") 
        propuesto por Efron (1987) (ver Sección 5.3.2 de Davison y Hinkley, 1997).
    
    - `"all"`: calcula los cinco tipos de intervalos anteriores.


Como ya se comentó, la función `boot()` admite estadísticos multivariantes 
(haciendo que la función `statistic` devuelva un vector en lugar de un escalar),
pero por defecto las funciones anteriores consideran el primer componente
como el estadístico principal. 
Para obtener resultados de otros componentes del vector de estadísticos
habrá que establecer el parámetro `index` igual al índice deseado.
Además, en algunos casos (por ejemplo para la obtención de intevalos de confianza
estudentizados con la función `boot.ci()`) se supone, por defecto, que el segundo
componente del vector de estadísticos contiene estimaciones de la varianza del
estadístico para cada réplica boostrap.

::: {.example #media-dt-desconocida-boot name="Inferencia sobre la media con varianza desconocida continuación"}
<br> 

Continuando con el Ejemplo \@ref(exm:media-dt-desconocida) de
inferencia sobre la media con varianza desconocida. 
Para obtener la estimación por intervalo de confianza del tiempo de vida medio 
de los microorganismos con el paquete `boot`, podríamos emplear
el siguiente código:


```r
library(boot)
muestra <- lifetimes

statistic <- function(data, i){
  remuestra <- data[i]
  c(mean(remuestra), var(remuestra)/length(remuestra))
}

set.seed(1)
res.boot <- boot(muestra, statistic, R = 1000)
res.boot
```

```
 ## 
 ## ORDINARY NONPARAMETRIC BOOTSTRAP
 ## 
 ## 
 ## Call:
 ## boot(data = muestra, statistic = statistic, R = 1000)
 ## 
 ## 
 ## Bootstrap Statistics :
 ##      original       bias    std. error
 ## t1* 0.8053333  0.003173267 0.158330646
 ## t2* 0.0259338 -0.002155755 0.007594682
```

```r
boot.ci(res.boot)
```

```
 ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
 ## Based on 1000 bootstrap replicates
 ## 
 ## CALL : 
 ## boot.ci(boot.out = res.boot)
 ## 
 ## Intervals : 
 ## Level      Normal              Basic             Studentized     
 ## 95%   ( 0.4918,  1.1125 )   ( 0.4825,  1.0980 )   ( 0.4715,  1.2320 )  
 ## 
 ## Level     Percentile            BCa          
 ## 95%   ( 0.5127,  1.1282 )   ( 0.5384,  1.1543 )  
 ## Calculations and Intervals on Original Scale
```

El intervalo marcado como `Studentized` se obtuvo empleando el mismo estadístico del Ejemplo \@ref(exm:media-dt-desconocida).

:::

### Extensiones del bootstrap uniforme con el paquete `boot` 

Estableciendo parámetros adicionales de la función `boot` se pueden llevar 
a cabo modificaciones del bootstrap uniforme (Capítulo \@ref(modunif)). 
Algunos de estos parámetros son los siguientes:

* `strata`: permite realizar remuestreo estratificado estableciendo este parámetro
  como un vector numérico o factor que defina los grupos.

* `sim = c("ordinary" , "parametric", "balanced", "permutation", "antithetic")`:
  permite establecer distintos tipos de remuestreo. 
  Por defecto es igual a `"ordinary"` que se corresponde con el bootstrap uniforme,
  descrito anteriormente. Entre el resto de opciones destacaríamos 
  `sim = "permutation"`, que permite realizar contrastes de
  permutaciones (remuestreo sin reemplazamiento), y `sim = "parametric"`,
  que permite realizar bootstrap paramétrico (Sección \@ref(modunif-boot-par)). 
  En este último caso también habrá que establecer los parámetros `ran.gen` y
  `mle`, y la función `statistics` no empleará el segundo parámetro de índices.

* `ran.gen`: función que genera los datos. El primer argumento será el conjunto de datos
  original y el segundo un vector de parámetros adicionales 
  (normalmente los valores de los parámetros de la distribución).

* `mle`: parámetros de la distribución (típicamente estimados por máxima verosimilitud)
  o parámetros adicionales para `ran.gen` ó `statistics`.

Además hay otros parámetros para el procesamiento en paralelo: `parallel = c("no", "multicore", "snow")`, `ncpus`, `cl`. 
En el Apéndice \@ref(intro-hpc) se incluye una pequeña introducción al procesamiento en paralelo y se muestran algunos ejemplos sobre el uso de estos parámetros.
También se puede consultar la ayuda de la función `boot()` (`?boot`).

El paquete `boot` también incluye otras funciones que implementan métodos boostrap para otros tipos de datos, como la función `tsboot()` para series de tiempo (ver p.e. Cao y Fernández-Casal, 2022, [Capítulo 9](https://rubenfcasal.github.io/book_remuestreo/bootdep.html)) o la función `censboot()` para datos censurados (ver p.e. Cao y Fernández-Casal, 2022, [Capítulo 8](https://rubenfcasal.github.io/book_remuestreo/bootcen.html)).

Finalmente destacar que hay numerosas extensiones implementadas en otros paquetes utilizando el paquete `boot` (ver *Reverse dependencies* en la [web de CRAN](https://cran.r-project.org/package=boot)).
Por ejemplo en la Sección \@ref(boot-reg) se ilustrará el uso de la función `Boot()` del paquete `car` para hacer inferencia sobre modelos de regresión.


### Ejemplo: Bootstrap uniforme multidimensional {#boot-unif-multi}

Como ya se mostró en las Secciones \@ref(boot-unif) y \@ref(intro-paquetes) podemos implementar el bootstrap uniforme en el caso multidimensional (denominado también *remuestreo de casos* o *bootstrap de las observaciones*) de modo análogo al unidimensional.

Como ejemplo realizamos el Ejercicio \@ref(exr:unif-multi) empleando el paquete `boot` para estudiar la correlación lineal entre `prestige` e `income` del conjunto de datos `Prestige`:


```r
library(boot)

statistic <- function(data, i){
  remuestra <- data[i, ]
  cor(remuestra$income, remuestra$prestige)
}

set.seed(1)
B <- 1000
res.boot <- boot(Prestige, statistic, R = B)
res.boot
```

```
 ## 
 ## ORDINARY NONPARAMETRIC BOOTSTRAP
 ## 
 ## 
 ## Call:
 ## boot(data = Prestige, statistic = statistic, R = B)
 ## 
 ## 
 ## Bootstrap Statistics :
 ##      original      bias    std. error
 ## t1* 0.7149057 0.006306905  0.04406473
```

```r
plot(res.boot)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.9\linewidth]{09-Bootstrap_files/figure-latex/plot-boot-multi-1} 

}

\caption{Gráficos de diagnóstico para la distribución bootstrap del coeficiente de correlación.}(\#fig:plot-boot-multi)
\end{figure}

En este caso podemos observar que la distribución bootstrap del estimador es asimétrica, por lo que asumir que su distribución es normal podría no ser adecuado (por ejemplo para la construcción de intervalos de confianza, que se tratarán en la Sección \@ref(icboot-trans)).

Como comentario final, nótese que en principio el paquete boot está diseñado para obtener réplicas bootstrap de un estimador, por lo que si lo que nos interesa es emplear otro estadístico habría que construirlo a partir de ellas (como hacen otras funciones secundarias como `boot.ci()`).
Por ejemplo, si queremos emplear el estadístico $R = \hat \theta - \theta$
(bootstrap percentil básico o natural), podemos obtener la correspondiente distribución bootstrap (aproximada por Monte Carlo) con el siguiente código:


```r
estadistico_boot <- res.boot$t - res.boot$t0 
hist(estadistico_boot)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.75\linewidth]{09-Bootstrap_files/figure-latex/plot-boot-basico-1} 

}

\caption{Distribución bootstrap del estadístico percentil básico para el coeficiente de correlación.}(\#fig:plot-boot-basico)
\end{figure}

Podemos emplear la distribución empírica del estadístico bootstrap $R^{\ast} = \hat \theta^{\ast} - \hat \theta$ para aproximar la característica de interés de la distribución en el muestreo de $R = \hat \theta - \theta$.
Por ejemplo, para aproximar $\psi(u) = P\left( R \leq u \right)$ podríamos emplear la frecuencia relativa: 
$$\hat{\psi}_{B}(u) =
\frac{1}{B}\sum_{i=1}^{B}\mathbf{1}\left\{ R^{\ast (i)}\leq u\right\}.$$


```r
u <- 0
sum(estadistico_boot <= u)/B
```

```
 ## [1] 0.427
```

```r
# Equivalentemente:
mean(estadistico_boot <= u)
```

```
 ## [1] 0.427
```


