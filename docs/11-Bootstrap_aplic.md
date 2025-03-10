# Aplicaciones del remuestreo {#boot-aplic}

<!-- Capítulo \@ref(boot-aplic) -->



<!-- 
---
title: "Aplicaciones del remuestreo"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    pandoc_args: ["--number-offset", "10,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("11-Bootstrap_aplic.Rmd")
knitr::purl("11-Bootstrap_aplic.Rmd", documentation = 2)
knitr::spin("11-Bootstrap_aplic.R",knit = FALSE)

ver p.e. Cao y Fernández-Casal, 2021, [Sección 3.6](https://rubenfcasal.github.io/book_remuestreo/validez-de-la-aproximaci%C3%B3n-bootstrap.html)
-->

Como se ha visto en capítulos anteriores, los métodos bootstrap permiten aproximar la distribución en el muestreo de un estadístico.
En particular nos pueden interesar determinadas características de esta distribución, como el sesgo y la precisión de un estimador (Sección \@ref(sesgo-prec)).
Otro ejemplo es la construcción de intervalos de confianza, donde las características de interés son los cuantiles del estadístico pivotal empleado (Sección \@ref(boot-ic)).
También los contrastes de hipótesis, donde interesa aproximar cuantiles de la distribución muestral del estadístico del contraste bajo la hipótesis nula (Sección \@ref(boot-test)).


## Estimación del sesgo y la precisión de un estimador {#sesgo-prec}

Como ya se comentó, una de las principales aplicaciones de los métodos bootstrap es la aproximación del sesgo y de la varianza de un estimador.
Al igual que en capítulos anteriores, supongamos que $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$ es una m.a.s. de una población con distribución $F$ y que tenemos interés en realizar inferencia sobre un parámetro de la población $\theta =\theta \left( F \right)$ empleando un estimador $\hat{\theta}=T\left( \mathbf{X} \right)$.

Podemos considerar el estadístico
$$R=R\left( \mathbf{X}, F \right) = T\left( \mathbf{X} \right) - \theta \left( F \right) = \hat{\theta} - \theta,$$ 
para aproximar características de la distribución del estimador, ya que:
$$\begin{aligned}
Sesgo\left( \hat{\theta} \right) &= E\left( \hat{\theta}-\theta \right)
=E\left( R \right), \\
Var\left( \hat{\theta} \right) &= Var\left( \hat{\theta}-\theta \right)
=Var\left( R \right), \\
MSE\left( \hat{\theta} \right) &= E\left[ \left( \hat{\theta}-\theta \right)
^2\right] =E\left( R^2 \right).
\end{aligned}$$
Para ello, como se mostró en capítulos anteriores, consideraríamos una aproximación $\hat{F}$ de la distribución de probabilidad (por ejemplo, $\hat{F}=F_n$ si se considera que el bootstrap uniforme es adecuado) y emplearíamos Monte Carlo para aproximar la correspondiente distribución bootstrap:

1. Para cada $i=1,\ldots ,n$ generar $X_i^{\ast}$ a partir de
$\hat{F}$ y obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast}, \ldots ,X_n^{\ast} \right)$.

2. Calcular $R^{\ast}=R\left( \mathbf{X}^{\ast},\hat{F} \right)
= T\left( \mathbf{X}^{\ast} \right) -\theta \left( \hat{F} \right) =
\hat{\theta}^{\ast} - \hat{\theta}$.

3. Repetir $B$ veces los pasos 1-2 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$.

4. Utilizar las réplicas bootstrap para aproximar las características de interés de la distribución en el muestreo de $R$. 

    * Estimación bootstrap del sesgo:
      $$Sesgo^{\ast}\left( \hat{\theta}^{\ast} \right) =\bar{R}^{\ast}=\frac{1}{B}\sum_{b=1}^{B}R^{\ast (b)}.$$

    * Estimación bootstrap de la varianza:
      $$Var^{\ast}\left( \hat{\theta}^{\ast} \right) =\frac{1}{B} \sum_{b=1}^{B}\left( R^{\ast (b)}-\bar{R}^{\ast} \right)^2.$$

    * Estimación bootstrap del error cuadrático medio: 
      $$MSE^{\ast}\left( \hat{\theta}^{\ast} \right) =\frac{1}{B}\sum_{b=1}^{B}R^{\ast (b) 2}.$$


::: {.example #boot-simul-sesgo-var name="media y media truncada"}
<br>

Como ejemplo consideramos la aproximación del sesgo y de la varianza de la media y la media truncada al 10\% como estimadores de la media teórica del tiempo de vida de microorganismos mediante bootstrap uniforme (ver Ejercicio \@ref(exr:boot-simul-mediana)). 


```r
library(boot)
muestra <- simres::lifetimes

statistic <- function(data, i){
  remuestra <- data[i]
  c(mean(remuestra), mean(remuestra, trim = 0.1))
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
##     original    bias    std. error
## t1*  0.80533 0.0031733     0.15833
## t2*  0.75823 0.0119027     0.17394
```

Lamentablemente la función `print.boot()` calcula las aproximaciones bootstrap del sesgo y de la precisión pero no las almacena. 
En el caso más simple podríamos obtenerlas con el siguiente código:


```r
op <- with(res.boot, cbind(
  t0, apply(t, 2, mean, na.rm = TRUE) -  t0,
  apply(t, 2, sd, na.rm = TRUE)
  ))
rownames(op) <- paste0("t", 1:ncol(res.boot$t), "*")
colnames(op) <- c("original", "bias  ", " std. error")
op
```

```
##     original    bias    std. error
## t1*  0.80533 0.0031733     0.15833
## t2*  0.75823 0.0119027     0.17394
```

:::

::: {.exercise #unif-multi-sesgo-var}
<br>

Como continuación del ejemplo mostrado en la Sección \@ref(boot-unif-multi), emplear el paquete `boot` para aproximar el sesgo y la precisión del coeficiente de correlación lineal $r$ entre `prestige` e `income` del conjunto de datos `Prestige` mediante bootstrap uniforme multidimensional.
Comparar los resultados con los obtenidos mediante la aproximación asintótica normal, que como se indicó en el Ejercicio \@ref(exr:unif-multi), consideraría que el sesgo es nulo y como estimación de su varianza:
$$\widehat{Var}(r) = \frac{1 - r^2}{n - 2}.$$

::: 

---

Como también se comentó en la introducción del Capítulo \@ref(bootstrap), se pueden emplear otros tipos de remuestreo para la aproximación del sesgo y varianza de un estimador. 
El más conocido es el *jackknife*, que es uno de los métodos de remuestreo más antiguos (propuesto inicialmente por Quenouille, 1949).
De hecho el bootstrap surgió (Efron, 1979) como una alternativa a este tipo de remuestreo.
En el jackknife se consideran las $n$ remuestras obtenidas al ir eliminando cada una de las observaciones:
$$\mathbf{X}^{\ast}= \mathbf{X}_{(i)}=\left( X_1,\ldots ,X_{i-1},X_{i+1},\ldots
,X_n \right) \text{,  } i=1,\ldots ,n.$$
Para cada una de estas remuestras se obtienen las correspondientes réplicas del estadístico
$$\hat \theta_{(i)} = T \left( \mathbf{X}_{(i)} \right) \text{,  } i=1,\ldots ,n.$$
A partir de las cuales se aproxima el sesgo y la varianza, considerando un factor de elevación $n-1$ para que sean insesgadas:
$$\begin{aligned}
Sesgo_{jackk}^{\ast}\left( \hat{\theta}^{\ast} \right) 
&= \left( n-1 \right)\left( \overline{\theta_{(\cdot)}}-\hat{\theta} \right) 
= \frac{n-1}{n}\sum_{i=1}^{n}\left( \hat \theta_{(i)} - \hat{\theta} \right), \\
Var_{jackk}^{\ast}\left( \hat{\theta}^{\ast} \right) 
&= \frac{n-1}{n} \sum_{i=1}^{n}\left[ \hat \theta_{(i)} 
- \overline{\theta_{(\cdot)}}\right]^2,
\end{aligned}$$
donde $\overline{\theta_{(\cdot)}} = \frac{1}{n}\sum_{j=1}^{n}\hat \theta_{(j)}$.
Para más detalles, ver [Sección 2.2](https://rubenfcasal.github.io/book_remuestreo/jackknife.html) de Cao y Fernández-Casal (2021).


## Intervalos de confianza bootstrap {#boot-ic}

En esta sección consideraremos el problema de construcción, mediante bootstrap, de un intervalo de confianza bilateral, con nivel de confianza $1-\alpha$, para un parámetro $\theta$ de la distribución $F$.
Una vez elegido el método bootstrap adecuado, teniendo en cuenta la información disponible en el contexto del que se trate, otro aspecto importante es el método para la construcción del intervalo de confianza bootstrap de forma que la probabilidad de cobertura sea lo más parecida posible al nivel nominal $1-\alpha$.

Las diferencias entre los distintos métodos dependen del estadístico $R$ empleado y de las suposiciones sobre su distribución.
Como se comentó en la Sección \@ref(intro-pkgboot), la función `boot.ci()` del paquete `boot` permite construir distintos tipos de intervalos de confianza dependiendo del parámetro `type`.
En el Ejemplo \@ref(exm:media-dt-desconocida-boot) se ilustra la obtención de estimaciones por intervalo de confianza para la media empleando los distintos métodos bajo bootstrap uniforme (en el Capítulo \@ref(modunif) se incluyen ejemplos adicionales empleando bootstrap paramétrico y suavizado).

En esta sección se describirán brevemente los distintos métodos implementados en la función `boot.ci()`.
Para un tratamiento más detallado, incluyendo los órdenes de los errores de cobertura, ver por ejemplo el [Capítulo 4](https://rubenfcasal.github.io/book_remuestreo/icboot.html) de Cao y Fernández-Casal (2021) o el Capítulo 5 de Davison y Hinkley (1997).


### Aproximación normal {#boot-ic-norm}

Este método emplea las aproximaciones bootstrap del sesgo $Sesgo^{\ast}\left( \hat{\theta}^{\ast} \right)$ y de la varianza $Var^{\ast}\left( \hat{\theta}^{\ast} \right)$, y asume que la distribución del correspondiente estadístico studentizado es una normal estándar
$$\frac{\hat{\theta} - Sesgo^{\ast}\left( \hat{\theta}^{\ast} \right) - \theta}{\sqrt{Var^{\ast}\left( \hat{\theta}^{\ast} \right)}} \underset{aprox}{\sim }\mathcal{N}\left( 0, 1 \right).$$
De esta forma se obtiene la estimación por intervalo de confianza:
$$\hat{I}_{norm}=\left( \hat{\theta} - Sesgo^{\ast}\left( \hat{\theta}^{\ast} \right) - z_{1-\alpha /2}\sqrt{Var^{\ast}\left( \hat{\theta}^{\ast} \right)},\hat{\theta} - Sesgo^{\ast}\left( \hat{\theta}^{\ast} \right) + z_{1 - \alpha /2}\sqrt{Var^{\ast}\left( \hat{\theta}^{\ast} \right)} \right).$$

Podemos obtener este intervalo de confianza estableciendo `type = "norm"` (o `type = "all"`) en la llamada a la función `boot.ci()` (ver Ejemplo \@ref(exm:media-dt-desconocida-boot)).
        

### Método percentil directo {#boot-ic-perc}

Este método se basa en la construcción del intervalo de confianza, mediante bootstrap, empleando como estadístico el estimador $$R = \hat{\theta}.$$

Una vez elegido el método de remuestreo, empleando un estimador, $\hat{F}\,$, de la
distribución poblacional, $F$, la distribución en el muestreo de $R = \hat{\theta}$ se aproxima directamente mediante la distribución bootstrap de $R^{\ast}= \hat{\theta}^{\ast}$.
A partir de las réplicas bootstrap del estimador aproximamos los cuantiles $x_{\alpha /2}$ y $x_{1-\alpha /2}$ (denotando por $x_{\beta }$ el valor verificando $P^{\ast}\left( R^{\ast }\leq x_{\beta } \right) =\beta$), de forma que 
$$\begin{aligned}
1-\alpha &= 1-\frac{\alpha }{2}-\frac{\alpha }{2}
= P^{\ast}\left( \hat{\theta}^{\ast}<x_{1-\alpha /2} \right) - P^{\ast}\left( \hat{\theta}^{\ast}\leq x_{\alpha /2} \right) \\
&= P^{\ast}\left( x_{\alpha /2}<\hat{\theta}^{\ast}<x_{1-\alpha /2} \right),
\end{aligned}$$
y asumimos que esto aproxima lo que ocurre con la distribución poblacional
$$P\left( x_{\alpha /2} < \hat{\theta} < x_{1-\alpha /2} \right) \approx 1-\alpha.$$
De donde se obtiene el intervalo de confianza bootstrap calculado 
por el método percentil directo
$$\hat{I}_{perc}=\left( x_{\alpha /2}, x_{1-\alpha /2}  \right).$$

Una ventaja de los intervalos construidos con este método es que son invariantes frente a transformaciones del estimador (en el caso de que fuese más adecuado trabajar en otra escala, no sería necesario conocer la transformación).
Sin embargo, como se comentó en la Sección \@ref(boot-intro), la precisión puede verse seriamente afectada en el caso de estimadores sesgados.

Podemos obtener este intervalo de confianza estableciendo `type = "perc"` (o `type = "all"`) en la llamada a la función `boot.ci()` (ver Ejemplo \@ref(exm:media-dt-desconocida-boot)).

### Método percentil básico {#boot-ic-basic}

En este método se emplea como estadístico el estimador centrado (no estandarizado)
$$R = \hat{\theta}-\theta.$$
De forma análoga, la distribución en el muestreo de $R$ se aproxima mediante la distribución bootstrap de
$$R^{\ast}= \hat{\theta}^{\ast}-\theta \left( \hat{F} \right) = \hat{\theta}^{\ast}-\hat{\theta}.$$
A partir de las réplicas bootstrap del estadístico se aproximan los cuantiles $x_{\alpha /2}$ y $x_{1-\alpha /2}$ tales que
$$1-\alpha = P^{\ast}\left( x_{\alpha /2}<R^{\ast}<x_{1-\alpha /2} \right),$$
tomándolo como aproximación de lo que ocurre con la distribución poblacional
$$\begin{aligned}
1-\alpha &\approx P\left( x_{\alpha /2}<R<x_{1-\alpha /2} \right) \\
&= P\left( x_{\alpha /2} < \hat{\theta}-\theta < x_{1-\alpha /2} \right) \\
&= P\left( \hat{\theta} - x_{1-\alpha /2} < \theta <\hat{\theta} -x_{\alpha /2} \right).
\end{aligned}$$
De donde se obtiene el intervalo de confianza bootstrap calculado 
por el método percentil básico
$$\hat{I}_{basic}=\left( \hat{\theta} - x_{1-\alpha /2},\hat{\theta} - x_{\alpha /2} \right).$$

Podemos obtener este intervalo de confianza estableciendo `type = "basic"` (o `type = "all"`) en la llamada a la función `boot.ci()` (ver Ejemplo \@ref(exm:media-dt-desconocida-boot)).


### Método percentil-*t* {#boot-ic-stud}

Este método bootstrap, construye un intervalo de confianza bootstrap a partir del estadístico studentizado:
$$R = \frac{\hat \theta - \theta}{\sqrt{\widehat{Var}(\hat \theta)}}.$$
Procediendo de modo análogo a los casos anteriores, se aproximan los cuantiles $x_{\alpha /2}$ y $x_{1-\alpha /2}$ tales que
$$1-\alpha = P^{\ast}\left( x_{\alpha /2}<R^{\ast}<x_{1-\alpha /2} \right),$$
a partir de los cuales se obtiene el intervalo de confianza bootstrap calculado 
por el método percentil-*t* (o percentil studentizado)
$$\hat{I}_{stud}=\left( \hat{\theta} - x_{1-\alpha /2}\sqrt{\widehat{Var}(\hat \theta)},\hat{\theta} - x_{\alpha /2}\sqrt{\widehat{Var}(\hat \theta)} \right).$$

Si uno de los componentes del vector de estadísticos (por defecto el segundo) es una estimación de la varianza del estimador (por defecto el primer componente), podemos obtener este intervalo de confianza estableciendo `type = "stud"` (o `type = "all"`) en la llamada a la función `boot.ci()` (ver Ejemplo \@ref(exm:media-dt-desconocida-boot)).
En caso de que el primer y segundo componente del vector de estadísticos no sean el estimador y su varianza estimada, respectivamente, habrá que emplear el argumento `index`, que debe ser un vector de longitud 2 con las posiciones correspondientes.

Hay una variante de este método, denominada percentil-*t* simetrizado, en la que se asume que la distribución del estadístico es simétrica (aunque no está implementado en `boot.ci()`).
Asumiendo que esta suposición es correcta, podemos calcular los cuantiles de forma más eficiente (ya que dispondríamos del doble de información sobre las colas de la distribución).
En lugar de tomar cuantiles que dejen colas iguales ($\frac{\alpha }{2}$) a la izquierda y a la derecha, respectivamente, se considera el valor $x_{1-\alpha }$ que verifica $P^{\ast}\left( \left\vert R^{\ast}\right\vert \leq x_{1-\alpha } \right) =1-\alpha$. 
Así se obtiene el intervalo de confianza bootstrap
$$\hat{I}_{simstud}=\left( \hat{\theta} - x_{1-\alpha}\sqrt{\widehat{Var}(\hat \theta)},\hat{\theta} + x_{1-\alpha}\sqrt{\widehat{Var}(\hat \theta)} \right).$$

::: {.example #media-ic-stu-sim name="IC bootstrap para la media mediante el método percentil-*t* simetrizado"}
<br> 

Continuando con el Ejemplo \@ref(exm:media-dt-desconocida) de inferencia sobre la media con varianza desconocida, podríamos obtener una estimación por intervalo de confianza del tiempo de vida medio de los microorganismos empleando el método bootstrap percentil-*t* simetrizado con el siguiente código:


```r
muestra <- simres::lifetimes
n <- length(muestra)
alfa <- 0.05
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
  estadistico_boot[k] <- sqrt(n) * abs(x_barra_boot - x_barra)/cuasi_dt_boot
}

# Aproximación bootstrap del pto crítico
pto_crit <- quantile(estadistico_boot, 1 - alfa)

# Construcción del IC
ic_inf_boot <- x_barra - pto_crit * cuasi_dt/sqrt(n)
ic_sup_boot <- x_barra + pto_crit * cuasi_dt/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
##    2.5%   97.5% 
## 0.43347 1.17719
```

::: 


### Método BCa {#boot-ic-bca}

El método $BCa$ (bias-corrected and accelerated) propuesto por Efron (1987) considera una transformación de forma que la distribución se aproxime a la normalidad, construye el intervalo en esa escala asumiendo normalidad y transforma el resultado a la escala original empleando la distribución bootstrap.
El intervalo obtenido es de la forma:
$$\hat{I}_{bca}=\left( x_{\alpha /2}, x_{1-\alpha /2}  \right),$$
donde        
$$x_u = \hat G^{-1}\left(\Phi\left(z + \frac{z + z_u}{1-a(z+z_u)}\right)  \right),$$
siendo $\hat G$ la distribución empírica de $\hat{\theta}^{\ast}$, $\Phi(z)$ la función de distribución de la normal estándar, $z_u = \Phi^{-1}(u)$ el correspondiente cuantil y:

- $z = \Phi^{-1}(\hat G(\hat\theta))$ un factor de corrección de sesgo.
- $a$ la denominada constante aceleradora (o corrección de asimetría), que suele ser aproximada mediante jackknife.

Podemos obtener este intervalo de confianza estableciendo `type = "bca"` (o `type = "all"`) en la llamada a la función `boot.ci()` (ver Ejemplo \@ref(exm:media-dt-desconocida-boot)).
Para más detalles ver Sección 5.3.2 de Davison y Hinkley (1997).

<!-- 
Normal score transformation (NST)
Rank-based inverse normal transformation
https://rdrr.io/cran/rcompanion/man/blom.html
-->



::: {.exercise #boot-unif-multi-ic}
<br>

Como continuación del ejemplo mostrado en la Sección \@ref(boot-unif-multi), y de los Ejercicios \@ref(exr:unif-multi) y \@ref(exr:unif-multi-sesgo-var), emplear el paquete `boot` para obtener estimaciones por intervalo de confianza del coeficiente de correlación lineal $r$ entre `prestige` e `income` del conjunto de datos `Prestige` (mediante bootstrap uniforme multidimensional).
En el caso del método percentil-*t*, como se indicó en el Ejercicio \@ref(exr:unif-multi), considerar el estimador de la varianza:
$$\widehat{Var}(r) = \frac{1 - r^2}{n - 2}.$$

Comparar los resultados con la aproximación paramétrica implementada en la función `cor.test` y descrita en la siguiente sección.

::: 


### Ejemplo: IC bootstrap para el coeficiente de correlación {#icboot-trans}

Supongamos de nuevo que queremos estudiar la correlación entre dos variables $X$ e $Y$ a partir del coeficiente de correlación lineal de Pearson:
$$\rho =\frac{ Cov \left( X, Y \right) }
{ \sigma \left( X \right) \sigma \left( Y \right) },$$
empleando como estimador el coeficiente de correlación muestral:
$$r=\frac{\sum_{i=1}^{n}(x_i-\overline{x})(y_i-\overline{y})}
{\sqrt{ \sum_{i=1}^{n}(x_i-\overline{x})^{2}} 
\sqrt{\sum_{i=1}^{n}(y_i-\overline{y})^{2}}}.$$

Para realizar inferencias sobre el coeficiente de correlación, como aproximación más simple, se puede considerar que la distribución muestral de $r$ es aproximadamente normal (ver Ejercicio \@ref(exr:unif-multi)) y emplear el estadístico:

\begin{equation} 
\frac{r -\rho}{\sqrt{\frac{1 - r^2}{n - 2}}} \underset{aprox}{\sim } t_{n-2}
(\#eq:cor-t)
\end{equation} 

<!-- \mathcal{t}_{n-2} error en LaTeX-->

Pero esta aproximación solo sería válida en el caso de muestras grandes (o si la distribución bivariante de $(X, Y)$ es aproximadamente normal) cuando la correlación entre las variables es débil o moderada. 
En caso contrario la distribución muestral de $r$ puede ser muy asimétrica y los resultados obtenidos con el estadístico anterior no ser muy adecuados (esto concuerda con lo observado en la Sección \@ref(boot-unif-multi), al emplear bootstrap uniforme multidimensional para hacer inferencia sobre $R = r -\rho$).
Para evitar este problema se suelen obtener intervalos de confianza para $\rho$ empleando la transformación $Z$ de Fisher (1915):
$$Z = \frac{1}{2}\ln \left( \frac{1+r}{1-r} \right) = \operatorname{arctanh}(r),$$
que es una transformación (aprox.) normalizadora y estabilizadora de la varianza.
Suponiendo que $(X, Y)$ es normal bivariante y que hay independencia entre las observaciones:
$$Z \sim \mathcal{N}\left( \frac{1}{2}\ln \left( \frac{1+\rho}{1-\rho} \right), \frac{1}{n-3} \right).$$
El intervalo de confianza asintótico se obtiene empleando la aproximación normal tradicional en la escala $Z$ y aplicando posteriormente la transformación inversa:
$$r = \frac{\exp(2Z)-1}{\exp(2Z)+1} = \operatorname{tanh}(Z).$$

Esta aproximación está implementada en la función `cor.test()` del paquete base `stat` de R^[Se puede obtener el código tecleando en la consola `stats:::cor.test.default`.], además de que también realiza el contraste $H_0: \rho = 0$ empleando el estadístico \@ref(eq:cor-t).

Continuando con el ejemplo de la Sección \@ref(boot-unif-multi) (y de los Ejercicios \@ref(exr:unif-multi), \@ref(exr:unif-multi-sesgo-var) y \@ref(exr:boot-unif-multi-ic)), para obtener un intervalo de confianza para el coeficiente de correlación lineal entre las variables `income` y `prestige` del conjunto de datos `Prestige`, podríamos emplear el siguiente código:


```r
data(Prestige, package="carData")
# with(Prestige, cor.test(income, prestige))
cor.test(Prestige$income, Prestige$prestige)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  Prestige$income and Prestige$prestige
## t = 10.2, df = 100, p-value <2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.60447 0.79838
## sample estimates:
##     cor 
## 0.71491
```

También es de esperar que mejore la precisión de los intervalos de confianza bootstrap si se emplea una transformación que estabilice la varianza del estimador, especialmente en el caso del método basado en la aproximación normal y del bootstrap percentil básico. 
La función `boot.ci()` del paquete `boot` permite obtener intervalos de confianza calculados en una escala transformada del estadístico, mediante los parámetros:

- `h`: función vectorial que define la transformación. 
  Los intervalos se calculan en la escala de $h(t)$ y se aplica la función inversa (si se especifica) para transformarlos a la escala original.

- `hinv`: (opcional) función inversa de la transformación 
  (si no se especifica solo se calculan los intervalos en la escala transformada). 

- `hdot`: (opcional en el método percentil o básico) función derivada de la transformación 
  (empleada por algunos métodos para aproximar la varianza en la escala transformada mediante el método delta).

Por ejemplo, para considerar la transformación $Z$ de Fisher en este caso, se podría emplear el siguiente código:

```r
library(boot)

statistic <- function(data, i){
  remuestra <- data[i, ]
  cor(remuestra$income, remuestra$prestige)
}

set.seed(1)
res.boot <- boot(Prestige, statistic, R = 1000)

h <- function(t) atanh(t)
hdot <- function(t) 1/(1 - t^2)
hinv <- function(t) tanh(t)

# boot.ci(res.boot, type = "norm", h = h)
boot.ci(res.boot, type = "norm", h = h, hdot = hdot, hinv = hinv)
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = res.boot, type = "norm", h = h, hdot = hdot, 
##     hinv = hinv)
## 
## Intervals : 
## Level      Normal        
## 95%   ( 0.6016,  0.7858 )  
## Calculations on Transformed Scale;  Intervals on Original Scale
```

Esto sería en principio preferible a trabajar en la escala original, ya que la distribución bootstrap en la escala transformada se aproximaría más a la normalidad:


```r
ht <- h(res.boot$t)
hist(ht, freq = FALSE, breaks = "FD", main = "")
curve(dnorm(x, mean=mean(ht), sd=sd(ht)), lty = 2, add = TRUE)
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/icboot-trans-plot-1.png" alt="Distribución bootstrap en la escala transformada." width="70%" />
<p class="caption">(\#fig:icboot-trans-plot)Distribución bootstrap en la escala transformada.</p>
</div>


## Contrastes de hipótesis bootstrap {#boot-test}

El objetivo de los contrastes de hipótesis es, a partir de la información 
que proporciona una muestra, decidir (tratando de controlar el riesgo de 
equivocarse al no disponer de toda la información)
entre dos hipótesis sobre alguna característica de interés de la población: 
hipótesis nula ($H_{0}$) e hipótesis alternativa ($H_{1}$).

Entre los distintos tipos de contrastes de hipótesis (e.g. paramétricos, 
no paramétricos, ...), nos centraremos principalmente en los contrastes 
de bondad de ajuste. En este caso interesará distinguir principalmente 
entre hipótesis nulas simples (especifican un único modelo) y compuestas 
(especifican un conjunto/familia de modelos).

Para realizar el contraste se emplea un estadístico $D\left( X_1,\ldots ,X_n;H_0\right)$,
que mide la discrepancia entre la muestra observada y la hipótesis nula, 
con distribución conocida (o que se puede aproximar) bajo $H_0$.
Por ejemplo, en el caso de una hipótesis nula paramétrica 
es habitual emplear un estadístico studentizado de la forma:
$$D\left( X_1,\ldots ,X_n;H_0\right) =
\frac{\hat{\theta}-\theta _0}{\hat\sigma_{\hat\theta}}$$
(o algún tipo de razón de verosimilitudes).

La regla de decisión depende de la hipótesis alternativa y 
del riesgo asumible al rechazar $H_0$ siendo cierta:
$$P\left( \text{rechazar }H_0\mid H_0\text{ cierta}\right) =\alpha,$$ 
denominado nivel de significación. 
Se determina una región de rechazo (RR) a partir de los valores que tiende 
a tomar el estadístico cuando $H_1$ es cierta, 
de forma que^[Aunque cuando la hipótesis nula es compuesta:
$P\left( D\in RR \mid H_0\text{ cierta}\right) \leq \alpha$.]:
$$P\left( D\in RR \mid H_0\text{ cierta}\right) =\alpha.$$
Se rechaza la hipótesis nula cuando el valor observado del 
estadístico $\hat{d}=D\left( x_1,\ldots ,x_n;H_0\right)$ pertenece
a la región de rechazo.

Para medir el nivel de evidencia en contra de $H_0$ se emplea el
$p$-valor del contraste (también denominado valor crítico o 
tamaño del contraste), el menor valor del nivel 
de significación para el que se rechaza $H_0$ 
(que se puede interpretar también como la 
probabilidad de obtener una discrepancia mayor o igual que
$\hat{d}$ cuando $H_0$ es cierta).

El cálculo del $p$-valor dependerá por tanto de la hipótesis altervativa.
Por ejemplo, si el estadístico del contraste tiende a tomar valores
grandes cuando $H_0$ es falsa (contraste unilateral derecho):
$$p = P\left( D \geq \hat{d} \middle| H_0\right).$$
En otros casos (contrastes bilaterales) hay evidencias en contra de 
$H_0$ si el estadístico toma valores significativamente grandes o pequeños.
En estos casos la distribución del estadístico del contraste bajo $H_0$ 
suele ser simétrica en torno al cero, por lo que:
$$p = 2P\left( D \geq \vert \hat{d} \vert \middle| H_0 \right).$$
Pero si esta distribución es asimétrica:
$$p = 2 \min \left\{ P\left( D \leq \hat{d} \middle| H_0 \right),
P\left( D \geq \hat{d} \middle| H_0\right) \right\}.$$

La regla de decisión a partir del $p$-valor es siempre la misma.
Rechazamos $H_0$, al nivel de significación $\alpha$, si $p \leq \alpha$,
en cuyo caso se dice que el contraste es estadísticamente significativo
(rechazamos $H_0$ con mayor seguridad cuanto más pequeño es el $p$-valor).
Por tanto, la correspondiente variable aleatoria $\mathcal{P}$ debería verificar:
$$P\left( \mathcal{P} \leq \alpha \middle| H_0\right)= \alpha.$$
Es decir, la distribución del $p$-valor bajo $H_0$ debería ser $\mathcal{U}(0,1)$
(si la distribución del estadístico del constrate es continua).


En los métodos tradicionales de contrastes de hipótesis se conoce o se 
puede aproximar la distribución del estadístico del contraste bajo $H_0$.
Muchas de estas aproximaciones están basadas en resultados asintóticos 
y pueden no ser adecuadas para tamaños muestrales pequeños.
En ese caso, o si no se dispone de estas herramientas, 
se puede recurrir a métodos de remuestreo para aproximar el $p$-valor.
Uno de los procedimientos más antiguos es el denominado
*contraste de permutaciones* (Fisher, 1935; Pitman, 1937; Welch, 1937).
Aunque el bootstrap paramétrico y el semiparamétrico son los 
procedimientos de remuestreo más empleados para aproximar 
la distribución del estadístico de contraste bajo la hipótesis nula.

La idea es obtener remuestras de una aproximación de la distribución del
estadístico bajo $H_0$. 
En el bootstrap paramétrico y semiparamétrico se estima la distribución 
de los datos bajo la hipótesis nula, $\hat{F}_0$, y se obtienen réplicas del 
estadístico a partir de remuestras de esta distribución (no sería adecuado 
emplear directamente la distribución empírica). 
En el caso de los contrastes de permutaciones las remuestras se obtienen 
directamente de los datos, remuestreando sin reemplazamiento los valores 
de la respuesta (y manteniendo fijas las covariables).

Finalmente, se emplean las réplicas bootstrap
del estadístico $d_1^{\ast},\ldots, d_B^{\ast}$ para aproximar el $p$-valor.
Por ejemplo, en el caso de un contraste unilateral en el que el estadístico del
contraste tiende a tomar valores grandes si la hipótesis nula es falsa,
se podría emplear como aproximación:
$$p_{boot} = \frac{1}{B}\#\left\{ d_i^{\ast} \geq \hat{d} \right\}.$$
Mientras que en el caso bilateral, asumiendo que la distribución del estadístico
no es necesariamente simétrica, habría que emplear:
$$p_{boot} = \frac{2}{B} \min \left(\#\left\{ d_i^{\ast} \leq \hat{d} \right\},
\#\left\{ d_i^{\ast} \geq \hat{d} \right\}\right).$$


### Contrastes bootstrap paramétricos {#contrastes-parametricos}

En los casos en los que la hipótesis nula especifica por completo la distribución
(hipótesis nula simple) o solo desconocemos los valores de algunos parámetros 
(hipótesis nula paramétrica compuesta) podemos emplear
bootstrap paramétrico para obtener las remuestras bootstrap de los datos 
(realmente en el primer caso se trataría de simulaciones Monte Carlo).
Siempre hay que tener en cuenta que las réplicas bootstrap del estadístico se 
deberían obtener empleando el mismo procedimiento utilizado en la muestra
(p.e. reestimando los parámetros si es el caso).

<!-- 
Pendiente: 
Referencia test ji-cuadrado 
Mover parte a sección 7.4.3 y reescribir

-->


::: {.example #ks-test-sim-exp name="Contraste de Kolmogorov-Smirnov"}
<br>

Se trata de un contraste de bondad de ajuste (similar a la prueba de 
Cramer-von Mises o a la de Anderson-Darling, implementadas en el paquete 
`goftest` de R, que son en principio mejores). 
A partir de $X_1,\ldots ,X_n$ m.a.s. de $X$ con función de distribución $F$, 
se pretende contrastar:
$$\left \{ 
\begin{array}{l}
H_0 : F = F_0 \\ 
H_1 : F \neq F_0 
\end{array}
\right. $$
siendo $F_0$ una función de distribución continua. 
El estadístico empleado para ello compara la función de distribución bajo 
$H_0$ ($F_0$) con la empírica ($F_n$):
$$\begin{aligned}
    D_n=&\sup_{x}|F_n(x)-F_0(x)| \\
    =&\max_{1 \leq i\leq n}\left \{
    |F_n(X_{(i)})-F_0(X_{(i)})|,|F_n(X_{(i-1)})-F_0(X_{(i)})|\right \} \\
    =&\max_{1 \leq i\leq n}\left \{ \frac{i}{n}-F_0(X_{(i)}), \ F_0(X_{(i)})-\frac{i-1}{n}\right \} \\
    =&\max_{1 \leq i\leq n}\left \{ D_{n,i}^{+},\ D_{n,i}^{-}\right \},
\end{aligned}$$
y su distribución bajo $H_0$ no depende $F_0$ (es de distribución libre), 
si $H_0$ es simple y $F_0$ es continua. 
Esta distribución está tabulada (para tamaños muestrales grandes se utiliza 
la aproximación asintótica) y se rechaza $H_0$ si el valor observado $d$ 
del estadístico es significativamente grande:
$$p = P \left( D_n \geq d \right) \leq \alpha.$$
Este método está implementado en la función `ks.test()` del paquete base de R:

```r
ks.test(x, y, ...)
```
donde `x` es un vector que contiene los datos, `y` es una función de distribución 
(o una cadena de texto que la especifica; también puede ser otro vector de datos 
para el contraste de dos muestras) y `...` representa los parámetros de la distribución.

Si $H_0$ es compuesta, el procedimiento habitual es estimar los parámetros desconocidos 
por máxima verosimilitud y emplear $\hat{F}_0$ en lugar de $F_0$. 
Sin embargo, al proceder de esta forma es de esperar que $\hat{F}_0$ se aproxime más 
que $F_0$ a la distribución empírica, por lo que los cuantiles de la distribución de 
$D_n$ pueden ser demasiado conservativos (los $p$-valores tenderán a ser mayores de 
lo que deberían) y se tenderá a aceptar la hipótesis nula. 

Para evitar este problema, en el caso de contrastar normalidad se desarrolló el test 
de Lilliefors, implementado en la función `lillie.test()` del paquete `nortest` 
(también hay versiones en este paquete para los métodos de Cramer-von Mises y
Anderson-Darling). Como ejemplo analizaremos el comportamiento de ambos métodos
para contrastar normalidad considerando 1000 pruebas con muestras de tamaño 30 de 
una $\mathcal{N}(0,1)$ (estudiaremos el *tamaño de los contrastes*). 



```r
# Valores iniciales
library(nortest)
set.seed(1)
nx <- 30
mx <- 0
sx <- 1
nsim <- 1000
# Realizar contrastes
pvalor.ks <- numeric(nsim)
pvalor.lil <- numeric(nsim)
for(isim in 1:nsim) {
  rx <- rnorm(nx, mx, sx)
  pvalor.ks[isim] <- ks.test(rx, "pnorm", mean(rx), sd(rx))$p.value
  pvalor.lil[isim] <- lillie.test(rx)$p.value
}
```

Bajo la hipótesis nula el $p$-valor debería de seguir una distribución uniforme,
por lo que podríamos generar el correspondiente histograma para estudiar el
tamaño del contraste. Alternativamente podríamos representar su función de 
distribución empírica, que se correspondería con la proporción de rechazos
para los distintos niveles de significación. 



```r
old.par <- par(mfrow=c(2, 2))
# Test de KS
  # Histograma
hist(pvalor.ks, freq=FALSE, main = "")
abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
  # Distribución empírica
curve(ecdf(pvalor.ks)(x), type = "s", lwd = 2, main = '', 
      ylab = 'Proporción de rechazos', xlab = 'Nivel de significación')
abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE)
# Test de Lilliefors
  # Histograma
hist(pvalor.lil, freq=FALSE, main = "")
abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
  # Distribución empírica
curve(ecdf(pvalor.lil)(x), type = "s", lwd = 2, main = '', 
      ylab = 'Proporción de rechazos',  xlab = 'Nivel de significación')
abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE)
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/ks-lil-plot-1.png" alt="Distribución del p-valor (izquierda) y tamaño (proporción de rechazos bajo la hipótesis nula; derecha), aproximados por Monte Carlo, para el contraste de Kolmogorov-Smirnov (arriba) y el de Lilliefors (abajo)." width="90%" />
<p class="caption">(\#fig:ks-lil-plot)Distribución del p-valor (izquierda) y tamaño (proporción de rechazos bajo la hipótesis nula; derecha), aproximados por Monte Carlo, para el contraste de Kolmogorov-Smirnov (arriba) y el de Lilliefors (abajo).</p>
</div>

```r
par(old.par)
```

En el caso del contraste de Kolmogorov-Smirnov (KS) se observa que el $p$-valor 
tiende a tomar valores grandes y por tanto se rechaza la hipótesis nula 
muchas menos veces de las que se debería.

En el caso de otras distribuciones se puede emplear bootstrap paramétrico para 
aproximar la distribución del estadístico del contraste.
Es importante recordar que el bootstrap debería imitar el procedimiento
empleado sobre la muestra, por lo que en este caso también habría que estimar
los parámetros en cada remuestra 
(en caso contrario aproximaríamos la distribución de $D_n$).

Por ejemplo, la siguiente función implementaría el contraste KS de
bondad de ajuste de una variable exponencial aproximando el
$p$-valor mediante bootstrap paramétrico:


```r
ks.exp.boot <- function(x, nboot = 10^3) {
  DNAME <- deparse(substitute(x))
  METHOD <- "Kolmogorov-Smirnov Test of pexp by bootstrap" 
  n <- length(x)
  RATE <- 1/mean(x)
  ks.exp.stat <- function(x, rate = 1/mean(x)) { # se estima el parámetro
    DMinus <- pexp(sort(x), rate=rate) - (0:(n - 1))/n
    DPlus <- 1/n - DMinus
    Dn = max(c(DMinus, DPlus))
  }  
  STATISTIC <- ks.exp.stat(x, rate = RATE) 
  names(STATISTIC) <- "Dn"
  # PVAL <- 0
  # for(i in 1:nboot) {
  #   rx <- rexp(n, rate = RATE)
  #   if (STATISTIC <= ks.exp.stat(rx)) PVAL <- PVAL + 1
  # }
  # PVAL <- PVAL/nboot
  # PVAL <- (PVAL + 1)/(nboot + 1) # Alternativa para aproximar el p-valor
  rx <- matrix(rexp(n*nboot, rate = RATE), ncol=n)
  PVAL <- mean(STATISTIC <= apply(rx, 1, ks.exp.stat))
  return(structure(list(statistic = STATISTIC, alternative = "two.sided", 
                   p.value = PVAL, method = METHOD, data.name = DNAME), 
                   class = "htest"))
}
```

Como ejemplo estudiaremos el caso de contrastar una distribución exponencial
considerando 500 pruebas con muestras de tamaño 30 de una $Exp(1)$ 
y 200 réplicas bootstrap (para disminuir el tiempo de computación).


``` r
# Valores iniciales
set.seed(1)
nx <- 30
ratex <- 1
nsim <- 500
# Realizar contrastes
pvalor.ks <- numeric(nsim)
pvalor.ks.boot <- numeric(nsim)
for(isim in 1:nsim) {
  rx <- rexp(nx, ratex)
  pvalor.ks[isim] <- ks.test(rx, "pexp", 1/mean(rx))$p.value
  pvalor.ks.boot[isim] <- ks.exp.boot(rx, nboot = 200)$p.value
}
# Generar gráficos
old.par <- par(mfrow=c(2, 2))
# Test de KS
  # Histograma
hist(pvalor.ks, freq=FALSE, main = "")
abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
  # Distribución empírica
curve(ecdf(pvalor.ks)(x), type = "s", lwd = 2, main = '', 
      ylab = 'Proporción de rechazos', xlab = 'Nivel de significación')
abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE)
# Contraste bootstrap paramétrico
  # Histograma
hist(pvalor.ks.boot, freq=FALSE, main = "")
abline(h=1, lty=2)   # curve(dunif(x,0,1), add=TRUE)
  # Distribución empírica
curve(ecdf(pvalor.ks.boot)(x), type = "s", lwd = 2, main = '', 
      ylab = 'Proporción de rechazos',  xlab = 'Nivel de significación')
abline(a=0, b=1, lty=2)   # curve(punif(x, 0, 1), add = TRUE)
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/ks-boot-plot-1.png" alt="Distribución del p-valor (izquierda) y tamaño (proporción de rechazos bajo la hipótesis nula; derecha), aproximados por  Monte Carlo, para el contraste de Kolmogorov-Smirnov (arriba) y el correspondiente contraste boostrap paramétrico (abajo)." width="90%" />
<p class="caption">(\#fig:ks-boot-plot)Distribución del p-valor (izquierda) y tamaño (proporción de rechazos bajo la hipótesis nula; derecha), aproximados por  Monte Carlo, para el contraste de Kolmogorov-Smirnov (arriba) y el correspondiente contraste boostrap paramétrico (abajo).</p>
</div>

``` r
par(old.par)
```

El estadístico de Kolmogorov-Smirnov `Dn = max(c(DMinus, DPlus))` tiene ventajas desde el
punto de vista teórico, pero puede no ser muy potente para detectar diferencias entre la
distribución bajo la hipótesis nula y la distribución de los datos. 
La ventaja de la aproximación por simulación es que no estamos atados a resultados teóricos
y podemos emplear el estadístico que se considere oportuno 
(la principal desventaja es el tiempo de computación). 
Por ejemplo, podríamos pensar en utilizar como estadístico la suma de los errores en
valor absoluto del correspondiente gráfico PP, y solo habría que cambiar el estadístico 
`Dn` en la función `ks.exp.sim` por `Dn = sum(abs( (1:n - 0.5)/n -  pexp(sort(x), rate=rate) ))`.

:::


### Contrastes de permutaciones

Supongamos que a partir de una muestra 
$\left\{ \left( \mathbf{X}_i, Y_i\right): i=1,\ldots, n \right\}$
estamos interesados en contrastar la hipótesis nula de independencia
entre $\mathbf{X}$ e $Y$:
$$H_0: F_{Y \mid \mathbf{X}} = F_Y$$
o equivalentemente que $\mathbf{X}$ no influye en la distribución de $Y$.

En este caso los valores de la respuesta serían intercambiables bajo la hipótesis nula,
por lo que podríamos obtener las remuestras manteniendo fijos los valores^[Nótese que 
no se hace ninguna suposición sobre el tipo de covariables,
podrían ser categóricas, numéricas o una combinación de ambas.]
$\mathbf{X}_i$ y permutando los $Y_i$. Es decir:

1. Generar $Y^{\ast}_i$, con $i=1,\ldots, n$, mediante muestreo
    sin reemplazamiento de $\left\{ Y_i: i=1,\ldots, n \right\}$.

2. Considerar la remuestra bootstrap
   $\left\{ \left( \mathbf{X}_i, Y^{\ast}_i\right): i=1,\ldots, n \right\}$.
   
Se pueden realizar contrastes de este tipo con el paquete `boot` estableciendo 
el parámetro `sim = "permutation"` al llamar a la función `boot()` (el argumento
`i` de la función `statistic` contendrá permutaciones del vector de índices).
Puede ser también de interés el paquete [`coin`](https://cran.r-project.org/web/packages/coin/index.html), 
que implementa muchos contrastes de este tipo.

::: {.example #perm-test-cor name="Inferencia sobre el coeficiente de correlación lineal"}
<br>

Continuando con el ejemplo de las secciones \@ref(boot-unif-multi) y \@ref(icboot-trans) (y de los Ejercicios \@ref(exr:unif-multi), \@ref(exr:unif-multi-sesgo-var) y \@ref(exr:boot-unif-multi-ic)), estamos interesados en hacer inferencia sobre el coeficiente de correlación lineal $\rho$ empleando el coeficiente de correlación muestral $r$ como estimador.
En este caso sin embargo, consideraremos como ejemplo el conjunto de datos `dogs` 
del paquete `boot`, que contiene observaciones sobre el consumo de 
oxígeno cardíaco (`mvo`) y la presión ventricular izquierda (`lvp`) 
de 7 perros domésticos.


```r
library(boot)
data('dogs', package = "boot")
# plot(dogs)
cor(dogs$mvo, dogs$lvp)
```

```
## [1] 0.85369
```

```r
# with(dogs, cor(mvo, lvp))
```

Como ya se comentó, para realizar inferencias sobre $\rho$ podemos emplear la función `cor.test()`:

```r
cor.test(dogs$mvo, dogs$lvp)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  dogs$mvo and dogs$lvp
## t = 3.67, df = 5, p-value = 0.015
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.28180 0.97801
## sample estimates:
##     cor 
## 0.85369
```

```r
# with(dogs, cor.test(mvo, lvp))
```
Esta función realiza el contraste $H_0: \rho = 0$ empleando el estadístico:
$$\frac{r\sqrt{n - 2}}{\sqrt{1 - r^2}} \underset{aprox}{\sim } t_{n-2},$$
<!-- \mathcal{t}_{n-2} error en LaTeX-->
bajo la hipótesis nula de que la verdadera correlación es cero.
Alternativamente se pueden realizar contrastes unilaterales estableciendo
el parámetro `alternative` igual a `"less"` o `"greater"`.
Por ejemplo, para contrastar $H_0: \rho \leq 0$ podríamos emplear:

```r
cor.test(dogs$mvo, dogs$lvp, alternative = "greater")
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  dogs$mvo and dogs$lvp
## t = 3.67, df = 5, p-value = 0.0073
## alternative hypothesis: true correlation is greater than 0
## 95 percent confidence interval:
##  0.41959 1.00000
## sample estimates:
##     cor 
## 0.85369
```

Para realizar el contraste con la función `boot` podríamos
emplear el siguiente código:

```r
library(boot)

statistic <- function(data, i) cor(data$mvo, data$lvp[i])

set.seed(1)
res.boot <- boot(dogs, statistic, R = 1000, sim = "permutation")
# res.boot
```

Posteriormente emplearíamos las réplicas (almacenadas en `res.boot$t`) y el valor
observado del estadístico del contraste (almacenado en `res.boot$t0`)
para aproximar el $p$-valor:

```r
hist(res.boot$t, freq = FALSE, main ="")
abline(v = res.boot$t0, lty = 2)
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/perm-test-cor-1.png" alt="Distribución del estadístico del contraste bajo la hipótesis nula aproximada mediante permutación de las observaciones." width="70%" />
<p class="caption">(\#fig:perm-test-cor)Distribución del estadístico del contraste bajo la hipótesis nula aproximada mediante permutación de las observaciones.</p>
</div>

Por ejemplo, para el contraste unilateral $H_0: \rho \leq 0$ 
(`alternative = "greater"`), obtendríamos:

```r
pval.greater <- mean(res.boot$t >= res.boot$t0)
pval.greater
```

```
## [1] 0.009
```

Mientras que para realizar el contraste bilateral $H_0: \rho = 0$
(`alternative = "two.sided"`), sin asumir que 
la distribución del estadístico de contraste es simétrica:

```r
pval.less <- mean(res.boot$t <= res.boot$t0)
pval <- 2*min(pval.less, pval.greater)
pval
```

```
## [1] 0.018
```

:::

### Contrastes bootstrap semiparamétricos {#contrastes-semiparametricos}

Este tipo de aproximación se emplearía en el caso de que la hipótesis nula 
(o la alternativa) especifique un modelo semiparamétrico, 
con una componente paramétrica y otra no paramétrica.
Típicamente se incluye el error en la componente no paramétrica, y podríamos emplear
el bootstrap residual (también denominado semiparamétrico o basado en modelos) 
descrito en la Sección \@ref(boot-residual).

En esta sección nos centraremos en inferencia sobre modelos de regresión lineales (aunque el procedimiento sería análogo en el caso de modelos más generales), empleando como ejemplo el conjunto de datos `Prestige` del paquete `carData`, considerando como variable respuesta `prestige` (puntuación de ocupaciones obtenidas a partir de una encuesta) y como variables explicativas: `income` (media de ingresos en la ocupación) y `education` (media de los años de educación).


```r
data(Prestige, package = "carData")
# ?Prestige
```

En la mayoría de los casos nos interesa contrastar un **modelo reducido**
frente a un **modelo completo** (que generaliza el modelo reducido).
Por ejemplo, en el caso de modelos lineales (estimados por mínimos cuadrados) 
se dispone del test $F$ para realizar los contrastes de este tipo, 
que emplea el estadístico:
$$F=\frac{n - q}{q - q_0}\frac{RSS_0 - RSS}{RSS},$$
<!-- $$F=\frac{\frac{RSS_0 - RSS}{q - q_0}}{\frac{RSS}{n - q}},$$ -->
siendo $n$ el número de observaciones, $RSS$ y $q$ la suma de cuadrados residual y 
el número de parámetros distintos del modelo completo, 
y $RSS_0$ y $q_0$ los correspondientes al modelo reducido.
Este estadístico sigue una distribución $\mathcal{F}_{q - q_0, n - q}$ 
bajo $H_0$ y las hipótesis habituales del modelo lineal 
($\varepsilon_i$ i.i.d. $\mathcal{N}(0, \sigma^2)$).

El contraste de regresión sería un caso particular. Por ejemplo, 
para contrastar si `income` y `education` influyen linealmente en `prestige` 
podemos emplear el siguiente código:


```r
modelo <- lm(prestige ~ income + education, data = Prestige)
summary(modelo)
```

```
## 
## Call:
## lm(formula = prestige ~ income + education, data = Prestige)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -19.404  -5.331   0.015   4.980  17.689 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -6.847779   3.218977   -2.13    0.036 *  
## income       0.001361   0.000224    6.07  2.4e-08 ***
## education    4.137444   0.348912   11.86  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.81 on 99 degrees of freedom
## Multiple R-squared:  0.798,	Adjusted R-squared:  0.794 
## F-statistic:  196 on 2 and 99 DF,  p-value: <2e-16
```

También podemos obtener el valor observado del estadístico $F$ 
a partir de los resultados del método `summary.lm()`:


```r
res <- summary(modelo)
# names(res)
stat <- res$fstatistic[1]
df <- res$fstatistic[2]
dfr <- res$fstatistic[3]
res$fstatistic
```

```
##  value  numdf  dendf 
## 195.55   2.00  99.00
```

o haciendo los cálculos a mano:


```r
n <- nrow(Prestige)
q <- 3
q0 <- 1
rss0 <- with(Prestige, sum((prestige - mean(prestige))^2))
rss <- sum(residuals(modelo)^2)
inc.mse <- (rss0 - rss)/(q - q0)  # Incremento en varibilidad explicada
msr <-  rss/(n - q)               # Variabilidad residual
inc.mse/msr
```

```
## [1] 195.55
```

Desde el punto de vista de comparación de modelos, el modelo
reducido bajo la hipótesis nula es:


```r
modelo0 <- lm(prestige ~ 1, data = Prestige)
```

y podemos realizar el contraste mediante la función `anova()`


```r
anova(modelo0, modelo)
```

```
## Analysis of Variance Table
## 
## Model 1: prestige ~ 1
## Model 2: prestige ~ income + education
##   Res.Df   RSS Df Sum of Sq   F Pr(>F)    
## 1    101 29895                            
## 2     99  6039  2     23857 196 <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Para aproximar la distribución de este estadístico bajo $H_0$ podríamos adaptar 
el bootstrap semiparamétrico^[En este caso también podríamos emplear un contraste 
de permutaciones.] descrito en la Sección \@ref(boot-residual):


```r
library(boot)

pres.dat <- Prestige
# pres.dat$fit0 <- mean(Prestige$prestige)
# pres.dat$fit0 <- predict(modelo0)
pres.dat$res0 <- with(Prestige, prestige - mean(prestige))
# pres.dat$res0 <- residuals(modelo0)

mod.stat <- function(data, i) {
    data$prestige <- mean(data$prestige) + data$res0[i]
    fit <- lm(prestige ~ income + education, data = data)
    summary(fit)$fstatistic[1]
}

set.seed(1)
boot.mod <- boot(pres.dat, mod.stat, R = 1000)
boot.mod
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = pres.dat, statistic = mod.stat, R = 1000)
## 
## 
## Bootstrap Statistics :
##     original  bias    std. error
## t1*   195.55 -194.49      1.0963
```

```r
hist(boot.mod$t, freq = FALSE, breaks = "FD", main = "")
curve(pf(x, df, dfr, lower.tail = FALSE), lty = 2, add = TRUE)
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/semipar-test-plot-1.png" alt="Distribución del estadístico del contraste (bajo la hipótesis nula) aproximada mediante bootstrap semiparamétrico." width="70%" />
<p class="caption">(\#fig:semipar-test-plot)Distribución del estadístico del contraste (bajo la hipótesis nula) aproximada mediante bootstrap semiparamétrico.</p>
</div>

```r
# pval <- mean(boot.mod$t >= boot.mod$t0)
pval <- mean(boot.mod$t >= stat)
pval
```

```
## [1] 0
```

Procediendo de esta forma sin embargo estaríamos sobreestimando la variabilidad
del error cuando la hipótesis nula es falsa (la variabilidad no explicada por la 
tendencia es asumida por el error), lo que disminuirá la potencia del contraste. 
Para mejorar la potencia, siguiendo la idea propuesta por González-Manteiga 
y Cao (1993), se pueden remuestrear los residuos del modelo completo.
De esta forma reproduciríamos la variabilidad del error de forma consistente 
tanto bajo la hipótesis alternativa como bajo la nula.


```r
old.par <- par(mfrow=c(1,2))
# Variabilidad residual con el modelo reducido
hist(residuals(modelo0), xlim = c(-50, 50), main = "")
# Variabilidad residual con el modelo completo
hist(residuals(modelo), xlim = c(-50, 50), main = "")
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/semipar-var-plot-1.png" alt="Variabilidad residual con el modelo reducido (izquierda) y con el modelo completo (derecha)." width="90%" />
<p class="caption">(\#fig:semipar-var-plot)Variabilidad residual con el modelo reducido (izquierda) y con el modelo completo (derecha).</p>
</div>

```r
par(old.par)
```

Adicionalmente, como se mostró en la Sección \@ref(boot-residual), se puede emplear
la modificación propuesta en Davison y Hinkley (1997, Alg. 6.3, p. 271)
y remuestrear los residuos reescalados y centrados.


```r
pres.dat <- Prestige
# pres.dat$fit0 <- mean(Prestige$prestige)
# pres.dat$fit0 <- predict(modelo0)
# pres.dat$res <- residuals(modelo)
pres.dat$sres <- residuals(modelo)/sqrt(1 - hatvalues(modelo))
pres.dat$sres <- pres.dat$sres - mean(pres.dat$sres)

mod.stat <- function(data, i) {
    # data$prestige <- mean(data$prestige) + data$res[i]
    data$prestige <- mean(data$prestige) + data$sres[i]
    fit <- lm(prestige ~ income + education, data = data)
    summary(fit)$fstatistic[1]
}

set.seed(1)
boot.mod <- boot(pres.dat, mod.stat, R = 1000)
boot.mod
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = pres.dat, statistic = mod.stat, R = 1000)
## 
## 
## Bootstrap Statistics :
##     original  bias    std. error
## t1* 0.011644  1.0297      1.0297
```

En la aproximación del $p$-valor hay que tener en cuenta que al modificar los residuos
`boot.mod$t0` no va a coincidir con el valor observado del estadístico, 
almacenado en `stat` (por tanto habría que ignorar `original` y `bias`
en `Bootstrap Statistics`; 
la función `Boot()` del paquete `car` corrige este problema).


```r
hist(boot.mod$t, freq = FALSE, breaks = "FD", main = "")
curve(pf(x, df, dfr, lower.tail = FALSE), lty = 2, add = TRUE)
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/boot-semipar-test-plot-1.png" alt="Distribución del estadístico del contraste (bajo la hipótesis nula) aproximada mediante bootstrap semiparamétrico." width="70%" />
<p class="caption">(\#fig:boot-semipar-test-plot)Distribución del estadístico del contraste (bajo la hipótesis nula) aproximada mediante bootstrap semiparamétrico.</p>
</div>

```r
pval <- mean(boot.mod$t >= stat)
pval
```

```
## [1] 0
```

En el caso de modelos no lineales (o otros tipos de modelos lineales) puede ser 
complicado aproximar los grados de libertad para el cálculo del estadístico $F$, 
pero si empleamos bootstrap, vamos a obtener los mismos resultados considerando 
como estadístico:
$$\tilde F =\frac{RSS_0 - RSS}{RSS},$$
que se puede interpretar también como una medida del incremento en la variabilidad residual 
al considerar el modelo reducido (ya que únicamente difieren en una constante).
En este caso también se suelen emplear los residuos sin reescalar, ya que también puede ser
difícil encontrar la transformación adecuada.

::: {.exercise #test-semipar-cuadratico}
<br>

Al estudiar el efecto de las variables explicativas en el modelo 
anterior, podríamos pensar que no es adecuado asumir un efecto lineal
de alguna de las variables explicativas. Por ejemplo, si generamos los gráficos
parciales de residuos obtendríamos:


```r
car::crPlots(modelo)
```

<div class="figure" style="text-align: center">
<img src="11-Bootstrap_aplic_files/figure-html/semipar-crplot-1.png" alt="Efecto de las variables explicativas en el modelo (gráficos de residuos)." width="90%" />
<p class="caption">(\#fig:semipar-crplot)Efecto de las variables explicativas en el modelo (gráficos de residuos).</p>
</div>

En este caso podría ser razonable considerar un efecto cuadrático
de la variable `income`^[Para ajustar un modelo polinómico
puede ser recomendable, especialmente si el grado del polinomio es alto, 
emplear la función `poly()` ya que utiliza polinomios ortogonales. 
En el caso cuadrático, al emplear `y ~ x + I(x^2)`
estaremos considerando $1, x, x^2$, mientras que `y ~ poly(x, 2)` considerará 
polinomios de Legendre de la forma $1, x, \frac{1}{2}(3x^2-1)$. 
En este caso concreto, obtendríamos una parametrización equivalente 
empleando `modelo <- lm(prestige ~ poly(income, 2) + education, data = Prestige)`.]


```r
modelo <- lm(prestige ~ income + I(income^2) + education, data = Prestige)
summary(modelo)
```

```
## 
## Call:
## lm(formula = prestige ~ income + I(income^2) + education, data = Prestige)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -15.731  -4.900  -0.057   4.598  18.459 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -1.14e+01   3.27e+00   -3.47  0.00078 ***
## income       3.29e-03   5.67e-04    5.81  7.8e-08 ***
## I(income^2) -7.97e-08   2.17e-08   -3.67  0.00039 ***
## education    3.81e+00   3.41e-01   11.18  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.36 on 98 degrees of freedom
## Multiple R-squared:  0.822,	Adjusted R-squared:  0.817 
## F-statistic:  151 on 3 and 98 DF,  p-value: <2e-16
```

Para comparar el ajuste de este modelo respecto al del anterior, podemos
realizar un contraste empleando la función `anova()`:


```r
modelo0 <- lm(prestige ~ income + education, data = Prestige)
anova(modelo0, modelo)
```

```
## Analysis of Variance Table
## 
## Model 1: prestige ~ income + education
## Model 2: prestige ~ income + I(income^2) + education
##   Res.Df  RSS Df Sum of Sq    F  Pr(>F)    
## 1     99 6039                              
## 2     98 5308  1       731 13.5 0.00039 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Contrastar si el efecto de `income` es lineal mediante bootstrap residual, 
empleando como estadístico el incremento en la variabilidad residual con el 
modelo reducido y remuestreando los residuos del modelo completo (sin reescalar).
Aproximar el nivel crítico del contraste y el valor que tendría que superar el
estadístico para rechazar $H_0$ con un nivel de significación $\alpha = 0.05$.

:::



