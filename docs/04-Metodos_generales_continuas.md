# Simulación de variables continuas {#continuas}




<!-- 
---
title: "Simulación de variables continuas"
author: "Simulación Estadística (UDC)"
date: "Máster en Técnicas Estadísticas"
output: 
  bookdown::html_document2:
    pandoc_args: ["--number-offset", "3,0"]
    toc: yes 
    # mathjax: local            # copia local de MathJax, hay que establecer:
    # self_contained: false     # las dependencias se guardan en ficheros externos 
  bookdown::pdf_document2:
    keep_tex: yes
    toc: yes 
---

bookdown::preview_chapter("04-Metodos_generales_continuas.Rmd")
knitr::purl("04-Metodos_generales_continuas.Rmd", documentation = 2)
knitr::spin("04-Metodos_generales_continuas.R",knit = FALSE)

PENDENTE:
Empregar entorno algoritmo
Empregar exemplos en lugar de exercicios
  Exemplo inferencia bayesiana
-->

En este capítulo se expondrán métodos generales para simular distribuciones continuas: el método de inversión (siguiente sección), los basados en aceptación-rechazo (secciones \@ref(AR) y \@ref(modAR)) y el método de composición (Sección \@ref(composicion)). 
En todos los casos como punto de partida es necesario disponer de un método de generación de números pseudoaleatorios uniformes en $(0,1)$.


## Método de inversión {#inversion}

En general sería el método preferible para la simulación de una variable continua (siempre que se disponga de la función cuantil). 
Está basado en los siguientes resultados:

Si $X$ es una variable aleatoria con función de distribución $F$ continua y estrictamente monótona (invertible), entonces:
$$U = F\left( X \right) \sim \mathcal{U}(0, 1)$$
ya que:
$$G\left( u \right) = P\left( Y \leq u \right) 
= P\left( F(X) \leq u \right) \\
= P\left( X \leq F^{-1}(u) \right) 
= F\left( F^{-1}(u) \right) = u$$

El recíproco también es cierto, si $U \sim \mathcal{U}(0, 1)$ entonces: 
$$F^{-1}\left( U \right) \sim X$$

A partir de este resultado se deduce el siguiente algoritmo genérico para simular una variable continua con función de distribución $F$ invertible:

::: {.conjecture #inversion name="Método de inversión"}
<br>
  
1. Generar $U \sim \mathcal{U}(0, 1)$.

2. Devolver $X = F^{-1}\left( U \right)$.

:::
<!-- \@ref(cnj:inversion) -->

::: {.example #exp-inv name="simulación de una distribución exponencial"}
<br>

La distribución exponencial $\exp \left( \lambda \right)$ de parámetro $\lambda>0$
tiene como función de densidad $f(x) =\lambda e^{-\lambda x}$, si $x\geq 0$,
y como función de distribución:
$$F(x)=\left\{ \begin{array}{ll}
1-e^{-\lambda x} & \text{si } x \ge 0 \\
0 & \text{si } x < 0\\
\end{array} \right.$$

Teniendo en cuenta que:
$$1-e^{-\lambda x}=u \Leftrightarrow x=-\frac{\ln \left( 1-u\right) }{ \lambda }$$
el algoritmo para simular esta variable mediante el método de inversión es:

1. Generar $U \sim \mathcal{U}(0, 1)$.

2. Devolver $X=-\dfrac{\ln \left( 1-U\right) }{\lambda }$.

En el último paso podemos emplear directamente $U$ en lugar de $1-U$, ya que $1 - U \sim \mathcal{U}(0, 1)$.
Esta última expresión para acelerar los cálculos es la que denominaremos *forma simplificada*. 

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/inv-movie-1} 

}

\caption{Ilustración de la simulación de una distribución exponencial por el método de inversión.}(\#fig:inv-movie)
\end{figure}


El código para implementar este algoritmo en R podría ser el siguiente:

```r
tini <- proc.time()

lambda <- 2
nsim <- 10^5
set.seed(1)
u <- runif(nsim)
x <- -log(u)/lambda # -log(1-u)/lambda

tiempo <- proc.time() - tini
tiempo
```

```
##    user  system elapsed 
##    0.03    0.00    0.03
```

```r
hist(x, breaks = "FD", freq = FALSE, 
        main = "", xlim = c(0, 5), ylim = c(0, 2.5))
curve(dexp(x, lambda), lwd = 2, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/exp-inv-plot-1} 

}

\caption{Distribución de los valores generados de una exponencial mediante el método de inversión.}(\#fig:exp-inv-plot)
\end{figure}

Como se observa en la Figura \@ref(fig:exp-inv-plot) se trata de un método exacto (si está bien implementado) y la distribución de los valores generados se aproxima a la distribución teórica como cabría esperar con una muestra de ese tamaño.

:::


### Algunas distribuciones que pueden simularse por el método de inversión

A continuación se incluyen algunas distribuciones que se pueden simular
fácilmente mediante el método de inversión. Se adjunta una forma
simplificada del método que tiene por objeto evitar cálculos
innecesarios (tal y como se hizo en el ejemplo de la exponencial).

| Nombre  | Densidad  | $F(x)$  | $F^{-1}\left( U\right)$  | Forma simplificada |
| ------------- | :-----------: | :-----------: | :-----------: | :-----------: |
| $\exp\left( \lambda\right)$ ($\lambda>0$)  | $\lambda e^{-\lambda x}$, si $x\geq0$ | $1-e^{-\lambda x}$  | $-\dfrac{\ln\left( 1-U\right)  }\lambda$  | $-\dfrac{\ln U}\lambda$ |
| Cauchy  | $\dfrac1{\pi\left( 1+x^{2}\right)}$  | $\dfrac12+\dfrac{\arctan x}\pi$ | $\tan\left( \pi\left( U-\dfrac12\right) \right)$ | $\tan\pi U$ |
| Triangular en $\left( 0,a\right)$ | $\dfrac2a\left( 1-\dfrac xa\right)$, si $0\leq x\leq a$  | $\dfrac2a\left(x-\dfrac{x^{2}}{2a}\right)$  | $a\left( 1-\sqrt{1-U}\right)$  | $a\left( 1-\sqrt{U}\right)$ |
| Pareto ($a,b>0$) | $\dfrac{ab^{a}}{x^{a+1}}$, si $x\geq b$ | $1-\left( \dfrac bx\right)^{a}$ | $\dfrac b{\left( 1-U\right)^{1/a}}$ | $\dfrac b{U^{1/a}}$ |
| Weibull ($\lambda,\alpha>0$) | $\alpha\lambda^{\alpha}x^{\alpha-1}e^{-\left( \lambda x\right)  ^{\alpha}}$, si $x\geq0$ | $1-e^{-\left( \lambda x\right)^{\alpha}}$ | $\dfrac{\left( -\ln\left(1-U\right)  \right)^{1/\alpha}}\lambda$ | $\dfrac{\left( -\ln U\right)^{1/\alpha}}\lambda$ | 


::: {.exercise #ddexp name="distribución doble exponencial"}
<br>

La distribución doble exponencial^[Esta distribución también se puede generar fácilmente simulando una distribución exponencial y asignando un signo positivo o negativo con equiprobabilidad (ver Ejemplo \@ref(exm:dexp-mix)) y función [`simres::rdexp()`](https://rubenfcasal.github.io/simres/reference/ddexp.html) (fichero [*ddexp.R*](R/ddexp.R)).] (o distribución de Laplace) de
parámetro $\lambda$ tiene función de densidad:
$$f(x)  =\frac{\lambda}{2}e^{-\lambda\left\vert x\right\vert
}\text{, }x\in\mathbb{R}$$ 

y función de distribución:

$$F(x)  =\int_{-\infty}^{x}f\left( t\right)  dt=\left\{
\begin{array}{ll}
\frac{1}{2}e^{\lambda x} & \text{si } x<0\\
1-\frac{1}{2}e^{-\lambda x} & \text{si } x\geq0
\end{array}
\ \right.$$

a)  Escribir una función que permita generar, por el método de
    inversión, una muestra de $n$ observaciones de esta distribución.

    
    ```r
    ddexp <- function(x, lambda = 1){
    # Densidad doble exponencial
      lambda*exp(-lambda*abs(x))/2
    }
    
    rdexp <- function(lambda = 1){
    # Simulación por inversión
    # Doble exponencial
      u <- runif(1)
      if (u<0.5) {
        return(log(2*u)/lambda)
      } else {
        return(-log(2*(1-u))/lambda)
      }
    }
    
    rdexpn <- function(n = 1000, lambda = 1) {
    # Simulación n valores de doble exponencial
        x <- numeric(n)
        for(i in 1:n) x[i]<-rdexp(lambda)
        return(x)
    }
    ```


b)  Generar $10^{4}$ valores de la distribución doble exponencial de
    parámetro $\lambda=2$ y obtener el tiempo de CPU que tarda en
    generar la secuencia.

    
    ```r
    set.seed(54321)
    system.time(x <- rdexpn(10^4, 2))
    ```
    
    ```
    ##    user  system elapsed 
    ##    0.05    0.00    0.04
    ```


c)  Representar el histograma y compararlo con la densidad teórica.

    
    ```r
    hist(x, breaks = "FD", freq = FALSE, main="")
    # lines(density(x), col = 'blue')
    curve(ddexp(x, 2), add = TRUE)
    ```
    
    \begin{figure}[!htb]
    
    {\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/ddexp-inv-1} 
    
    }
    
    \caption{Distribución de los valores generados de una doble exponencial mediante el método de inversión.}(\#fig:ddexp-inv)
    \end{figure}
    
    Como se trata de un método exacto de simulación, si está bien implementado, la distribución de los valores generados debería comportarse como una muestra genuina de la distribución objetivo.
    
:::


### Ventajas e inconvenientes

La principal ventaja de este método es que, en general, sería aplicable a cualquier distribución continua (como se muestra en la Sección \@ref(transcuant), se puede extender al caso de que la función de distribución no sea invertible, incluyendo distribuciones discretas).

Uno de los principales problemas es que puede no ser posible encontrar una expresión explícita para $F^{-1}\left( u\right)$ (en ocasiones, como en el caso de la distribución normal, ni siquiera se dispone de una expresión explícita para la función de distribución).
Además, aún disponiendo de una expresión explícita para $F^{-1}\left( u\right)$, su evaluación directa puede requerir de mucho tiempo de computación.

Como alternativa a estos inconvenientes se podrían emplear métodos numéricos para resolver $F(x) - u = 0$ de forma aproximada, aunque habría que resolver numéricamente esta ecuación para cada valor aleatorio que se desea generar.
Otra posibilidad, en principio preferible, sería emplear una aproximación a $F^{-1}\left( u\right)$, dando lugar al *método de inversión aproximada* (como se indicó en la Sección \@ref(oprrng), R emplea por defecto este método para la generación de la distribución normal).
  

### Inversión aproximada

En muchos casos en los que no se puede emplear la expresión exacta de la función
cuantil $F^{-1}\left( u\right)$, se dispone de una aproximación suficientemente
buena que se puede emplear en el algoritmo anterior (se obtendrían simulaciones
con una distribución aproximada a la deseada).

Por ejemplo, para aproximar la función cuantil de la normal estándar, 
Odeh y Evans (1974) consideraron la siguiente función auxiliar^[R emplea una aproximación similar, basada en el algoritmo de Wichura (1988) más preciso, y que está implementado en el fichero fuente [qnorm.c](https://svn.r-project.org/R/trunk/src/nmath/qnorm.c).]:
$$ g\left( v\right)  =\sqrt{-2\ln v}\frac{A\left( \sqrt{-2\ln v}\right)
}{B\left( \sqrt{-2\ln v}\right)  },$$
siendo $A(x) =\sum_{i=0}^{4}a_{i}x^{i}$
y $B(x)  =\sum_{i=0}^{4}b_{i}x^{i}$ con:

$$\begin{array}{ll}
a_{0}=-0.322232431088 &  b_{0}=0.0993484626060 \\
a_{1}=-1 &  b_{1}=0.588581570495 \\
a_{2}=-0.342242088547 & b_{2}=0.531103462366 \\
a_{3}=-0.0204231210245 & b_{3}=0.103537752850 \\
a_{4}=-0.0000453642210148 & b_{4}=0.0038560700634
\end{array}$$

La aproximación consiste en utilizar $g\left( 1-u\right)$ en lugar de
$F^{-1}\left( u\right)$ para los valores de $u\in[10^{-20},\frac12]$
y $-g\left( u\right)$ si $u\in[\frac12,1-10^{-20}]$. Para $u\notin
[10^{-20},1-10^{-20}]$ (que sólo ocurre con una probabilidad de
$2\cdot10^{-20}$) la aproximación no es recomendable.


::: {.conjecture #Odeh-Evans name="de Odeh y Evans"}
<br>

1. Generar $U \sim U(0, 1)$.

2. Si $U<10^{-20}$ ó $U>1-10^{-20}$ entonces volver a 1.

3. Si $U<0.5$ entonces hacer $X=g\left(1-U\right)$ 
   en caso contrario hacer $X=-g\left( U\right)$.

4. Devolver $X$.

:::
<!-- \@ref(cnj:Odeh-Evans) -->


En manuales de funciones matemáticas, como [Abramowitz y Stegun (1964)](https://www.math.ubc.ca/~cbm/aands/frameindex.htm), 
se tienen aproximaciones de la función cuantil de las principales distribuciones
(por ejemplo en la página [993](https://www.math.ubc.ca/~cbm/aands/page_933.htm)
las correspondientes a la normal estándar). 


## Método de aceptación rechazo {#AR}

Se trata de un método universal alternativo al de inversión para
el caso de que no se pueda emplear la función cuantil, 
pero se dispone de una expresión (preferiblemente sencilla) para la
función de densidad objetivo $f\left( x \right)$. 

La idea es simular una variable aleatoria bidimensional $\left( X, Y\right)$ con distribución uniforme en el hipografo de $f$ (el conjunto de puntos del plano comprendidos entre el eje OX y $f$):
$$A_{f}=\left\{ \left( x,y\right) \in \mathbb{R}^{2}:0\leq y\leq f(x) \right\}.$$
De esta forma la primera componente tendrá la distribución deseada (Figura \@ref(fig:var-gen)):

\begin{figure}[!htb]

{\centering \includegraphics[width=0.8\linewidth]{04-Metodos_generales_continuas_files/figure-latex/var-gen-1} 

}

\caption{Puntos con distribución uniforme en el hipografo de una función de densidad.}(\#fig:var-gen)
\end{figure}

$$ P\left( a<X<b\right) = \frac{\text{Area de }\left\{ \left( x,y\right) \in 
\mathbb{R}^{2}:a<x<b;~0\leq y\leq f(x) \right\} }{\text{Area de }
A_{f}} \\
= \int_{a}^{b}f(x) dx $$

El resultado anterior es también válido para una cuasi-densidad $f^{\ast}$ (no depende de la constante normalizadora): 

* Si $\left( X,Y\right) \sim \mathcal{U}\left(A_{f^{\ast}}\right)$ entonces^[Emplearemos también $X\sim f$ para indicar que $X$ es una variable aleatoria con función de densidad $f$.] $X\sim f$.

Para simular una distribución uniforme en el hipografo $A_{f}$ (o en $A_{f^{\ast}}$), lo que se hace es utilizar una variable aleatoria auxiliar $T$ con función de densidad $g$, similar a $f$ y fácil de simular, y una constante $c > 0$ verificando:
$$f(x) \leq c\cdot g(x) \text{, }\forall x\in \mathbb{R}.$$
Podemos generar valores en $A_{cg} \supset A_{f}$ empleando el resultado:

* Si $T$ es una variable aleatoria con función de densidad $g$ 
  y $U \sim \mathcal{U}(0, 1)$ entonces
  $$\left( T,c\cdot U\cdot g(x) \right) \sim \mathcal{U}\left(
  A_{cg}\right)$$
  siendo
  $A_{cg}=\left\{ \left( x, y \right) \in \mathbb{R}^{2} : 0 \leq y \leq
  cg\left( x \right) \right\}$.

Teniendo en cuenta además que:

* Si $\left( T,Y\right) \sim \mathcal{U}\left( A\right)$ y 
  $B \subset A\Rightarrow \left. \left( T,Y\right) \right\vert _{B}
  \sim \mathcal{U}\left(B\right)$.
  
Entonces, si $\left( T, Y \right)$ sigue una distribución uniforme en $A_{cg}$, aceptando los valores de $\left( T, Y\right)$ que pertenezcan a $A_{f}$ (o a $A_{f^{\ast}}$) se obtendrán generaciones con distribución uniforme sobre $A_{f}$ (o $A_{f^{\ast}}$) y la densidad de la primera componente será $f$.


### Algoritmo

Supongamos que $f$ es la densidad objetivo y $g$ es una densidad
auxiliar (fácil de simular y similar a $f$), de forma que
existe una constante $c>0$ tal que:
$$f(x) \leq c\cdot g(x) 
\text{, }\forall x\in \mathbb{R},$$
(de donde se deduce que el soporte de $g$ debe contener el de $f$).

::: {.conjecture #aceptacion-rechazo name="Método de aceptación-rechazo; Von Neuman 1951"}
<br>

1.  Generar $U \sim \mathcal{U}(0, 1)$.

2.  Generar $T\sim g$.

3.  Si $c\cdot U\cdot g\left( T\right) \leq f\left( T\right)$
    devolver $X=T$,

    en caso contrario volver al paso 1.

:::
<!-- \@ref(cnj:aceptacion-rechazo) -->


### Densidades acotadas en un intervalo cerrado

Sea $f$ una función de densidad cualquiera con soporte en un intervalo cerrado $[a,b]$ (es decir, $\{x : f(x) > 0\}=[a,b]$) de tal forma que existe una constante $M>0$ tal que $f(x) \leq M$ $\forall x$ (es decir, $f$ es acotada superiormente). 
En este caso puede tomarse como densidad auxiliar $g$, la de una $\mathcal{U}(a,b)$. 
En efecto, tomando $c = M\left( b-a\right)$ y teniendo en cuenta que
$$g(x) = \left\{
\begin{array}{ll}\frac{1}{b-a} & \text{si } x \in [a,b]\\
0 & \text{en caso contrario}
\end{array} \right.$$ 
se tiene que $f(x) \leq M = \frac{c}{b-a}=c \cdot g(x)$, 
$\forall x \in [a,b]$. 
Así pues, el algoritmo quedaría como sigue:


1.  Generar $U,V\sim \mathcal{U}(0, 1)$.

2.  Hacer $T = a + \left( b-a \right) V$.

3.  Si $M \cdot U\leq f\left( T \right)$
    devolver $X = T$,

    en caso contrario volver al paso 1.


::: {.remark}
No confundir $M$ con $c = M \left( b - a \right)$.
:::

<!-- 
Pendiente: cambiar beta por triangular en [0, 2] 
Emplear código dtri?
-->

::: {.example #dbeta-dunif-ar name="simulación de distribución beta a partir de la uniforme"}
<br>

Para simular una variable con función de densidad $\mathcal{Beta}(\alpha, \beta)$:
$$f(x)=\frac{\Gamma (\alpha + \beta )}{\Gamma (\alpha )\Gamma (\beta )}
x^{\alpha -1}(1-x)^{\beta -1}\text{ si }0 \leq x \leq 1,$$
(siguiendo la notación de la función `dbeta(x, shape1, shape2)` de `R`), podemos considerar como distribución auxiliar una $\mathcal{U}(0,1)$,
con $g(x) = 1$ si $0 \leq x \leq 1$.


:::

Esta distribución está acotada y es unimodal, si $\alpha$ y $\beta$ son mayores^[Si $\alpha$ o $\beta$ son iguales a 1 puede simularse fácilmente por el método de inversión y si alguno es menor que 1 esta densidad no está acotada.] que 1, y su moda es $\frac{\alpha - 1} {\alpha + \beta - 2}$, por lo que: 
$$c = M = \max_{0 \leq x \leq 1}f(x) = f\left( \frac{\alpha - 1} {\alpha + \beta - 2} \right).$$

Por ejemplo, considerando $\alpha = 2$ y $\beta = 4$, si comparamos la densidad objetivo con la auxiliar reescalada (Figura \@ref(fig:dbeta-dunif-plot)), confirmamos que esta última está por encima (y llegan a tocarse, por lo que se está empleando la cota óptima; ver siguiente sección).


```r
# densidad objetivo: dbeta
# densidad auxiliar: dunif
s1 <- 2
s2 <- 4
curve(dbeta(x, s1, s2), -0.1, 1.1, lwd = 2)
m <- dbeta((s1 - 1)/(s1 + s2 - 2), s1, s2)
# abline(h = m, lty = 2)
segments(0, m, 1, m, lty = 2, lwd = 2)
abline(v = 0, lty = 3)
abline(v = 1, lty = 3)
abline(h = 0, lty = 3)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/dbeta-dunif-plot-1} 

}

\caption{Densidad objetivo beta (línea continua) y densidad auxiliar unifome reescalada (línea discontinua).}(\#fig:dbeta-dunif-plot)
\end{figure}

El siguiente código implementa el método de aceptación-rechazo para simular valores de la densidad objetivo (se incluye una variable "global" `ngen` para contar el número de generaciones de la distribución auxiliar):


```r
ngen <- 0 

rbeta2 <- function(s1 = 2, s2 = 2) {
  # Simulación por aceptación-rechazo
  # Beta a partir de uniforme
  m <- dbeta((s1 - 1)/(s1 + s2 - 2), s1, s2)
  while (TRUE) {
    u <- runif(1)
    x <- runif(1)
    ngen <<- ngen+1
    if (m*u <= dbeta(x, s1, s2)) return(x)
  }
}

rbeta2n <- function(n = 1000, s1 = 2, s2 = 2) {
  # Simulación n valores Beta(s1, s2)
  x <- numeric(n)
  for(i in 1:n) x[i]<-rbeta2(s1, s2)
  return(x)
}
```

Empleando estas funciones podemos generar una muestra de $10^3$ observaciones de una $\mathcal{Beta}(2, 4)$ (calculando de paso el tiempo de CPU):


```r
set.seed(1)
nsim <- 1000
ngen <- 0
system.time(x <- rbeta2n(nsim, s1, s2))
```

```
##    user  system elapsed 
##    0.03    0.00    0.04
```

Para analizar la eficiencia podemos emplear el número de generaciones de la distribución auxiliar (siguiente sección):


```r
{cat("Número de generaciones = ", ngen)
cat("\nNúmero medio de generaciones = ", ngen/nsim)
cat("\nProporción de rechazos = ", 1-nsim/ngen, "\n")}
```

```
## Número de generaciones =  2121
## Número medio de generaciones =  2.121
## Proporción de rechazos =  0.5285243
```

Finalmente podemos representar la distribución de los valores generados y compararla con la densidad teórica:


```r
hist(x, breaks = "FD", freq = FALSE, main = "")
curve(dbeta(x, s1, s2), col = 2, lwd = 2, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/dbeta-ar-1} 

}

\caption{Distribución de los valores generados mediante el método de aceptación-rechazo.}(\#fig:dbeta-ar)
\end{figure}

Al ser un método exacto de simulación (si está bien implementado), la distribución de los valores generados debería comportarse como una muestra genuina de la distribución objetivo.


::: {.exercise #dacotada-ar}
<br>

Dar un algoritmo para simular la función de densidad dada por 
$f\left(x\right) = \frac{1}{16} \left( 3x^{2}+2x+2 \right)$ si
$0 \le x \le 2$, cero en otro caso. Estudiar su eficiencia.

:::

### Eficiencia del algoritmo

Como medida de la eficiencia del algoritmo de aceptación-rechazo podríamos considerar el número de iteraciones del algoritmo, es decir, el número de generaciones de la densidad auxiliar y de comparaciones para aceptar un valor de la densidad objetivo.
Este número $N$ es aleatorio y sigue una distribución geométrica (número de pruebas necesarias hasta obtener el primer éxito) con parámetro $p$ (probabilidad de éxito) la probabilidad de aceptación en el paso 3:
$$p = \frac{\text{area}\left(A_{f}\right)}{\text{area}\left( A_{cg}\right)} = \frac{1}{c}.$$
Por tanto:
$$E\left( N \right) = \frac1p = c$$ 
es el número medio de iteraciones del algoritmo (el número medio de pares de variables $\left( T,U\right)$ que se necesitan generar, y de comparaciones, para obtener una simulación de la densidad objetivo).

Es obvio, por tanto, que cuanto más cercano a 1 sea el valor de $c$ más eficiente será el algoritmo (el caso de $c=1$ se correspondería con $g=f$ y no tendría sentido emplear este método).
Una vez fijada la densidad $g$, el valor óptimo será:
$$c_{\text{opt}}=\max_{\{x : g(x) > 0\}} \frac{f(x)}{g(x)}.$$

<br>

::: {.remark}
Hay que tener en cuenta que la cota óptima es el número medio de iteraciones $c$ solo si conocemos las constantes normalizadoras. 
Si solo se conoce la cuasidensidad $f^{\ast}$ de la distribución objetivo (o de la auxiliar), la correspondiente cota óptima:
$$\tilde{c} = \max_{\{x : g(x) > 0\}} \frac{f^{\ast}(x)}{g(x)}$$ 
asumirá la constante desconocida, aunque siempre podemos aproximar por simulación el verdadero valor de $c$ y a partir de él la constante normalizadora (ver Ejercicio \@ref(exr:post-pri-ar)).
Basta con tener en cuenta que, si $f(x) = f^{\ast}(x)/k$:
$$\frac{1}{c} = \frac{\text{area}\left(A_{f^{\ast}}\right)}{\text{area}\left( A_{\tilde{c}g}\right)} = \frac{k}{\tilde{c}},$$
y por tanto $k= \tilde{c}/c$.
:::


::: {.example #dnorm-ddexp-ar name="simulación de la normal a partir de la doble exponencial"}
<br>

Se trata de simular la distribución normal estándar, con función de densidad:
$$f(x)  =\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}} \text{, } x\in\mathbb{R}\text{, }$$
empleando el método de aceptación-rechazo considerando como distribución auxiliar una doble exponencial con $\lambda=1$ (o distribución de Laplace):
$$g(x)  =\frac{1}{2}e^{-\left| x \right|} \text{, } x\in\mathbb{R}.$$
Esta distribución se utilizó en el Ejercicio \@ref(exr:ddexp), donde se definió la densidad auxiliar `ddexp(x, lambda)` y la función `rdexp(lambda)` para generar un valor aleatorio de esta distribución.

En este caso el soporte de ambas densidades es la recta real y el valor óptimo para $c$ es:
$$c_{\text{opt}} = \max_{x\in\mathbb{R}}\frac{f(x)}{g(x) } = \max_{x\in\mathbb{R}} \frac{\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}}{\frac{1}{2}e^{-\left| x\right|  }} = \sqrt{\frac{2}{\pi}}\max_{x\in\mathbb{R}}e^{\varphi(x)} = \sqrt{\frac{2}{\pi}}e^{\max_{x\in\mathbb{R}}\varphi(x)},$$
donde $\varphi(x) = -\frac{x^{2}}{2}+\left| x \right|$.
Dado que esta función es simétrica, continua en toda la recta real y diferenciable tantas veces como se desee salvo en $x=0$, bastará encontrar su máximo absoluto en el intervalo $[0,\infty]$:
$$\begin{aligned}
x  & >0\Rightarrow\varphi^{\prime}(x)  =-x+1,\varphi
^{\prime\prime}(x)  =-1;\\
\{x  & >0,\varphi^{\prime}(x)  =0\}\Leftrightarrow x=1.
\end{aligned}$$ 
Por tanto, como $\varphi^{\prime\prime}(1) <0$, $\varphi$ alcanza un máximo relativo en $x=1$ y otro de idéntico valor en $x=-1$. 
Resulta fácil demostrar que ambos son máximos absolutos
(por los intervalos de crecimiento y decrecimiento de la función).
Como consecuencia:
$$c_{\text{opt}} = \sqrt{\frac{2}{\pi}}e^{\varphi(1)}
=\sqrt{\frac{2}{\pi}}e^{1/2} =\sqrt{\frac{2e}{\pi}} \approx 1.3155.$$ 

Si comparamos la densidad objetivo con la auxiliar reescalada con los parámetros óptimos (Figura \@ref(fig:dnorm-ddexp-plot)), vemos que esta última está por encima, como debería ocurrir, pero llegan a tocarse (lo que validaría el cálculo para la obtención de la cota óptima).


```r
# densidad objetivo: dnorm
# densidad auxiliar: ddexp
c.opt <- sqrt(2*exp(1)/pi)
lambda.opt <- 1
curve(c.opt * ddexp(x), xlim = c(-4, 4), lty = 2)
curve(dnorm, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/dnorm-ddexp-plot-1} 

}

\caption{Densidad objetivo (normal estándar, línea continua) y densidad auxiliar (doble exponencial, línea discontinua) reescalada.}(\#fig:dnorm-ddexp-plot)
\end{figure}

Alternativamente, en lugar de obtener la cota óptima de modo analítico, podríamos aproximarla numéricamente:


```r
# NOTA: Cuidado con los límites
# optimize(f = function(x) dnorm(x)/ddexp(x), maximum = TRUE, interval = c(-0.5,0.5))
optimize(f = function(x) dnorm(x)/ddexp(x), maximum = TRUE, interval = c(0, 2))
```

```
## $maximum
## [1] 1
## 
## $objective
## [1] 1.315489
```

Vemos que la aproximación numérica coincide con el valor óptimo real $c_{\text{opt}} \approx$  1.3154892 (que se alcanza en $x = \pm 1$).

Para establecer la condición de aceptación o rechazo es recomendable emplear 
$c\cdot U\cdot g\left( T\right) \leq f\left( T\right)$.
Aunque, en este caso concreto, se puede tener en cuenta que:
$$c\cdot U\cdot\frac{g\left( T\right)  }{f\left( T\right)  }=\sqrt{\frac
{2e}{\pi}}U\sqrt{\frac{\pi}{2}}\exp\left( \frac{T^{2}}{2}-\left\vert
T\right\vert \right)  =U\cdot\exp\left( \frac{T^{2}}{2}-\left\vert
T\right\vert +\frac{1}{2}\right).$$

Teniendo en cuenta los resultados anteriores, podríamos emplear el siguiente código para generar los valores de la densidad objetivo:


```r
ngen <- 0

rnormAR <- function() {
# Simulación por aceptación-rechazo
# Normal estandar a partir de doble exponencial
  c.opt <- sqrt(2*exp(1)/pi)
  lambda.opt <- 1
  while (TRUE) {
    u <- runif(1)
    x <- rdexp(lambda.opt) # rdexpn(1, lambda.opt)
    ngen <<- ngen + 1 # Comentar esta línea para uso normal
    # if (u*exp((x^2+1)*0.5-abs(x)) <= 1) return(x)
    if (c.opt * u * ddexp(x, lambda.opt) <= dnorm(x)) return(x)
  }
}

rnormARn <- function(n = 1000) {
# Simulación n valores N(0,1)
    x <- numeric(n)
    for(i in 1:n) x[i] <- rnormAR()
    return(x)
}
```


Generamos una muestra de $10^4$ observaciones:


```r
set.seed(1)
nsim <- 10^4
ngen <- 0
system.time(x <- rnormARn(nsim))
```

```
##    user  system elapsed 
##    0.18    0.00    0.20
```

Evaluamos la eficiencia:


```r
{cat("Número de generaciones = ", ngen)
cat("\nNúmero medio de generaciones = ", ngen/nsim)
cat("\nProporción de rechazos = ", 1-nsim/ngen, "\n")}
```

```
## Número de generaciones =  13178
## Número medio de generaciones =  1.3178
## Proporción de rechazos =  0.2411595
```

Estos valores serían aproximaciones por simulación de los correspondientes valores teóricos (valor medio $c \approx 1.3155$ y probabilidad de rechazo $1 - p = 1 - 1/c \approx 0.23983$). 
A partir de ellos podríamos decir que el algoritmo es bastante eficiente.

Finalmente comparamos la distribución de los valores generados con la densidad teórica:
    

```r
hist(x, breaks = "FD", freq = FALSE, main = "")
curve(dnorm, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/dnorm-ar-1} 

}

\caption{Distribución de los valores generados mediante el método de aceptación-rechazo.}(\#fig:dnorm-ar)
\end{figure}

Podemos observar que la distribución de los valores generados es la que cabría esperar de una muestra de tamaño `nsim` de la distribución objetivo (lo que nos ayudaría a confirmar que el algoritmo está bien implementado, al ser un método exacto de simulación).

:::


### Elección de la densidad auxiliar

El principal problema con este método es encontrar una densidad auxiliar $g$ de forma que $c_{\text{opt}}$ sea próximo a 1.
Una solución intermedia consiste en seleccionar una familia paramétrica de densidades $\{g_{\theta} : \theta \in \Theta\}$ entre las que haya alguna que se parezca bastante a $f$, 
encontrar el valor de $c$ óptimo para cada densidad de esa familia:
$$c_{\theta}=\max_{x}\frac{f(x)  }{g_{\theta}(x) }$$ 
y, finalmente, elegir el mejor valor $\theta_{0}$ del parámetro, en el sentido de ofrecer el menor posible $c_{\theta}$:
$$c_{\theta_{0}}=\min_{\theta\in\Theta}\max_{x}\frac{f(x) }{g_{\theta}(x)}.$$


::: {.example #dnorm-ddexp-arb name="simulación de la normal mediante la doble exponencial continuación"}
<br>

Continuando con el Ejemplo \@ref(exm:dnorm-ddexp-ar) anterior sobre la simulación de una normal estándar mediante el método de aceptación-rechazo, en lugar de fijar la densidad auxiliar a una doble exponencial con $\lambda=1$, consideraremos el caso general de $\lambda>0$:
$$g_{\lambda}(x)  = \frac{\lambda}{2}e^{-\lambda\left| x \right|} \text{, } x\in\mathbb{R}.$$

Si pretendemos encontrar el mejor valor de $\lambda$, en términos de eficiencia del algoritmo, debemos calcular:
$$c_{\lambda_{0}} = \min_{\lambda>0}\max_{x\in\mathbb{R}}\frac{f(x)}{g_{\lambda(x)}}
=\min_{\lambda>0}\max_{x\in\mathbb{R}}\frac{\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}}{\frac{\lambda}{2}e^{-\lambda \left| x \right| }}.$$
De forma totalmente análoga a la vista para el caso $\lambda=1$, se tiene que:
$$c_{\lambda}=\max_{x\in\mathbb{R}}\frac{\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}}{\frac{\lambda}{2}e^{-\lambda \left| x \right|  }}=\frac{1}{\lambda}\sqrt{\frac{2}{\pi}}\max_{x\in\mathbb{R}}e^{\varphi
_{\lambda(x)}}=\frac{1}{\lambda}\sqrt{\frac{2}{\pi}}e^{\max_{x\in\mathbb{R}}\varphi_{\lambda(x)} },$$
donde $\varphi_{\lambda}(x)  =-\frac{x^{2}}{2}+\lambda\left| x \right|$.
De forma totalmente similar también puede probarse que $\varphi_{\lambda}$ alcanza su máximo absoluto en los puntos $x=\pm\lambda$, siendo dicho valor máximo
$\varphi_{\lambda}\left( \pm\lambda \right) = \frac{\lambda^{2}}{2}$.
Como consecuencia:
$$c_{\lambda}=\frac{1}{\lambda}\sqrt{\frac{2}{\pi}}e^{\varphi_{\lambda}\left( \pm\lambda \right)} =\frac{e^{\frac{\lambda^{2}}{2}}}{\lambda}\sqrt{\frac{2}{\pi}}.$$

Finalmente debemos encontrar $\lambda_{0}$ tal que $c_{\lambda_{0}}=\min_{\lambda>0}c_{\lambda}$.
Como: 
$$\frac{\partial c_{\lambda}}{\partial \lambda} =\sqrt{\frac{2}{\pi}}\frac{e^{\frac{\lambda^{2}}{2}}\left( \lambda^{2}-1\right)  }{\lambda^{2}},$$
entonces $\frac{\partial c_{\lambda}}{\partial \lambda} = 0\Leftrightarrow\lambda=1$, ya que $\lambda>0$.
Además:
$$\left. \frac{\partial ^{2}c_{\lambda}}{\partial \lambda^{2}}\right|_{\lambda=1} =\left.\sqrt{\frac
{2}{\pi}}\frac{e^{\frac{\lambda^{2}}{2}}\left(  \lambda^{5}-\lambda
^{3}+2\lambda\right)  }{\lambda^{4}}\right|_{\lambda=1}=2\sqrt{\frac{2e}{\pi}}>0,$$
luego en $\lambda=1$ se alcanza el mínimo.


```r
curve(exp(x^2/2)/x*sqrt(2/pi), 0.1, 2.5,
       xlab = expression(lambda), ylab = expression(c[lambda]))
abline(v = 1, lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/cop-lambda-1} 

}

\caption{Representación de la cota óptima dependiedo del valor del parámetro.}(\#fig:cop-lambda)
\end{figure}

De esto se deduce que la mejor densidad auxiliar doble exponencial es la correspondiente a $\lambda=1$. 
Por tanto el algoritmo más eficiente, con esta familia de densidades auxiliares, es el expuesto en el Ejemplo \@ref(exm:dnorm-ddexp-ar).

Alternativamente también podríamos aproximar simultáneamente el parámetro óptimo y la cota óptima de la densidad auxiliar numéricamente:


```r
# Obtención de valores c y lambda óptimos aproximados
fopt <- function(lambda) {
  # Obtiene c fijado lambda
  optimize(f = function(x) dnorm(x)/ddexp(x,lambda),
           maximum = TRUE, interval = c(0, 2))$objective
}

# Encontrar lambda que minimiza
res <- optimize(fopt, interval = c(0.5, 2))
lambda.opt2 <- res$minimum
c.opt2 <- res$objective
lambda.opt2 
```

```
## [1] 0.9999987
```

```r
c.opt2
```

```
## [1] 1.315489
```

:::


### Ejemplo: inferencia bayesiana {#bayes-ar}

El algoritmo de aceptación-rechazo se emplea habitualmente en inferencia bayesiana.
Denotando por:

* $f(x|\theta)$ la densidad muestral.

* $\pi(\theta)$ la densidad a priori.

* $\mathbf{x}=(x_1,...,x_n)^{\top}$ la muestra observada.

El objetivo sería simular la distribución a posteriori de $\theta$:
$$\pi(\theta|\mathbf{x}) = \frac{L(\mathbf{x}|\theta)\pi(\theta)}{\int L(\mathbf{x}|\theta)\pi(\theta)d\theta},$$
siendo $L(\mathbf{x}|\theta)$ la función de verosimilitud ($L(\mathbf{x}|\theta) = \prod\limits_{i=1}^{n}f(x_i|\theta)$ suponiendo i.i.d.). 
Es decir:
$$\pi(\theta | \mathbf{x})\propto L(\mathbf{x}| \theta)\pi(\theta).$$

Como esta distribución cambia al variar la muestra observada, puede resultar difícil encontrar una densidad auxiliar adecuada para simular valores de la densidad a posteriori $\pi(\theta|\mathbf{x})$.
Por ejemplo, podríamos emplear la densidad a priori $\pi(\theta)$ como densidad auxiliar.
Teniendo en cuenta que: 

* $\pi(\theta |\mathbf{x})/\pi(\theta)\propto L(\mathbf{x}|\theta)$

* $L(\mathbf{x}|\theta)\leq \tilde{c}=L(\mathbf{x}|\hat{\theta})$ siendo
  $\hat{\theta}$ el estimador máximo verosímil de $\theta$.

El algoritmo sería el siguiente:

1. Generar $U \sim \mathcal{U}(0, 1)$.

2. Generar $\tilde{\theta}\sim \pi(\theta)$.

3.  Si $L(\mathbf{x}|\hat{\theta})\cdot U \leq 
    L(\mathbf{x}|\tilde{\theta})$ devolver $\tilde{\theta}$,

    en caso contrario volver al paso 1.
    
Aunque, como se muestra en el siguiente ejercicio, esta elección de densidad auxiliar puede ser muy poco adecuada, siendo preferible en la práctica emplear un método adaptativo que construya la densidad auxiliar de forma automática (Sección \@ref(ars)).  

::: {.exercise #post-pri-ar name="Simulación de la distribución a posteriori a partir de la distribución a priori"}
<br>


<!-- 
Basado en el Ejercicio 2.23 de Robert y Casella, 2009 
Cambiar a ejemplo y apartado b) a ejercicio $n = 100$, $\theta_{0}=3$
Comentar resultados (intervalos de probabilidad a posteriori) y eficiencia del método
-->


Para la estimación Bayes de la media de una normal se suele utilizar
como distribución a priori una Cauchy.

a)  Generar una muestra i.i.d. $X_{i}\sim N(\theta_{0},1)$ de tamaño
    $n=10$ con $\theta_{0}=1$. Utilizar una $Cauchy(0,1)$
    (`rcauchy()`) como distribución a priori y como densidad auxiliar
    para simular por aceptación-rechazo una muestra de la densidad a
    posteriori (emplear `dnorm()` para construir la verosimilitud).
    Obtener el intervalo de probabilidad/credibilidad al 95%.


    
    ```r
    mu0 <- 1
    n <- 10
    nsim <- 10^4
    set.seed(54321)
    x <- rnorm(n, mean = mu0)
    
    # Función de verosimilitud
    # lik1 <- function(mu) prod(dnorm(x, mean = mu)) # escalar
    lik <- Vectorize(function(mu) prod(dnorm(x, mean = mu))) # vectorial
    
    # Cota óptima
    # Estimación por máxima verosimilitud
    emv <- optimize(f = lik, int = range(x), maximum = TRUE)
    emv
    ```
    
    ```
    ## $maximum
    ## [1] 0.7353805
    ## 
    ## $objective
    ## [1] 3.303574e-08
    ```
    
    ```r
    c <- emv$objective
    ```

    En este caso concreto, ya sabríamos que el estimador máximo verosímil es la media muestral:

    
    ```r
    mean(x)
    ```
    
    ```
    ## [1] 0.7353958
    ```
    
    y por tanto:
    
    
    ```r
    c <- lik(mean(x))
    c   
    ```
    
    ```
    ## [1] 3.303574e-08
    ```
    
    ```r
    # f.cuasi <- function(mu) sapply(mu, lik1)*dcauchy(mu)
    f.cuasi <- function(mu) lik(mu)*dcauchy(mu)    
    curve(c * dcauchy(x), xlim = c(-4, 4), ylim = c(0, c/pi), lty = 2,
          xlab = "mu", ylab = "cuasidensidad")
    curve(f.cuasi, add = TRUE)
    ```
    
    \begin{figure}[!htb]
    
    {\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/post-pri-plot-1} 
    
    }
    
    \caption{Comparación de la cuasidensidad a posteriori (línea contínua) con la densidad a priori reescalada (línea discontinua).}(\#fig:post-pri-plot)
    \end{figure}
    
    Por ejemplo, podríamos emplear el siguiente código para generar simulaciones de la distribución a posteriori mediante aceptación-rechazo a partir de la distribución de Cauchy:
    
    
    ```r
    ngen <- nsim
    mu <- rcauchy(nsim)
    ind <- c*runif(nsim) > lik(mu) # TRUE si no verifica condición
    # Volver a generar si no verifica condición
    while (sum(ind)>0){
      le <- sum(ind)
      ngen <- ngen + le
      mu[ind] <- rcauchy(le)
      ind[ind] <- c*runif(le) > lik(mu[ind]) # TRUE si no verifica condición
    }
    
    { # Número generaciones
      cat("Número de generaciones = ", ngen)
      cat("\nNúmero medio de generaciones = ", ngen/nsim)
      cat("\nProporción de rechazos = ", 1-nsim/ngen,"\n")
    }
    ```
    
    ```
    ## Número de generaciones =  59422
    ## Número medio de generaciones =  5.9422
    ## Proporción de rechazos =  0.8317122
    ```
    
    A partir de la aproximación del número medio de generaciones podemos aproximar la constante normalizadora:

    
    ```r
    cte <- c*nsim/ngen
    # integrate(f.cuasi, -Inf, Inf)
    f.aprox <- function(mu) f.cuasi(mu)/cte
    ```

    Finalmente, a partir de los valores generados podemos aproximar el intervalo de probabilidad al 95% (intervalo de credibilidad bayesiano):

    
    ```r
    q <- quantile(mu, c(0.025, 0.975))
    q
    ```
    
    ```
    ##       2.5%      97.5% 
    ## 0.05001092 1.26026227
    ```
    
    ```r
    # Representar estimador e IC Bayes
    hist(mu, freq=FALSE, breaks = "FD", main="")
    # abline(v = mean(x), lty = 3) # Estimación frecuentista
    abline(v = mean(mu), lty = 2, lwd = 2)  # Estimación Bayesiana
    abline(v = q, lty = 2)
    curve(f.aprox, col = "blue", add = TRUE)
    ```
    
    \begin{figure}[!htb]
    
    {\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/ic-bayes-1} 
    
    }
    
    \caption{Distribución de los valores generados y aproximación del intervalo de credibilidad.}(\#fig:ic-bayes)
    \end{figure}


b)  Repetir el apartado anterior con $n=100$.

:::



## Modificaciones del método de aceptación-rechazo {#modAR}

En el tiempo de computación del algoritmo de aceptación-rechazo influye:

* La proporción de aceptación (debería ser grande).

* La dificultad de simular con la densidad auxiliar.

* El tiempo necesario para hacer la comparación en el paso 3.

En ciertos casos el tiempo de computación necesario para evaluar $f(x)$ puede ser alto.
Para evitar evaluaciones de la densidad se puede emplear una función "squeeze" que aproxime la densidad por abajo (una envolvente inferior):
$$s(x)\leq f(x) \text{, }\forall x\in \mathbb{R}.$$

::: {.conjecture #marsaglia name="Marsaglia 1977"}
<br>

1.  Generar $U \sim \mathcal{U}(0, 1)$ y $T\sim g$.

2.  Si $c\cdot U\cdot g\left( T\right) \leq s\left( T\right)$ devolver $X=T$,

    en caso contrario

    2.a.  si $c\cdot U\cdot g\left( T\right) \leq f\left( T\right)$ 
          devolver $X=T$,
    
    2.b.  en caso contrario volver al paso 1.

:::

\begin{figure}[!htb]

{\centering \includegraphics[width=0.8\linewidth]{images/squeeze} 

}

\caption{Ilustración del algoritmo de aceptación-rechazo con envolvente inferior (función "squeeze").}(\#fig:squeeze)
\end{figure}

Cuanto mayor sea el área bajo $s(x)$ (más próxima a 1)
más efectivo será el algoritmo.

Se han desarrollado métodos generales para la construcción de las
funciones $g$ y $s$ de forma automática
(cada vez que se evalúa la densidad se mejoran las aproximaciones).
Estos métodos se basan principalmente en que una transformación de
la densidad objetivo es cóncava o convexa.

### Muestreo por rechazo adaptativo (ARS) {#ars}

Supongamos que $f$ es una cuasi-densidad log-cóncava 
(i.e. $\frac{\partial ^{2}}{\partial x^{2}}\log f(x) <0,
~\forall x$).

Sea $S_n=\left\{ x_{i}:i=0,\cdots ,n+1\right\}$ con
$f(x_{i})$ conocidos.

Denotamos por $L_{i,i+1}(x)$ la recta pasando por $\left( x_{i},\log
f(x_{i})\right)$ y $\left( x_{i+1},\log f(x_{i+1})\right)$

* $L_{i,i+1}(x)\leq \log f(x)$ en el intervalo 
  $I_{i}=(x_{i},x_{i+1}]$

* $L_{i,i+1}(x)\geq \log f(x)$ fuera de $I_{i}$

En el intervalo $I_{i}$ se definen las envolventes de $\log f\left(
x\right)$:

* $\underline{\phi}_n(x)=L_{i,i+1}(x)$

* $\overline{\phi}_n(x)=\min \left\{L_{i-1,i}(x),L_{i+1,i+2}(x)\right\}$

Las envolventes de $f(x)$ en $I_{i}$ serán:

* $s_n(x)=\exp \left( \underline{\phi}_n(x)\right)$

* $G_n(x)=\exp \left( \overline{\phi}_n(x)\right)$


Tenemos entonces que:
$$s_n(x)\leq f(x) \leq G_n(x)=c\cdot g_n(x)$$
donde $g_n(x)$ es una mixtura discreta de distribuciones tipo exponencial truncadas (las tasas pueden ser negativas), que se puede simular fácilmente combinando el método de composición (Sección \@ref(composicion)) con el método de inversión.

::: {.conjecture #gilks name="Gilks 1992"}
<br>

1.  Inicializar $n$ y $s_n$.

2.  Generar $U \sim \mathcal{U}(0, 1)$ y
    $T\sim g_n$.

3.  Si $U\cdot G_n\left( T\right) \leq s_n\left( T\right)$
    devolver $X=T$,

    en caso contrario, 

    3.a Si $U\cdot G_n\left( T\right) \leq f\left( T\right)$
        devolver $X=T$.

    3.b Hacer $n=n+1$, añadir $T$ a $S_n$ 
        y actualizar $s_n$ y $G_n$.

4.  Volver al paso 2.

:::

Gilks y Wild (1992) propusieron una ligera modificación empleando tangentes para construir la cota superior, de esta forma se obtiene un método más eficiente pero requiere especificar la derivada de la densidad objetivo (ver Figura \@ref(fig:squeeze)).

La mayoría de las densidades de la familia exponencial de distribuciones son log-cóncavas.
Hörmann (1995) extendió esta aproximación al caso de densidades $T_{c}$-cóncavas:
$$T_{c}(x) = signo(c)x^{c} \\ T_{0}(x) = \log (x).$$
Aparte de la transformación logarítmica, la transformación $T_{-1/2}(x)=-1/\sqrt{x}$ es habitualmente la más empleada.

### Método del cociente de uniformes

Se puede ver como una modificación del método de aceptación-rechazo, de especial interés cuando el soporte no es acotado.

Si $(U,V)$ se distribuye uniformemente sobre:
$$C_{f^{\ast}} = \left\{ (u, v) \in \mathbb{R}^{2} : 
0<u\leq \sqrt{f^{\ast}(v/u)} \right\},$$
siendo $f^{\ast}$ una función no negativa integrable (cuasi-densidad), entonces $X=V/U$ tiene función de densidad proporcional a $f^{\ast}$ (Kinderman y Monahan, 1977). 
Además $C_{f^{\ast}}$ tiene área finita, por lo que pueden generarse fácilmente los valores $(U,V)$ con distribución $\mathcal{U}\left(C_{f^{\ast}}\right)$ a partir de componentes uniformes unidimensionales (aceptando los puntos dentro de $C_{f^{\ast}}$).

De modo análogo al método de aceptación-rechazo, hay modificaciones para acelerar los cálculos y automatizar el proceso, construyendo regiones mediante polígonos:
$$C_{i}\subset C_{f^{\ast}}\subset C_{s}.$$
También se puede extender al caso multivariante y considerar transformaciones adicionales. 
Ver por ejemplo el paquete [`rust`](https://paulnorthrop.github.io/rust/index.html).


::: {.example #cauchy-rou name="simulación de la distribución de Cauchy mediante cociente de uniformes"}
<br>

Si consideramos la distribución de Cauchy:
$$f(x) = \frac{1}{\pi (1 + x^2)} \text{, } x\in \mathbb{R},$$
eliminando la constante por comodidad $f(x) \propto 1/(1 + x^2)$, se tiene que:
$$\begin{aligned}
C_{f^{\ast}} & = \left\{ (u, v) \in \mathbb{R}^{2} : 0 <u \leq \frac{1}{\sqrt{1 + (v/u)^2}}  \right\} \\
& = \left\{ (u, v) \in \mathbb{R}^{2} : u > 0, u^2 \leq \frac{u^2}{u^2 + v^2}  \right\} \\
& = \left\{ (u, v) \in \mathbb{R}^{2} : u > 0, u^2 + v^2 \leq 1  \right\}, 
\end{aligned}$$
dando como resultando que $C_{f^{\ast}}$ es el semicírculo de radio uno, y podemos generar valores con distribución uniforme en esta región a partir de $\mathcal{U}\left([0,1]\times[-1,1] \right)$.


El correspondiente algoritmo está implementado en la función  [`rcauchy.rou()`](https://rubenfcasal.github.io/simres/reference/rcauchy.rou.html) del paquete [`simres`](https://rubenfcasal.github.io/simres) (fichero [*ar.R*](R/ar.R)):


```r
simres::rcauchy.rou
```

```
## function (n) 
## {
##     ngen <- n
##     u <- runif(n, 0, 1)
##     v <- runif(n, -1, 1)
##     x <- v/u
##     ind <- u^2 + v^2 > 1
##     while (le <- sum(ind)) {
##         ngen <- ngen + le
##         u <- runif(le, 0, 1)
##         v <- runif(le, -1, 1)
##         x[ind] <- v/u
##         ind[ind] <- u^2 + v^2 > 1
##     }
##     attr(x, "ngen") <- ngen
##     return(x)
## }
## <bytecode: 0x000000003f43bb88>
## <environment: namespace:simres>
```

```r
set.seed(1)
nsim <- 10^4
rx <- simres::rcauchy.rou(nsim)

hist(rx, breaks = "FD", freq = FALSE, main = "", xlim = c(-6, 6))
curve(dcauchy, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{04-Metodos_generales_continuas_files/figure-latex/rcauchy-rou-1} 

}

\caption{Distribución de los valores generados mediante el método de cociente de uniformes.}(\#fig:rcauchy-rou)
\end{figure}

```r
# Número generaciones
ngen <- attr(rx, "ngen")
{cat("Número de generaciones = ", ngen)
cat("\nNúmero medio de generaciones = ", ngen/nsim)
cat("\nProporción de rechazos = ", 1-nsim/ngen,"\n")}
```

```
## Número de generaciones =  12751
## Número medio de generaciones =  1.2751
## Proporción de rechazos =  0.2157478
```

:::


## Método de composición (o de simulación condicional) {#composicion}

En ocasiones la densidad de interés se puede expresar como una mixtura discreta de densidades:
$$f(x)=\sum_{j=1}^{k}p_{j}f_{j}(x)$$
con $\sum_{j=1}^{k}p_j=1$, $p_j\geq 0$ y $f_j$ densidades (sería también válido para funciones de distribución, incluyendo el caso discreto).

::: {.conjecture #mixtura-discreta name="simulación de una mixtura discreta"}
<br>

1. Generar $J$ con distribución $P\left( J=j \right) = p_j$.

2. Generar $X\sim f_J$.

:::

::: {.example #dexp-mix name="distribución doble exponencial"}
<br>

A partir de la densidad de la distribución doble exponencial:
$$f(x) =\frac{\lambda }{2}e^{-\lambda \left\vert x\right\vert } \text{, }\forall x\in \mathbb{R},$$
se deduce que:
$$f(x) = \frac{1}{2}f_{1}(x) + \frac{1}{2}f_{2}(x)$$
siendo:
$$f_{1}(x) = \left\{ 
\begin{array}{ll}
\lambda e^{-\lambda x} & \text{si } x\geq 0 \\ 
0 & \text{si } x<0
\end{array}
\right., \  
f_{2}(x) = \left\{ 
\begin{array}{ll}
\lambda e^{\lambda x} & \text{si } x<0 \\ 
0 & \text{si } x\geq 0
\end{array}
\ \right.$$

El algoritmo resultante sería el siguiente (empleando dos números pseudoaleatorios uniformes, el primero para seleccionar el índice y el segundo para generar un valor de la correspondiente componente mediante el método de inversión):

1. Generar $U,V \sim \mathcal{U}(0, 1)$.

2. Si $U<0.5$ devolver $X= -\ln( 1-V )/\lambda$.

3. En caso contrario devolver $X= \ln(V)/\lambda$.

Este método está implementado en la función [`rdexp()`](https://rubenfcasal.github.io/simres/reference/ddexp.html) del paquete [`simres`](https://rubenfcasal.github.io/simres) (fichero [*ddexp.R*](R/ddexp.R)).

:::


Observaciones:

* En ocasiones se hace un reciclado de los números aleatorios
  (solo se genera una uniforme, e.g. $V=2(U-0.5)$ si
  $U\in (0.5,1)$).

* En ciertas ocasiones por comodidad, para simular una muestra de
  tamaño $n$, se simulan muestras de tamaño $np_{i}$ con densidad
  $f_{i}$ y se combinan aleatoriamente.


Otro ejemplo de una mixtura discreta es el estimador tipo núcleo de la densidad (ver e.g. la ayuda de la función `density()` de R o la Sección \@ref(modunif-boot-suav)).
Simular a partir de una estimación de este tipo es lo que se conoce como *bootstrap suavizado*. 

En el caso de una mixtura continua tendríamos:
$$f(x)=\int g(x|y)h(y)dy$$

::: {.conjecture #mixtura-continua name="simulación de una mixtura continua"}
<br>

1. Generar $Y\sim h$.

2. Generar $X\sim g(\cdot |Y)$.

:::

Este algoritmo es muy empleado en Inferencia Bayesiana y en la simulación de algunas variables discretas (como la Binomial Negativa, denominada también distribución Gamma-Poisson, o la distribución Beta-Binomial; ver Sección \@ref(notables-disc)),
ya que el resultado sería válido cambiando las funciones de densidad $f$ y $g$ por funciones de masa de probabilidad.

## Métodos específicos para la generación de algunas distribuciones notables {#notables-cont}

En el pasado se ha realizado un esfuerzo considerable para desarrollar métodos eficientes para la simulación de las distribuciones de probabilidad más importantes.
Estos algoritmos se describen en la mayoría de los libros clásicos de simulación (e.g. Cao, 2002, Capítulo 5)^[Cuidado con la notación y la parametrización empleadas, puede variar entre referencias. Por ejemplo, en Cao (2002) la notación de la distribución Gamma es ligeramente distinta a la empleada en R y en el presente libro.], principalmente porque resultaba necesario implementar estos métodos durante el desarrollo de software estadístico.
Hoy en día estos algoritmos están disponibles en numerosas bibliotecas y no es necesario su implementación (por ejemplo, se puede recurrir a R o emplear su librería matemática disponible en <https://svn.r-project.org/R/trunk/src/nmath>).
Sin embargo, además de que muchos de ellos servirían como ilustración de la aplicación de los métodos generales expuestos en secciones anteriores, pueden servir como punto de partida para la generación de otras distribuciones.

Entre los distintos métodos disponibles para la generación de las distribuciones continuas más conocidas podríamos destacar:

* Método de Box-Müller para la generación de normales independientes (que se puede generalizar para otras distribuciones o incorporar dependencia).

* Algoritmos de Jöhnk (1963) y Cheng (1978) para la generación de la distribución beta (como ejemplo de la eficiencia de los métodos de aceptación-rechazo).

### Método de Box-Müller

Se basa en la siguiente propiedad. Dadas dos variables aleatorias independientes $E \sim \exp\left(  1\right)$ y
$U \sim \mathcal{U}( 0, 1 )$, las variables
$\sqrt{2E} \cos 2\pi U$ y $\sqrt{2E}\operatorname{sen} 2\pi U$ son
$\mathcal{N}( 0, 1 )$ independientes.

::: {.conjecture #box-muller name="de Box-Müller 1958"}
<br>

1. Generar $U,V\sim \mathcal{U}(0, 1)$.

2. Hacer $W_1=\sqrt{-2\ln U}$ y $W_2=2\pi V$.

3. Devolver $X_1=W_1\cos W_2$, $X_2=W_1\operatorname{sen}W_2$.

:::

Podemos hacer que la función `rnorm()` de R emplee este algoritmo estableciendo el parámetro `normal.kind` a `"Box-Muller"` en una llamada previa a `set.seed()` o `RNGkind()`.

Este método está relacionado con el denominado *método FFT* (transformada de Fourier; e.g. Davies y Harte, 1987) para la generación de una normal multidimensional con dependencia, que resulta ser equivalente al *Circular embedding* (Dietrich and Newsam, 1997).
La idea de estos métodos es que, considerando módulos exponenciales y fases uniformes generamos normales independientes, pero cambiando la varianza de los módulos ($W_1$) podemos inducir dependencia.
Adicionalmente, cambiando la distribución de las fases ($W_2$) se generan distribuciones distintas de la normal.

### Simulación de la distribución beta

Existen multitud de algoritmos para simular la distribución $\mathcal{Beta}(a, b)$. 
Probablemente, el más sencillo de todos es el que se obtiene a partir de la distribución gamma o de la chi-cuadrado, si se dispone de un algoritmo para generar estas distribuciones,  empleando la definición habitual de la distribución beta:

Si $Y \sim \mathcal{Gamma}(a, s)$ y $Z \sim  \mathcal{Gamma}(b, s)$ son independientes, entonces $$X=\frac{Y}{Y+Z} \sim \mathcal{Beta}(a, b).$$

Como la distribución resultante no depende de $s$  y $\chi^2_{n} \overset{d}{=} \mathcal{Gamma}\left(\tfrac{n}{2}, \tfrac{1}{2}\right)$, se podría considerar $Y \sim \chi^2_{2a}$ y $Z \sim \chi^2_{2b}$ independientes.

También se podrían emplear resultado conocidos relacionados con esta distribución, como por ejemplo que la distribución del estadístico de orden $k$ de una muestra de tamaño $n$ de una distribución uniforme tiene una distribución beta:
$$U_{(k)} \sim \mathcal{Beta}(k,n+1-k).$$
El resultado es el algoritmo de Fox (1963), que podría ser adecuado para simular esta distribución cuando $a, b \in \mathbb{N}$ y son valores pequeños.

::: {.conjecture #fox name="de Fox 1963"}
<br>
  
1. Generar $U_1, U_2, \ldots, U_{a+b-1} \sim \mathcal{U}(0, 1)$.

2. Ordenar: $U_{(1)}\leq U_{(2)}\leq\cdots\leq U_{(a+b-1)}$.

3. Devolver $X=U_{(a)}$.

:::

Es obvio que este algoritmo puede resultar muy lento si alguno de los dos parámetros es elevado (pues habrá que simular muchas uniformes para conseguir un valor simulado de la beta). 
Además, en función de cuál de los dos parámetros, $a$ ó $b$, sea mayor, resultará más eficiente, en el paso 2, comenzar a ordenar por el mayor, luego el segundo mayor, y así sucesivamente, o hacerlo empezando por el menor. 
En cualquier caso, es obvio que no es necesario ordenar todos los valores $U_{i}$ generados, sino tan sólo encontrar el que ocupa el lugar $a$-ésimo.

Un método válido aunque $a$ ó $b$ no sean enteros es el dado por el
algoritmo de Jöhnk (1964).

::: {.conjecture #johnk name="de Jöhnk 1964"}
<br>
  
1. Generar $U_1, U_2\sim \mathcal{U}(0, 1)$.
    
2. Hacer $V = U_1^{\frac1a}$, $W = U_2^{\frac1b}$ y $S = V+W$.

3. Si $S \leq 1$ devolver $X = \frac VS$,

    en caso contrario volver al paso 1.

:::

El método resulta extremadamente ineficiente para $a$ ó $b$ mayores que 1. 
Esto es debido a que la condición $S\leq1$ del paso 3 puede tardar muchísimo en verificarse. 
Por este motivo, el algoritmo de Jöhnk sólo es recomendable para $a<1$ y $b<1$. 
Como remedio a esto puede usarse el algoritmo de Cheng (1978) que es algo más complicado de implementar^[R implementa este algoritmo en el fichero fuente [rbeta.c](https://svn.r-project.org/R/trunk/src/nmath/rbeta.c).] pero mucho más eficiente.

::: {.conjecture #cheng name="de Cheng 1978"}
<br>
  
Inicialización:

1. Hacer $\alpha = a + b$.

2. Si $\min(a,b) \leq1$ entonces hacer $\beta=\frac1{\min(  a,b)}$, en otro caso hacer $\beta=\sqrt{\frac{\alpha-2}{2pq-\alpha}}$.

3. Hacer $\gamma=a+\frac1\beta$.

Simulación:

1. Generar $U_1, U_2\sim \mathcal{U}(0, 1)$.

2. Hacer $V=\beta\cdot\ln\left(  \frac{U_1}{1-U_1}\right)$ y $W=a\cdot e^{V}$.

3. Si $\alpha\cdot\ln\left( \frac\alpha{b+W}\right) +\gamma V-\ln4 \ge \ln\left(  U_1^{2}U_2\right)$ devolver $X=\frac W{b+W}$,

    en caso contrario volver al paso 1.

:::

