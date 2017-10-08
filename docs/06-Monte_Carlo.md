Integración y Optimización Montecarlo
=====================================




***Pendiente: ejemplos y ejercicios***

Integración Monte Carlo (clásica)
-----------------------

La integración Monte Carlo se emplea principalmente para aproximar
integrales multidimensionales de la
forma:
$$I = \int \cdots \int h\left( x_{1},\ldots ,x_{d}\right) 
dx_{1}\cdots dx_{d}$$ 
donde puede presentar ventajas respecto a los métodos tradicionales de
integración numérica.

Supongamos que nos interesa:
$$I = \int_{0}^{1}h\left( x\right) dx$$
Si $x_{1},x_{2},\ldots ,x_{n}$ *i.i.d.* $\mathcal{U}\left( 0,1\right)$
entonces:
$$I = E\left( h\left( \mathcal{U}\left( 0,1\right) \right) \right)
\approx \frac{1}{n}\sum\limits_{i=1}^{n}h\left( x_{i}\right)$$

Si el intervalo de integración es
genérico:
$$I = \int_{a}^{b}h\left( x\right) dx = 
(b-a)\int_{a}^{b}h\left( x\right) \frac{1}{(b-a)}dx.$$
Si $x_{1},x_{2},\ldots ,x_{n}$ *i.i.d.*
$\mathcal{U}\left( a, b\right)$:
$$I\approx \frac{1}{n}\sum\limits_{i=1}^{n}h\left( x_{i}\right) (b-a)$$

De forma más general, si nos interesa
aproximar:
$$\theta = E\left( h\left( X\right) \right) = \int h\left( x\right) f(x)dx$$siendo
$X\sim f$, entonces, si $x_{1},x_{2},\ldots ,x_{n}$ *i.i.d.*
$X$:
$$\theta \approx \frac{1}{n}\sum\limits_{i=1}^{n}h\left( x_{i}\right)$$

Muestreo por importancia
------------------------

Para aproximar la integral:
$$\theta = E\left( h\left( X\right) \right) = \int h\left( x\right) f(x)dx,$$
puede ser preferible generar observaciones de una densidad $g$ 
que tenga una forma similar al producto $hf$.

Si $Y\sim g$:
$$\theta  = \int h\left( x\right) f(x)dx 
 = \int \frac{h\left( x\right) f(x)}{g(x)}g(x)dx
 = E\left( q\left( Y\right) \right).$$
siendo
$q\left( x\right)  = \frac{h\left( x\right) f(x)}{g(x)}$.



Si $y_{1},y_{2},\ldots ,y_{n}$ *i.i.d.* $Y\sim g$:
$$\theta \approx \frac{1}{n}\sum\limits_{i=1}^{n}q\left( y_{i}\right) 
= \frac{1}{n}\sum\limits_{i=1}^{n}w(y_{i})h\left( y_{i}\right)  
= \hat{\theta}_{g}$$
con $w\left( x\right)  = \frac{f(x)}{g(x)}$.

En este caso
$Var(\hat{\theta}_{g}) = Var\left( q\left( Y\right) \right) /n$, 
pudiendo reducirse significativamente respecto al método
clásico si:
$$g(x)\underset{aprox.}{\propto }h\left( x\right) f(x).$$



Para aplicar el TCL, la varianza del estimador $\hat{\theta}_{g}$ 
es finita si:
$$E\left( q^{2}\left( Y\right) \right)  
= \int \frac{h^{2}\left( x\right)^{2}(x)}{g(x)}dx 
= \int h^{2}\left( x\right) \frac{f^{2}(x)}{g(x)}dx 
< \infty.$$

La idea básica es que si la densidad $g$ tiene colas más pesadas que
la densidad $f$ con mayor facilidad puede dar lugar a 
"simulaciones" con varianza finita 
(podría emplearse en casos en los que no existe
$E \left( h^{2} \left( X \right) \right)$;
ver Sección XX en el Capítulo XX).

La distribución de los pesos $w(y_{i})$ debería ser homogénea para
evitar datos influyentes (inestabilidad).

### Remuestreo (del muestreo) por importancia

Cuando $f$ y/o $g$ son cuasi-densidades, para evitar calcular ctes
normalizadoras, se emplea como
aproximación:
$$\theta \approx \frac{\sum\limits_{i=1}^{n}w(y_{i})h\left( y_{i}\right) }{ \sum\limits_{i=1}^{n}w(y_{i})}.$$

Adicionalmente, puede verse que con un muestreo de 
$\left\{y_{1},y_{2},\ldots ,y_{n}\right\}$ 
ponderado por $w(y_{i})$ 
(prob. $=w(y_{i})\left/ \sum\nolimits_{i=1}^{n}w(y_{i}) \right.$ ) 
se obtiene una simulación aproximada de $f$ 
(*Sample importance resampling*, Rubin, 1987).



Optimización Monte Carlo
------------------------

Supongamos que estamos interesados en la minimización de una
función:
$$\underset{\mathbf{x}\in D}{\min }f(\mathbf{x}).$$

Hay una gran cantidad de algoritmos numéricos para resolver
problemas de optimización no lineal multidimensional, por ejemplo
los basados en el método de Newton-Raphson.

La idea original consiste en buscar los ceros de su primera derivada
(o del gradiente) empleando una aproximación
iterativa:
$$\mathbf{x}_{i+1} = \mathbf{x}_{i}-[Hf(\mathbf{x}_{i})]^{-1}\nabla f(\mathbf{x} = 
_{i}),$$donde $Hf(\mathbf{x}_{i})$ es el hessiano de la función
(matriz de segundas derivadas) y $\nabla f(\mathbf{x}_{i})$ el
gradiente (vector de primeras derivadas).
Estos métodos normalmente funcionan muy bien cuando la función
objetivo no tiene mínimos locales (ideal $f$ cuadrática).
Los resultados obtenidos pueden ser muy malos en caso contrario
(especialmente en el caso multidimensional) y dependen en gran
medida del punto inicial.
Un ejemplo donde es habitual que aparezcan este tipo de problemas es
en la estimación por máxima verosimilitud (la función objetivo puede
ser multimodal).

Una alternativa sería tratar de generar valores aleatoriamente de
forma que las regiones donde la función objetivo es menor tuviesen
mayor probabilidad y menor probabilidad las regiones donde la
función objetivo es mayor.
Por ejemplo, se podría pensar en generar valores de acuerdo a una
densidad (tranformación
Boltzman-Gibbs):
$$g(x)\propto \exp \left( -f(x)/T\right) ,$$donde
$T>0$ es un parámetro (denominado temperatura) seleccionado de forma
que se garantize la integrabilidad.

Entre los métodos de optimización Monte Carlo podríamos destacar:

-   Métodos con gradiente aleatorio.

-   Temple simulado.

-   Algoritmos genéticos.

-   Montecarlo EM.

-   ...

Temple simulado
---------------

Método inspirado en el templado de un metal (se calienta el metal a
alta temperatura y se va enfriando lentamente).
En cada paso se reemplaza la aproximación actual por un valor
aleatorio “cercano”, elegido con una probabilidad que depende de la
mejora en la función objetivo y de un parámetro $T$
(denominado temperatura) que disminuye gradualmente durante
el proceso.

-   Cuando la temperatura es grande los cambios son bastante
    probables en cualquier dirección.

-   Al ir disminuyendo la temperatura los cambios tienden a ser
    siempre “cuesta abajo”.

Al tener una probabilidad no nula de aceptar una modificación
“cuesta arriba” se trata de evitar quedar atrapado en un
óptimo local.


Algoritmo:


```r
# Algoritmo
temp <- TEMP.INIT
place <- INIT.PLACEMENT()
cost.place <- COST(place)
while(temp < TEMP.FINAL) {
  while(LOOP.CRITERION()) {
    place.new <- PERTURB(place, temp)
    cost.new <- COST(place.new)
    cost.inc <- cost.new - cost.place
    temp <- SCHEDULE(temp)
    if ((cost.inc < 0) || (runif(1) > exp(-(cost.inc/temp)))) break
  }
  place <- place.new
  cost.place <- cost.new 
  # temp <- SCHEDULE(temp)
}
COST <- function(place, ...) {...}
SCHEDULE <- function(temp, ...) {...}
INIT.PLACEMENT <- function(...) {...}
LOOP.CRITERION <- function(...) {...}
```

Adaptado de [Premchand Akella (ppt)](http://www.ecs.umass.edu/ece/labs/vlsicad/ece665/slides/SimulatedAnnealing.ppt).

Este algoritmo se puede ver como una adaptación del método de
Metropolis-Hastings que veremos más adelante (Sección XX).


Algoritmos genéticos
--------------------

Los algoritmos genéticos tratan de encontrar la mejor solución
(entre un conjunto de soluciones posibles) imitando los procesos de
evolución biológica:

-   **Población**: formada por $n$ individuos $\mathbf{x}_{i}$
    codificados en **cromosomas**.

-   $f(\mathbf{x}_{i})$ ajuste/capacidad/**adaptación** del
    individuo $\mathbf{x}_{i}$.

-   **Selección**: los individuos mejor adaptados tienen mayor
    probabilidad de ser **padres**.

-   **Cruzamiento**: los cromosomas de dos padres se combinan para
    generar hijos.

-   **Mutación**: modificación al azar del cromosoma de los
    hijos (variabilidad).

-   **Elitismo**: el mejor individuo pasa a la siguiente generación.

Los paquetes de R `DEOptim` y `gafit` implementan algunos de estos
tipos de algoritmos.


