# Introducción a la simulación {#cap1}




Podríamos definir la *Simulación* como una técnica que consiste en realizar experimentos de muestreo sobre el modelo de un sistema, con el objetivo de recopilar información bajo determinadas condiciones. 

## Conceptos básicos

La experimentación directa sobre la realidad puede tener muchos inconvenientes:

-   Coste elevado.

    -   En ocasiones las pruebas son destructivas.
    
    -   Lentitud.

-   Puede no ser ética.

    -   Experimentación sobre seres humanos.

-   Puede resultar imposible.

    -   Acontecimientos futuros, alternativas en el pasado, ...

-   ...

Además la realidad puede ser demasiado compleja como para ser estudiada directamente y resultar preferible trabajar con un modelo del sistema real. 
Un modelo no es más que un conjunto de variables junto con ecuaciones matemáticas que las relacionan y restricciones sobre dichas variables. 
Habría dos tipos de modelos:

-   Modelos deterministas.

-   Modelos estocásticos (con componente aleatoria): tienen en cuenta la incertidumbre debida a que no se disponga de toda la información sobre las variables que influyen en el fenómeno en estudio (puede ser simplemente que haya errores de medida).

> "Nothing in Nature is random... a thing appears random only through the incompleteness of our knowledge."
>
> --- Spinoza, Baruch (Ethics, 1677)

La modelización es una etapa presente en la mayor parte de los trabajos de investigación (especialmente en las ciencias experimentales).
El modelo debería considerar las variables más relevantes para explicar el fenómeno en estudio y las principales relaciones entre ellas.
La inferencia estadística proporciona herramientas para estimar los parámetros y contrastar la validez de un modelo estocástico a partir de los datos observados.

La idea es emplear el modelo (suponiendo que es válido) para resolver el problema de interés. 
Si se puede obtener la solución de forma análitica, esta suele ser exacta (aunque en ocasiones solo se dispone de soluciones aproximadas, basadas en resultados asintóticos, o que dependen de supociones que pueden ser cuestionables) y a menudo la resolución también es rápida.
Cuando la solución no se puede obtener de modo analítico (o si la aproximación disponible no es adecuada) se puede recurrir a la simulación.

Nos centraremos en el caso de la *Simulación Estocástica*: las conclusiones se obtienen generando repetidamente simulaciones del modelo aleatorio.

<!-- 
Ejemplo: caballero de Meré 
Experimentación directa sobre la realidad
Modelo de probabilidad
Aproximación por simulación
-->

### Ventajas e inconvenientes de la simulación 

Ventajas (Shannon, 1975):

-   Cuando la resolución analítica no puede llevarse a cabo.

-   Cuando existen medios de resolver analíticamente el problema
    pero dicha resolución es complicada y costosa 
    (o solo proporciona una solución aproximada).

-   Si se desea experimentar antes de que exista el sistema 
    (pruebas para la construcción de un sistema).

-   Cuando es imposible experimentar sobre el sistema real 
    por ser dicha experimentación destructiva.

-   En ocasiones en las que la experimentación sobre el sistema es
    posible pero no ética.

-   En sistemas que evolucionan muy lentamente en el tiempo.


El principal incoveniente puede ser el tiempo de computación necesario, aunque gracias a la gran potencia de cálculo de los computadores actuales, se puede obtener rápidamente una solución aproximada en la mayor parte de los problemas susceptibles de ser modelizados.
Además siempre están presentes los posibles problemas debidos a emplear un modelo:

-   La construcción de un buen modelo puede ser una tarea muy costosa 
    (compleja, laboriosa y requerir mucho tiempo; 
    e.g. modelos climáticos).

-   Frecuentemente el modelo omite variables o relaciones importantes entre ellas
    (los resultados pueden no ser válidos para el sistema real).

-   Resulta difícil conocer la precisión del modelo formulado.

Otro problema de la simulación es que se obtienen resultados para unos valores concretos de los parámetros del modelo, por lo que en principio
resultaría complicado extrapolar las conclusiones a otras situaciones.


## Generación de números (pseudo)aleatorios

Aplicaciones:

-   Estadística:

    -   Muestreo, remuestreo, ...
    
    -   Aproximación de distribuciones (de estadísticos, estimadores, ...)
    
    -   Realización de contrastes, intervalos de confianza, ...
    
    -   Comparación de estimadores, contrastes, ...
    
    -   Validación teoría (distribución asintótica,...)
    
    -   Inferencia Bayesiana

-   Optimización: Algoritmos genéticos, ...

-   Análisis numérico: Aproximación de integrales, resolución de ecuaciones, ...

-   Computación: Diseño, verificación y validación de algoritmos,...

-   Criptografía: Protocolos de comunicación segura, ...

-   Física: Simulación de fenómenos naturales, ...

-   ...


En los capítulos \@ref(cap8) y \@ref(cap9) nos centraremos en algunas de las aplicaciones de utilidad en Estadística.


## Números aleatorios puros

Una sucesión de números aleatorios puros (*true random*), se caracteriza por que no existe ninguna regla o plan que nos permita conocer sus valores.

Normalmente son obtenidos por procesos físicos
(loterías, ruletas, ruidos,...) y se almacenaban en *tablas de dígitos aleatorios*. 
Por ejemplo, en 1955 la Corporación RAND publicó el libro [A Million Random Digits with 100,000 Normal Deviates](https://www.rand.org/pubs/monograph_reports/MR1418.html) que contenía números aleatorios generados mediante una ruleta electrónica conectada a una computadora  (ver Figura \@ref(fig:randbook)).


\begin{figure}[!htb]

{\centering \includegraphics[width=0.3\linewidth]{images/RAND} 

}

\caption{Líneas 10580-10594, columnas 21-40, del libro *A Million Random Digits with 100,000 Normal Deviates*.}(\#fig:randbook)
\end{figure}

El procediento para generar de forma manual números aleatorios 
en un rango de 1 a *m* era el siguiente:

-   Se selecciona al azar un pto de inicio en la tabla 
    y una dirección.

-   Se agrupan los dígitos de forma que “cubran” el valor de *m*.

-   Se seleccionan los valores menores o iguales que *m* 
    (se descarta el resto).

Hoy en día están disponibles generadores de números aleatorios “online”:

-   [http://www.random.org/integers](http://www.random.org/integers)
    (ver paquete `random` en R).

-   [http://www.fourmilab.ch/hotbits](http://www.fourmilab.ch/hotbits)

Aunque para un uso profesional puede ser recomendable emplear generadores implementados mediante hardware:

-   [http://software.intel.com](http://software.intel.com/en-us/articles/intel-digital-random-number-generator-drng-software-implementation-guide).

-   [http://spectrum.ieee.org](http://spectrum.ieee.org/semiconductors/devices/a-memristor-true-randomnumber-generator)


Sus principales aplicaciones hoy en día son en criptografía (impredictibilidad).


### Inconvenientes:

Para su aplicación en el campo de la Estadística:

-   Es necesario conocer su distribución.

-   Los datos generados deberían ser i.i.d.

(sería también muy recomendable para otros casos) y en general:

-   Reproductivilidad.

-   Pueden requerir del almacenamiento en tablas.

Además siempre está la posible presencia de sesgos, principalmente debidos a fallos del sistema o interferencias. 
Por ejemplo en el caso de la máquina RAND, fallos mecánicos en el sistema de grabación de los datos causaron problemas de aleatoriedad (Hacking, 1965, p. 129).


### Alternativas:

A partir de la década de 1960, al disponer de computadoras de mayor velocidad, empezó a resultar más eficiente generar valores mediante software en lugar de leerlos de las tablas. 
Se distingue entre dos tipos de secuencias:

-   *números pseudo-aleatorios*: simulan realizaciones de una variable aleatoria (uniforme).

-   números cuasi-aleatorios: secuencias determinísticas con una distribución más uniforme en el rango considerado (se podría pensar que son una única generación de una variable aleatoria).

Algunos problemas, como la integración numérica (en el Capítulo \@ref(cap9) se tratarán métodos de integración Monte Carlo), no dependen realmente de la aleatoriedad de la secuencia. Para evitar generaciones poco probables, se puede recurrir a secuencias cuasi-aleatorias, también denominadas *sucesiones de baja discrepancia* (hablaríamos entonces de métodos cuasi-Monte Carlo). La idea sería que la proporción de valores en una región cualquiera sea siempre aproximadamente proporcional a la medida de la región (como sucedería en media con la distribución uniforme, aunque no necesariamente para una realización concreta).

Por ejemplo, el paquete [`randtoolbox`](https://CRAN.R-project.org/package=randtoolbox) implementa métodos para la generación de secuencias cuasi-aleatorias (ver Figura \@ref(fig:randtoolbox)).


```r
library(randtoolbox)
n <- 2000
par.old <- par( mfrow=c(1,3))
plot(halton(n, dim = 2), xlab = 'x1', ylab = 'x2')
plot(sobol(n, dim = 2), xlab = 'x1', ylab = 'x2')
plot(torus(n, dim = 2), xlab = 'x1', ylab = 'x2')
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{01-Introduccion_files/figure-latex/randtoolbox-1} 

}

\caption{Secuencias cuasi-aleatorias bidimensionales obtenidas con los métodos de Halton (izquierda), Sobol (centro) y Torus (derecha).}(\#fig:randtoolbox)
\end{figure}

```r
par(par.old)
```

En este libro sólo consideraremos los números pseudoaleatorios y por comodidad se eliminará el prefijo "pseudo" en algunos casos.


Números pseudoaleatorios
------------------------

La mayoría de los métodos de simulación se basan en la posibilidad de generar números pseudoaleatorios que imiten las propiedades de generaciones independientes de una distribución $\mathcal{U}(0,1)$. 

El procedimiento habitual para obtiener estas secuencias es emplear un algoritmo recursivo denominado *generador*:

$$x_{i} = f\left( x_{i-1}, x_{i-2}, \cdots, x_{i-k}\right)$$

donde:

-   $k$ es el orden del generador.

-   $\left(  x_{0},x_{1},\cdots,x_{k-1}\right)$ es la *semilla*
  (estado inicial).

El *periodo* o *longitud del ciclo* es la longitud de la secuencia antes
de que vuelva a repetirse. Lo denotaremos por $p$.


Los números de la sucesión serán predecibles, conociendo el
algoritmo y la semilla.

-   Sin embargo, si no se conociesen, *no se debería poder distinguir* una serie de números pseudoaleatorios *de una sucesión de números verdaderamente aleatoria* (utilizando recursos computacionales razonables).

-   En caso contrario esta predecibilidad puede dar lugar a serios
problemas (e.g. [http://eprint.iacr.org/2007/419](http://eprint.iacr.org/2007/419)).

Como regla general, por lo menos mientras se está desarrollando un
programa, interesa *fijar la semilla de aleatorización*.

-   Permite la reproductibilidad de los resultados.

-   Facilita la depuración del código.

Todo generador de números pseudoaleatorios mínimamente aceptable
debe comportarse como si se tratase de una muestra genuina de datos
independientes de una $\mathcal{U}(0,1)$.
Otras propiedades de interés serían:

-   Reproducibilidad a partir de la semilla.

-   Periodo suficientemente largo.

-   Eficiencia (rapidez y requerimientos de memoria).

-   Portabilidad.

-   Generación de sub-secuencias (computación en paralelo).

-   Parsimonia, ...


> "... random numbers should not be generated with a method chosen at random." 
>
>
> --- Knuth, D.E. (TAOCP, 2002)


Gran cantidad de algoritmos:

-   Cuadrados medios, Lehmer,...

-   Congruenciales

-   Registros desfasados

-   Combinaciones

-   ...

Código fuente disponible en múltiples librerias:

-   GNU Scientific Library (GSL):
    [http://www.gnu.org/software/gsl/manual](http://www.gnu.org/software/gsl/manual/html\_node/Random-Number-Generation.html)

-   StatLib: [http://lib.stat.cmu.edu](http://lib.stat.cmu.edu)

-   Numerical recipes: [http://www.nrbook.com/nr3](http://www.nrbook.com/nr3)

-   [http://random.mat.sbg.ac.at/software](http://random.mat.sbg.ac.at/software)

-   KISS (Keep It Simple Stupid / Small and Simple):
    [http://www.fortran.com/kiss.f90](http://www.fortran.com/kiss.f90)

-   UNU.RAN (paquete `Runuran`):
    [http://statmath.wu.ac.at/unuran](http://statmath.wu.ac.at/unuran)

-   ...
  
Nos centraremos en los generadores congruenciales, descritos en la Sección \@ref(gen-cong).
