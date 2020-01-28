Introducción a la simulación
============================



## Introducción

  Problemas de la experimentación directa sobre la realidad:

  -   Coste elevado.

-   En ocasiones las pruebas son destructivas.

-   Lentitud.

-   Puede no ser ética.

-   experimentación sobre seres humanos.

-   Puede resultar imposible.

-   Acontecimientos futuros, ...

-   ...

Puede ser preferible trabajar con un modelo del sistema real. La realidad puede ser muy compleja por lo que es habitual emplear un
modelo para tratar de explicarla:

-   Modelos deterministas.

-   Modelos estocásticos (con componente aleatoria).

    -   Cuando no se dispone de la suficiente información sobre las
        variables que influyen en el fenómeno en estudio.
    
    -   Se tienen en cuenta esta incertidumbre.

La inferencia estadística proporciona herramientas para
estimar los parámetros y contrastar la validez de un modelo estocástico a
partir de los datos observados.

La idea es emplear el modelo para resolver el problema de interés. 
Cuando la solución no se puede obtener de modo
analítico (teórico) se puede recurrir a la simulación.

**Simulación**: realizar experimentos con un modelo con el objetivo
de recopilar información bajo determinadas condiciones. Nos centraremos en el caso de la **simulación estocástica**: 

-   Las conclusiones se obtienen habitualmente generando
    repetidamente simulaciones del modelo aleatorio.
    


### Ventajas de la simulación (Shannon, 1975):

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


### Inconvenientes de la simulación:

-   La construcción de un buen modelo puede ser una tarea muy costosa 
    (compleja, laboriosa y requerir mucho tiempo; e.g. modelos climáticos).

-   Frecuentemente el modelo omite variables o relaciones importantes entre ellas
    (los resultados pueden no ser válidos para el sistema real).

-   Resulta difícil conocer la precisión de la simulación
    (especialmente en lo relativo a la precisión del modelo formulado).

-   Problemas de extrapolación de las conclusiones.

-   Tiempo de computación.


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

-   Computación: Diseño, verificación y validación de algoritmos,...

-   Criptografía: Protocolos de comunicación segura, ...

-   Física: Simulación de fenómenos naturales, ...

-   Análisis numérico: Evaluación de expresiones, ...

-   ...


En el Capítulo XX nos centraremos en algunas de las aplicaciones de utilidad en Estadística.


## Números aleatorios puros

Una sucesión de números aleatorios puros (true random), se
caracteriza por que no existe ninguna regla o plan que nos permita
conocer sus valores.

Se almacena(ban) en *tablas de dígitos aleatorios* (*true random*)
y normalmente son obtenidos por procesos físicos
(loterías, ruletas, ruidos,...)


\begin{center}\includegraphics[width=0.7\linewidth]{images/RAND} \end{center}

Se emplean para seleccionar números aleatorios 
en un rango de 1 a *m*:

-   Se selecciona al azar un pto de inicio en la tabla 
    y una dirección.

-   Se agrupan los dígitos de forma que “cubran” el valor de *m*.

-   Se seleccionan los valores menores o iguales que *m* 
    (se descarta el resto).

Algunos enlaces:

-   [A Million Random Digits with 100,000 Normal Deviates](https://en.wikipedia.org/wiki/A_Million_Random_Digits_with_100,000_Normal_Deviates) .
    RAND Corporation. 1955.

-   Generadores de números aleatorios “online”:

    -   [http://www.random.org/integers](http://www.random.org/integers)
        (ver paquete **random** en R).
    
    -   [http://www.fourmilab.ch/hotbits](http://www.fourmilab.ch/hotbits)

-   Generadores mediante hardware:

    -   [http://software.intel.com](http://software.intel.com/en-us/articles/intel-digital-random-number-generator-drng-software-implementation-guide).
    
    -   [http://spectrum.ieee.org](http://spectrum.ieee.org/semiconductors/devices/a-memristor-true-randomnumber-generator)


  Sus principales aplicaciones hoy en día son en criptografía (impredictibilidad).

### Inconvenientes:

  -   Es necesario/recomendable conocer su distribución.

-   Los datos generados deberían ser i.i.d.

-   Reproductivilidad.

-   Pueden requerir del almacenamiento en tablas.

### Alternativas:

  -   **números pseudo-aleatorios**: simulan realizaciones de una
variable aleatoria (uniforme).

-   números cuasi-aleatorios: secuencias determinísticas con una
distribución más uniforme en el rango considerado (se podría
                                                   pensar que son una única generación de una variable aleatoria).


Números pseudoaleatorios
------------------------

### Generación de números pseudoaleatorios mediante software

La mayoría de los métodos de simulación se basan en la posibilidad
de generar números pseudoaleatorios con distribución
$\mathcal{U}(0,1)$. 
Se obtienen mediante un algoritmo recursivo denominado
**generador**:

  $$x_{i}=f\left(  x_{i-1},x_{i-2},\cdots,x_{i-k}\right)$$

  -   $k$ orden del generador.

-   $\left(  x_{0},x_{1},\cdots,x_{k-1}\right)$ **semilla**
  (estado inicial).

El **periodo** o *longitud del ciclo* es la longitud de la secuencia antes
de que vuelva a repetirse. Lo denotaremos por $p$.


Los números de la sucesión serán predecibles, conociendo el
algoritmo y la semilla.

-   Sin embargo, si no se conociesen, **no se debería poder
distinguir** una serie de números pseudoaleatorios **de una
sucesión de números verdaderamente aleatoria** (utilizando
                                                recursos computacionales razonables).

-   En caso contrario esta predecibilidad puede dar lugar a serios
problemas (e.g. [http://eprint.iacr.org/2007/419](http://eprint.iacr.org/2007/419)).

Como regla general, por lo menos mientras se está desarrollando un
programa, interesa **fijar la semilla de aleatorización**.

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
  
Nos centraremos en los generadores congruenciales, descritos en la [Sección 3.1][Generadores congruenciales].
