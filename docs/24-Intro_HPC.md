# Introducción al procesamiento en paralelo en R {#intro-hpc}

<!-- bookdown::preview_chapter("25-Intro_HPC.Rmd") -->

<!-- Apéndice \@ref(intro-hpc) -->




En este apéndice se pretenden mostrar las principales herramientas para el procesamiento en paralelo disponibles en R y dar una idea de su funcionamiento. 
Para más detalles se recomienda ver [CRAN Task View: High-Performance and Parallel Computing with R](https://cran.r-project.org/view=HighPerformanceComputing)
(además de HPC, High Performance Computing, también incluye herramientas para computación distribuida)^[También puede ser de interés la presentación [R y HPC (uso de R en el CESGA)](https://www.r-users.gal/sites/default/files/10_aurelio_rodriguez.pdf) de Aurelio Rodríguez en las [VI Xornada de Usuarios de R en Galicia](https://www.r-users.gal/pagina/programa-2018)].

## Introducción

Emplearemos la siguiente terminología:

-   **Núcleo**: término empleado para referirse a un procesador lógico de un equipo (un equipo puede tener un único procesador con múltiples núcleos que pueden realizar operaciones en paralelo). 
    También podría referirse aquí a un equipo (nodo) de una red (clúster de equipos; en la práctica cada uno puede tener múltiples núcleos). *Un núcleo puede ejecutar procesos en serie*.
    
-   **Clúster**: colección de núcleos en un equipo o red de equipos. 
    *Un clúster puede ejecutar varios procesos en paralelo*.


Por defecto la versión oficial de R emplea un único núcleo, aunque se puede compilar de forma que realice cálculos en paralelo (e.g. librería LAPACK).
También están disponibles otras versiones de R que ya implementan por defecto procesamiento en paralelo ([multithread](https://mran.revolutionanalytics.com/documents/rro/multithread)): 
-   [Microsoft R Open](https://mran.revolutionanalytics.com): versión de R con rendimiento mejorado.


Métodos simples de paralelización^[Realmente las herramientas estándar son *OpenMP* para el procesamiento en paralelo con memoria compartida en un único equipo y *MPI* para la computación distribuida en múltiples nodos.]:

-   *Forking*: Copia el proceso de R a un nuevo núcleo (se comparte el entorno de trabajo). 
    Es el más simple y eficiente pero **no está disponible en Windows**.

-   *Socket*: Lanza una nueva versión de R en cada núcleo, como si se tratase de un cluster de equipos comunicados a través de red (hay que crear un entorno de trabajo en cada núcleo).
    Disponible en todos los sistemas operativos.


## Paquetes en R

Hay varios paquetes que se pueden usar para el procesamiento paralelo en R,
entre ellos podríamos destacar:

-   `parallel`: forma parte de la instalación base de R y fusiona los paquetes `multicore` (forking) y `snow` (sockets; Simple Network of Workstations). 
    Además incluye herramientas para la generación de números aleatorios en paralelo (cada proceso empleara una secuencia y los resultados serán reproducibles). 

    Incluye versiones "paralelizadas" de la familia `*apply()`: 
    `mclapply()` (forking), `parLapply()` (socket), ...
    
-   `foreach`: permite realizar iteraciones y admite paralelización con el operador `%dopar%`, aunque requiere paquetes adicionales como `doSNOW` o `doParallel` (recomendado).

-   `rslurm`: permite la ejecución distribuida en clústeres Linux que implementen SLURM (Simple Linux Utility for Resource Management), un gestor de recursos de código abierto muy empleado.


## Ejemplos

Si se emplea el paquete `parallel` en sistemas tipo Unix (Linux, Mac OS X, ...), se podría evaluar en paralelo una función llamando directamente a `mclapply()`. 
Por defecto empleará todos los núcleos disponibles, pero se puede especificar un número menor mediante el argumento `mc.cores`. 


```r
library(parallel)
ncores <- detectCores()
ncores
```

```
## [1] 20
```

```r
func <- function(k) {
  i_boot <- sample(nrow(iris), replace = TRUE)
  lm(Petal.Width ~ Petal.Length, data = iris[i_boot, ])$coefficients
}

RNGkind("L'Ecuyer-CMRG") # Establecemos Pierre L'Ecuyer's RngStreams...
set.seed(1)

system.time(res.boot <- mclapply(1:100, func)) # En Windows llama a lapply() (mc.cores = 1)
```

```
##    user  system elapsed 
##    0.03    0.01    0.05
```

```r
# res.boot <- mclapply(1:100, func, mc.cores = ncores - 1) # En Windows genera un error si mc.cores > 1
```

En Windows habría que crear previamente un cluster, llamar a una de las funciones `par*apply()` y finalizar el cluster:


```r
cl <- makeCluster(ncores - 1, type = "PSOCK")
clusterSetRNGStream(cl, 1) # Establecemos Pierre L'Ecuyer's RngStreams con semilla 1...

system.time(res.boot <- parSapply(cl, 1:100, func))
```

```
##    user  system elapsed 
##    0.00    0.00    0.01
```

```r
# stopCluster(cl)

str(res.boot)
```

```
##  num [1:2, 1:100] -0.415 0.429 -0.363 0.42 -0.342 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:2] "(Intercept)" "Petal.Length"
##   ..$ : NULL
```

Esto también se puede realizar en Linux (`type = "FORK"`), aunque podríamos estar trabajando ya en un cluster de equipos...

También podríamos emplear balance de carga si el tiempo de computación es variable (e.g. `parLapplyLB()` o `clusterApplyLB()`) pero no sería recomendable si se emplean números pseudo-aleatorios (los resultados no serían reproducibles).

Además, empleando las herramientas del paquete `snow` se puede representar el uso del cluster (*experimental* en Windows):

```r
# library(snow)
ctime <- snow::snow.time(snow::parSapply(cl, 1:100, func))
ctime
plot(ctime)
```

Hay que tener en cuenta la sobrecarga adicional debida a la comunicación entre nodos al paralelizar (especialmente con el enfoque de socket).


### Procesamiento en paralelo con la función `boot()`

La función `boot::boot()` incluye parámetros para el procesamiento en paralelo: 
`parallel = c("no", "multicore", "snow")`, `ncpus`, `cl`.
Si `parallel = "snow"` se crea un clúster en la máquina local durante la ejecución, salvo que se establezca con el parámetro `cl`. 

Veamos un ejemplo empleando una muestra simulada:

```r
n <- 100
rate <- 0.01
mu <- 1/rate
muestra <- rexp(n, rate = rate)
media <- mean(muestra)
desv <- sd(muestra)

library(boot)
statistic <- function(data, i){
  remuestra <- data[i]
  c(mean(remuestra), var(remuestra)/length(remuestra))
}
B <- 2000
set.seed(1)

system.time(res.boot <- boot(muestra, statistic, R = B))
```

```
##    user  system elapsed 
##    0.03    0.00    0.03
```

```r
# system.time(res.boot <- boot(muestra, statistic, R = B, parallel = "snow"))
system.time(res.boot <- boot(muestra, statistic, R = B, parallel = "snow", cl = cl))
```

```
##    user  system elapsed 
##    0.03    0.00    0.03
```

### Estudio de simulación {#estudio-sim-boot}

Si se trata de un estudio más complejo, como por ejemplo un estudio de simulación en el que se emplea bootstrap, la recomendación sería tratar de paralelizar en el nivel superior para minimizar la sobrecarga debida a la comunicación entre nodos.

Por ejemplo, a continuación se realiza un estudio de simulación comparando las probabilidades de cobertura y las longitudes de los intervalos de confianza implementados en la función `boot.ci()`.


```r
t.ini <- proc.time()

nsim <- 500

getSimulation <- function(isim, B = 2000, n = 30, alfa = 0.1, mu = 100) {
    rate <- 1/mu # 0.01
    resnames <- c("Cobertura", "Longitud")
    # intervals <- c("Normal", "Percentil", "Percentil-t", "Percentil-t simetrizado")
    intervals <- c("Normal", "Basic", "Studentized", "Percentil", "BCa")
    names(intervals) <- c("normal","basic", "student", "percent", "bca")
    intervals <- intervals[1:4]
    resultados <- array(dim = c(length(resnames), length(intervals)))
    dimnames(resultados) <- list(resnames, intervals)
    # for (isim in 1:nsim) { # isim <- 1
    muestra <- rexp(n, rate = 0.01)
    media <- mean(muestra)
    desv <- sd(muestra)
    # boot()
    library(boot)
    statistic <- function(data, i){
      remuestra <- data[i]
      c(mean(remuestra), var(remuestra)/length(remuestra))
    }
    res.boot <- boot(muestra, statistic, R = B)
    res <- boot.ci(res.boot, conf = 1 - alfa)
    # Intervalos
    res <- sapply(res[names(intervals)], function(x) {
      l <- length(x)
      x[c(l-1, l)]
    })
    # resultados
    resultados[1, ] <- apply(res, 2,
                                   function(ic) (ic[1] < mu) && (mu < ic[2])) # Cobertura
    resultados[2, ] <- apply(res, 2, diff) # Longitud
    resultados
}

parallel::clusterSetRNGStream(cl)
result <- parLapply(cl, 1:nsim, getSimulation)
# stopCluster(cl)

# result
t.fin <- proc.time() - t.ini
print(t.fin)
```

```
##    user  system elapsed 
##    0.03    0.00    1.92
```

```r
resnames <- c("Cobertura", "Longitud")
intervals <- c("Normal", "Basic", "Studentized", "Percentil", "BCa")
names(intervals) <- c("normal","basic", "student", "percent", "bca")
intervals <- intervals[1:4]
resultados <- sapply(result, function(x) x)
dim(resultados) <- c(length(resnames), length(intervals), nsim)
dimnames(resultados) <- list(resnames, intervals, NULL)

res <- t(apply(resultados, c(1, 2), mean))
res
```

```
##             Cobertura Longitud
## Normal          0.858   57.788
## Basic           0.854   57.739
## Studentized     0.904   67.003
## Percentil       0.860   57.739
```

```r
knitr::kable(res, digits = 3)
```



|            | Cobertura| Longitud|
|:-----------|---------:|--------:|
|Normal      |     0.858|   57.788|
|Basic       |     0.854|   57.739|
|Studentized |     0.904|   67.003|
|Percentil   |     0.860|   57.739|


El último paso es finalizar el cluster:

```r
stopCluster(cl)
```


