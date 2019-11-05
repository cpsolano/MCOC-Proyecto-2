# MCOC-Proyecto-2

En este proyecto se busca poder simular el transporte de sedimento minero.

### Integrantes
 - Piedad Bull: https://github.com/piedadbull/MCOC-Proyecto-2-
 - Matias Echagüe: https://github.com/meechaguep/MCOC-Proyecto-2
 - Pedro Naretto: https://github.com/PedroNaretto/MCOC-Proyecto-2
 - Catalina Solano: https://github.com/cpsolano

## Entrega 4

En este proyecto se busca poder simular el transporte de sedimento minero de fondo de ríos. Esto se realizará en base a métodos langrangianos, es decir, que estudia el movimiento de cada partícula de forma individual.

Se tomarán como supuestos las siguientes condiciones:

- Diámetros de la partícula (d)= 0.15 mm
- Densidad de la partícula (rho_particula) = 2650 kg/(m^3)
- Densidad del agua (rho_agua) = 1000 kg/(m^3)
- Constante de drag (Cd) = 0.47   
- Constante de lifting (Cl) = 0.2    
- Constante de peso (Cm) = 0.5 
- Velocidad del flujo en x (vfx) = 5.0 m/s   
- Velocidad del flujo en y (vfy) = 0.1 m/s 


El intervalo de tiempo en los que grafica la posición de la partícula fue de (dt) = 0.001 s, con un tiempo máximo de simulación (tmax) de 1 s. Esto fué igual para las tres pruebas mostradas en los resultados.

### Validación

Al comparar los datos que se obtuvieron a partir del código con los del profesor, se puede ver una similitud para 2 particulas, con una diferencia de un decimal, lo cual se puede justificar por temas de unidades, pero se logra hacer una validación por los resultados obtenidos que se mueven en un rango parecido haciendo la excepción del decimal.

![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/2%20p.png)
![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/particle_positions%202.png)


Por otro lado, no se logra hacer una validación para una mayor cantidad de particular con respecto al movimiento que se espera entre ellas, pero se logran valores parecidos igual. El tener un movimiento distinto se puede deber a una modelación distinta del choque entre particulas.

![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/5%20p.png)

![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/10%20p.png)

### Comparación para 20 particulas

![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/20.1%20p.png)

![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/particle_positions%201.png)

### Computador Catalina

El computador utilizado por mi parte es un HP de 15.6" Intel Core i5 8va generación, 8 GB de Memoria RAM.

### Resultados

El código se corrió con tres cantidades distintas de partículas.
Primero se hizo la simulación con 4 partículas, con el cual se demoró 14.33 segundo en correr. Luego se volvió a realizar la simulación con 11 partículas, demorándose un tiempo de 120.06 segundos. Finalmente, el número de partículas fue 20, con el que demoró un tiempo mayor, el cual fue 823.03 segundos.

A continuación se presentan los gráficos para 4,11 y 20 partículas.

### 4 partículas
![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/4%20p.png)

### 11 partículas
![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/11%20p.png)

### 20 partículas

![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/20%20p.png)


## Entrega 6

En esta nueva entrega se solicitó optimizar el código de tal forma que al probar con distinta cantidad de partículas, se pudiera demostrar un comportamiento lineal entre el tiempo y la cantidad de partícula.

### Gráficos de movimientos para distinta cantidad de particulas


### Resultados Computador Catalina
El computador utilizado por mi parte es un HP de 15.6" Intel Core i5 8va generación, 8 GB de Memoria RAM, se mantuvo enchufado.


### Resultados Computador Matias


### Resultados Computador Pedro
Se utilizo un computador Acer Aspire con Intel Core i5 de 2.3 GHz con 8 GB de memoria RAM, el cual se mantuvo enchufado.

### Resultados Computador Piedad

El computador que se uso para esta entrega es un MacBook Pro de 13-inch con un procesador 2 GHz Intel Core i5. Tiene 8 GB de memoria ram y un almacenamiento flash de 251 GB. Además, tiene una tarjeta gráfica Intel Iris Graphic 540 de 1536 MB. El computador se mantuvo desenchufado.

![al text](https://github.com/piedadbull/MCOC-Proyecto-2-/blob/master/Entrega6/GraficoTiempoVsNparticulas.png)






