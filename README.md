# MCOC-Proyecto-2

En este proyecto se busca poder simular el transporte de sedimento minero.

### Integrantes
 - Piedad Bull: https://github.com/piedadbull
 - Matias Echagüe: https://github.com/meechaguep
 - Pedro Naretto: https://github.com/PedroNaretto
 - Catalina Solano: https://github.com/cpsolano

## Entrega 4

En este proyecto se busca poder simular el transporte de sedimento minero de fondo de ríos. Esto se realizará en base a métodos langrangianos, es decir, que estudia el movimiento de cada partícula de forma individual.

Se tomarán como supuestos las siguientes condiciones:

- Diámetros de la partícula (d)= 15 mm
- Densidad de la partícula (rho_particula) = 2650 kg/(m^3)
- Densidad del agua (rho_agua) = 1000 kg/(m^3)
- Constante de drag (Cd) = 0.47   
- Constante de lifting (Cl) = 0.2    
- Constante de peso (Cm) = 0.5 
- Velocidad del flujo en x (vfx) = 5.0 m/s   
- Velocidad del flujo en y (vfy) = 0.1 m/s 


El intervalo de tiempo en los que grafica la posición de la partícula fue de (dt) = 0.001 s, con un tiempo máximo de simulación (tmax) de 1 s. Esto fué igual para las tres pruebas mostradas en los resultados.

### Validación



### Computador Catalina

El computador utilizado por mi parte es un HP de 15.6" Intel Core i5 8va generación, 8 GB de Memoria RAM.

### Resultados

El código se corrió con tres cantidades distintas de partículas.
Primero se hizo la simulación con 4 partículas, con el cual se demoró 1.9 segundo en correr. Luego se volvió a realizar la simulación con 11 partículas, demorándose un tiempo de 20.1 segundos. Finalmente, el número de partículas fue 20, con el que demoró un tiempo mayor, el cual fue 121.77 segundos.

A continuación se presentan los gráficos para 4,11 y 20 partículas.

### 4 partículas
![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/Gráfico%204%20p.png)

### 11 partículas
![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/Gráfico%2011%20p.png)

### 20 partículas
![al text](https://github.com/cpsolano/MCOC-Proyecto-2/blob/master/Gráficos/Gráfico%2020%20p.png)

