For the c++ codes:

To compile:
g++ -o rosenbrock rosenbrock.cpp metodos_numericos.cpp optimizacion.cpp -lgsl -lgslcblas -lm
g++ -o wood wood.cpp metodos_numericos.cpp optimizacion.cpp -lgsl -lgslcblas -lm

To run (linux / MacOS:
./ rosenbrock
./ wood

You need to install GNU Scientific Library (GSL).

These two programs run the CSO algorithm to find the optimal argument of Rosenbrock and Wood functions. They print the norm of the gradiente in each iteration.

For the Python codes:
python dropWave_cso.py

You need to install matplotlib and pandas.

The program runs an animation of the agents in each iteration.


The report of the project is the file ChickenSwarmOptimization.pdf and the slides I used in class are presentacion-CSO.pdf



Spanish:
En la presente carpeta se incluye:
ChickenSwarmOptimization.pdf - El reporte final del proyecto.
presentacion-CSO.pdf - La presentación usada en la expocisión de la entrega de proyecto.
dropWave_cso.py - Código de python que recrea una animación del algoritmo al aplicarlo a la función Drop-Wave de dimensión 2.
rosenbrock_cso.py - Código de python con la animación de CSO aplicado en la función de Rosenbrock.
Códigos - Una carpeta donde se incluyen los códigos implementados en C++

Para compilar los códigos de C++ basta con hacer:
g++ -o programa rosenbrock.cpp metodos_numericos.cpp optimizacion.cpp -lgsl -lgslcblas -lm

La función cso() se encuentra en optimizacion.cpp y recibe los siguientes parámetros:
n: Número de agentes a usar en el algoritmo
*function: Función a optimizar
*x: Arreglo de tipo double donde almacenar el punto óptimo
*lb: Arreglo con los límites inferirores para cada variable
*ub: Arreglo con los límites superiores para cada variable
dimension: Entero con el valor del número de variables de que depende la función objetivo
iteration: Número de iteraciones a realizar
G: Número de iteraciones a pasar antre cada actualización de las jerarquías (por defecto es 5)
FL: Double definido en el reporte (Por defecto es 0.5)
