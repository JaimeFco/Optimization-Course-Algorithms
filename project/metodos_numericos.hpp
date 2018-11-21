//
//  metodos_numericos.hpp
//  proyecto00
//
//  Created by Jaime Francisco Aguayo on 10/09/17.
//  Copyright © 2017 Jaime Francisco Aguayo. All rights reserved.
//

#ifndef metodos_numericos_hpp
#define metodos_numericos_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <math.h>
#include <cmath>
#include <utility>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>

#define iterate(i, n) for(i = 0; i<n; ++i)

/** Funciones Elementales */

// Crea matriz de doble presición.
double **createMatrix(int nr, int nc);

// Libera memoria de matriz
void freeMatrix(double ** mat);

// Crea vector
double * createVector(int n);

// Libera memoria del vector
void freeVect(double * vect);

// Multiplica matriz por vector, Ma = B
void prodMatrizVector(double **M, double *A, double *B, int n);
void prodMatrizVector2(double **M, double *A, double *B, int n, int m);  // Longitud de vector es m

//Producto matricial entre P y Q que se guarda en R
void productoMatriz(double **P, double **Q, double **R, int n); // Modifica la matriz
void producto_matricial(double** P, double** Q, double** R, int n);
void producto_matricial2(double** P, double** Q, double** R, int n, int m, int p); // P de dim nxm, Q de dim mxp, R de dim nxp

//Obtiene norma infinito de la matriz A
double matrixInfNorm(double **A, int nr, int nc);

// Calcula la norma de la diferencia entre dos vectores u y v.
double vectorDifNorm(double *u, double *v, int n);

// Calcula la norma euclídea del vector
double vectorNorm(double *v, int n);

// Calcula la suma de dos vectores
void vectorSum(double *v1, double *v2, double *r, int n);

// Calcula la resta de dos vectores
void vectorDif(double *v1, double *v2, double *r, int n);

// Calcula el producto v1Av2, regresa por la izquierda el escalar resultante
double vectorMatrizVector(double *v1, double **A, double *v2, int n);

// Regresa el producto punto de A y B
double productoPunto(double *A, double *B, int n);

// Lee una matriz de archivo y la guarda con dimensiones nrxnc. Regresa 0 no se pudo almacenar.
double ** readMatrix(const char file[], int &nr, int &nc);

// Lee un vector de un archivo binario
double * readArray(const char file[], int &n);

// Guarda la matriz M en un archivo binario con el nombre especificado
void writeMatrix(const char file[], double ** M, int nr, int nc);

// Imprime la matriz
void printMatrix(double **M, int nr, int nc);

// Limpia la matriz llenandola de 0s
void clearMatrix(double **A, int nr, int nc);

// Limpia el vector llenandolo de 0s
void clearVector(double *v, int n);

// Copia la matriz A a la matriz B
inline void copyMatrix(double **A, double **B, int nr, int nc);

// Copia el vector a al vector b
void copyVector(double *a, double *b, int n);

// Crea la matriz de aa^t
void prodVectorVector(double *A, double *B, double ** C, int n);

/** Para soluciones de ecuaciones lineales */

// Resuelve Dx = B donde D es una matriz diagonal, regresa 0 si se resolvió sin ningún problema
// en otro caso regresa algo distinto de 0
int solveDiagonalMatrix(double *D, double *B, double *X, int n);

// Resuelve Mat*X = B donde Mat es una matriz triangular inferior, regrea 0 si no hubo problemas.
// algo distinto de 0 si hubo algún error al tratar de resolver la operación
int solveLowerMatrix(double **Mat, int nr, int nc, double *B, double *X);

// Resuelve Mat*X = B donde Mat es una matriz triangular superior,
// regresa 0 si no hubo problemas en la división
int solveUpperMatrix(double **Mat, int nr, int nc, double *B, double *X);

// Factoriza la matriz A en la multiplicación de L*U con el método Krout,
// OJO: No invade memoria de A por lo que desperdicia memoria
int fatorizacionLU(double **A, double **L, double **U, int n);

//Calcula el error que hay entre A y LU por memedio de la norma infinito
double errorLU(double **A, double **L, double **U, int size);

// Resuelve la ecuación MX = v, M matriz cuyas diagonales,
// de abajo hacia arriba son a, b y c. Regresa 0 si se pudo resolver y 1 en otro caso
int solveTridiagonalMatrix2(double *a, double *b, double *c, int n, double *v, double *X);

// Resuelve la ecuación MX = v, M matriz cuyas diagonales,
// de abajo a arriba son a, b y c. Inicialmente X = v. Regresa 0 si se pudo resolver y 1 en otro caso
int solveTridiagonalMatrix(double *a, double *b, double *c, int n, double *d, double *x);

// Trata de resolver Ax = B por distintos métodos, regresa 0 si lo logró
int solve(double ** A, int n, double *X, double *B);



/** Para el cálculo de eigenvalores y/o eigenvectores */

// Metodo de la potencia para encontrar en eigenvalor mas grande y su respectivo eigenvector
// Matriz A cuadrada, vector V para guardar el eigenvector, tol tolerancia para parar y
// MAX es el numero maximo de iteraciones. Al final regresa el eigenvalor encontrado
double metodoDeLaPotencia(double **A, int N, double *v = NULL, double tol = 0.001, double MAX = 1000000);

// Funcion que recibe una matriz A, un vector inicial X, un valor delta donde realizar la busqueda,
// un valor de tolerancia para el resultado y la direccion donde guardar el eigen vector encontrado.
// Al final deja el vector propio en X, y regresa el número de iteraciones que realizó.
double metodoPotenciaInversa(double **A, int n, double *X, double delta, double tol, double & eigenvalue, int N = 10000);

// Usa el método de Jacobi que usa rotaciones de Givens para regresa los eigenvalores y eigenvectores
// de la tariz A. Los eigenvectores se encuentran ordenados en la matriz eigenvectors y eigenvalores
// se encuentran en el mismo orden en eigenvalues. Se debe dar tolerancia y max de iteraciones
double metodoJacobi(double **A, int n, double ** &eigenvectors, double *eigenvalues, double tol, int N = 20000);

// Regresa el valor más grande en valor absoluto de la matriz mat fuera de la diagonal.
void greatValue(double **mat, int n, int &i, int &j);

// Rota matrices simétricas donde c = cos y s = sin del ángulo a girar con pivote i,j
// Se usa con el metodo de Jacobi
void rot(double **mat, int n, int i, int j, double s, double c);

// Metodo del gradiente conjugado, resuelve el sistema Ax = b para matrices simetricas.
// Es necesario dar un valor de tolerancia y un número maximo de iteraciones, regresa el error
double metodoGradienteConjugado(double **A, int n, double *x, double *b, double tol, int N = 0);


/** Para interpolación */

// Guarda en A la matriz de vandermonde de tamaño n+1 producida por x
void vandermondeMatrix(double *x, double **A, int n);

// Genera en C los coeficientes del polinomio que interpola los puntos x, y
void polinomioInterpolador(double *x, double *y, double *c, int n);

// Evalua el polinomio con coeficientes en C
double evaluaPolinomio(double x, double *c, int n);

// Crea el polinomio de lagrange Ln,i y lo evalua en x, regresa la evaluación
// OJO: n NO es el tamaño del arreglo, el arreglo va de 0 a n
double evaluaLagrangePol(double *X, int n, int i, double x);

// Interpola el punto x por medio de los polinomios de Lagrange
// OJO: n aquí SI es el tamaño de los vectores
double interpolaLagrange(double x, double * X, double * Y, int n);

// Esta funcion devuelve el producto de j=0 a i-1 de (x-X[j]).
// Es el polinomio base para la interpolacion de Newton.
double NewtonBasisPolynomial(double *X, int i, double x);

// Funcion que regresa el arreglo a, las diferencias divididas F[x0] , F[x0,x1], ... para usarse en el
// polinomio interpolador de Newton. n es el tamaño de los arreglos, sus indices son 0, 1, ..., n-1
void NewtonDividedDifference(double *X, double *Y, double *a, int n);

// Funcion que devuelve el valor del punto x en el polinomio interpolador de Newton
// n es el tamaño de los vectores, por lo que sus indices van de 0 a n-1
double interpolaNewton(double x, double *X, double *Y, int n);

// Funcion que devuelve en y las evaluaciones de los puntos x en el polinomio
// interpolador de Newton generado por X y Y, n es el tamaño de los vectores X y Y
// mientras que N es el tamaño de los vectores x y y.
void interpolaNewton(double *x, double *y, int N, double *X, double *Y, int n);

// Función que interpola un conjunto de puntos por medio de splines cúbicos.
// La entrada contiene un arreglo con las coordenadas de los puntos en x y otro con los puntos en y
// además del tamaño n de los arreglos. El arreglo X contiene un conjutno de puntos a interpolar
// y los valores de la interpolaciones se guardan en Y, estos dos arreglo deben tener tamaño m
// OJO: Se obvia que las X's estan en orden ascendente
void interpolaSplineCubico(double *x, double *y, int n, double *X, double *Y, int m);

// Interpola el punto x usando spline cúbico
double interpolaSplineCubico(double X, double *x, double *y, int n);

// Funcion que se usa en interpolaSplineCubico para calcular los coeficientes
void coeficientes_splines_cubicos(double *x, double *y, int n, double *M);

// Metodo de aproximación de elementos finitos
void interpolacionElementosFinitos(double *X, double *Y, int N, double *x, double *y, int n, int m, double lambda);

// Calcula coeficientes para metodo de elementos finitos
void calcula_aproximacion(double *phi, int m, double *x, double *y, int n, double lambda);


// Método de integración de Romberg, calcula la integral de func de a a b, con partición 2^n
double rombergIntegration(double (*func)(double), double a, double b, int n);

// Auxiliar que calcula el siguiente término de la expresión recursiva de la regla del trapecio
double nextTrapeziumRule(double previous, int n, double a, double b, double (*func)(double));


void saveTxt(double *a1, double *a2, int n, const char *cfile);


// Calcula la inversa de una matriz via Cholesky. Regresa 1 si no se pudo factorizar.
bool inversa(double ** M, double ** Mi, int n);

// Traspone una matriz M de nxn
void traspose(double ** M, int n);

// Traspone una matriz en otra
void transpose2(double ** O, double ** D, int n, int m);

// Calcula la Factorización de Cholesky ( M = LL^t )
double ** cholesky(double **M, int n);

void scaleVector(double s, double * v, int n);

bool inversaCholesky(double ** M, double ** Mi, int n);

void squareMatrixSum(double ** A, double ** B, double ** R, int n);



// GSL

// Lectura de la matriz en el archivo cfile.
// Devuelve NULL si no se pudo abrir el archivo.
gsl_matrix *readMatrix(char *cfile);

// Lectura del vector en el archivo cfile.
// Devuelve NULL si no se pudo abrir el archivo.
gsl_vector *readVector(char *cfile);

// Escritura de la matriz mat en el archivo cfile.
// Devuelve 0 en caso de exito y 1 si no.
int writeMatrix(gsl_matrix *mat, char *cfile);

// Imprime en la consola las entradas de la matriz
void printMatrixGSL(gsl_matrix *mat);

// Multiplica por dft a los elementos de la matriz
void scaleMatrixGSL(gsl_matrix *mat, double dft);

void scaleMatrix(double s, double **M, int r, int c);

// Devuelve un doble apuntador para poder trabajar con las
// entradas de las matriz mat con un arreglo bidimensional
double **gslmat2array2d(gsl_matrix *mat);

// Imprime en la consola las entradas del vector
void printVectorGSL(gsl_vector *vec);

// Hace que el vector GSL haga referencia a el arreglo v
// OJO: No hace copia del arreglo, solo lo referencía.
void initGslVector(gsl_vector      *vec, double *v, int dim);

// Inicializa gsl_matrix con los datos de un arreglo bidimensional
// OJO: No hace copia el arreglo 2-D, solo lo referencía.
void initGslMatrix(gsl_matrix  *mat, double **M, int nr, int nc);

// Copia el contenido de el vector de la GSL v al arreglo A, siempre y cuando A tenga memoria
void copyGsl_vector2Array(gsl_vector  *v, double *A);


#endif /* metodos_numericos_hpp */
