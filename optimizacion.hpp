//
//  optimizacion.hpp
//  proyecto00
//
//  Created by Jaime Francisco Aguayo on 08/03/18.
//  Copyright © 2018 Jaime Francisco Aguayo. All rights reserved.
//

#ifndef optimizacion_hpp
#define optimizacion_hpp

#include <stdio.h>
#include "metodos_numericos.hpp"


/** Optimización */


double * lineSearchMethod(int n, double (*f)(double*), void (*g)(double*, double*), void (*H)(double *, double**), double * x0, const char type[] = NULL, int M = 2048, double tolg = __FLT_EPSILON__, double tolx = __FLT_EPSILON__, double tolf= __FLT_EPSILON__, double alpha = 0.1);


// Función que recibe:
// xnext: vector donde vamos a guardar el nuevo valor x
// x: vector con el valor anteior (x actual)
// n: tamaño de los vectores
// H: Función que calcula el hessiano evaluado en su parametro
// G: Funcion que calcula el gradiente
//
// Devuelve:
// Paso para x en xnext
// un bool 0 = se pudo usar hessiano, 1 = se tuvo que usar un paso fijo con valor de 0.5
bool NewtonStep(double* xnext, double* x, int n, void (*H)(double *, double**), void (*G)(double*, double*));

// Función que recibe:
// xnext: vector donde vamos a guardar el nuevo valor x
// x: vector con el valor anteior (x actual)
// n: tamaño de los vectores
// H: Función que calcula el hessiano evaluado en su parametro
// G: Funcion que calcula el gradiente
// H0: Matriz Heesiana inicial
// Devuelve:
// Paso para x en xnext
// un bool 0 = se pudo usar hessiano, 1 = se tuvo que usar un paso fijo con valor de 0.5
bool aproxHessStep(double* xnext, double* x, int n, void (*G)(double*, double*), double (*f)(double *), double & alpha);


double line_search_steepest_descent(int n, double (*f)(double*), void (*g)(double*, double*), double * x0, double * xoptimo, const char type[] = "backtracking", int M = 1024, double alpha0 = 1.0, double tolg = DBL_EPSILON, double tolx = DBL_EPSILON, double tolf= DBL_EPSILON);

double backtrackingStep(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*), int n);

double quadraticInterpolationStep(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*), int n);

double cubicInterpolationStep(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*), int n);

// Método Dogleg de optimización basado en región de confianza
// Inputs:
// n: tamaño de los vectores
// f: apuntador a función objetivo
// g: spuntador a función gradiente
// x0 punto inicial
// xoptimo: apuntador vector donde se guaradará el punto óptimo
// Dmax: Radio máximo de la región de confianza
// M: número máximo de iteraciones
// tolg, tolx y tolf: valores de tolerancia para gradiente, diferencia entre x's y diferencia entre valores de f respectivamente
double dogleg(int n, double (*f)(double*), void (*g)(double*, double*), void (*H)(double*, double **), double * x0, double * xoptimo, double Dmax, int M, double tolg, double told);

// Método del gradiente conjugado lineal que converge en a lo más n pasos
// Input:
// n: tamaño de los vectores
// f: apuntador a función objetivo
// g: spuntador a función gradiente
// Q: Matriz del modelo cuadrático
// x0: punto inicial
// xoptimo: apuntador vector donde se guaradará el punto óptimo
// tolg, tolx y tolf: valores de tolerancia para gradiente, diferencia entre x's y diferencia entre valores de f respectivamente
// Output:
// Regresa por la izquierda el valor de la función en el óptimo
double conjGrad(int n, double (*f)(double*), void (*g)(double*, double*), double ** Q, double * x0, double * xoptimo, double tolg = DBL_EPSILON, double tolx = DBL_EPSILON, double tolf= DBL_EPSILON);

double conjGradNL(int n, double (*f)(double*), void (*g)(double*, double*), const char metodo[], double * x0, double * xoptimo, double tolg = DBL_EPSILON, double tolx = DBL_EPSILON, double tolf = DBL_EPSILON, double M = 5000);


// Método que resuelve un sitema de ecuaciones no lineales F(x)=0 usando el método de Newton.
// Input:
// x0: vector con el punto inicial
// x: vector donde se guardará el punto x*
// n: dimensión d los vectores
// F: función F(x) del sistema
// J: función que calculo el Jacobiano de F
// tol: tolerancia para el paro
// M: número máximo de iteraciones
void solveSENLNewton(double *x0, double *x, int n, void (*F)(double*, double*), void (*J)(double*, double**), double tol, int M);

void solveSENLBroyden2(double *x0, double *x, int n, void (*F)(double*, double*), double ** J0, double tol, int M);
void solveSENLBroyden(double *x0, double *x, int n, void (*F)(double*, double*), double ** J0, double tol, int M);


// Método cuasi-Newton BFGS para optimizar funciones
// Input:
// x0: vector con el punto inicial
// x: vector donde se guardará el punto x*
// n: dimensión del problema
// f: apuntador a la función a minimizar
// H0: Aproximación de la inverza de la matriz hessiana en x0
// tol: Tolerancia para la norma del gradiente
// M: máximo de iteraciones
void BFGS(double *x0, double *x, int n, double (*f)(double*), void (*g)(double*, double*), double ** H0, double tol, int M);

// Función que calcula la matriz hessiana en el punto x de la función f con un tamaño h
void hessianFiniteDif(double **H, double (*f)(double*), double *x, int n, double h);

// Método de Levenberg-Marquardt para resolver mínimos cuadrados no lineales
// Input:
// x0: vector con el punto inicial
// x: vector donde se guardará el punto x*
// n: dimensión del dominio
// m: dimensión del codominio
// R: función R(x) del sistema R: R^n -> R^m
// J: función que calculo el Jacobiano de R
// tol: tolerancia para el paro
// M: número máximo de iteraciones
void LM_method(double *x0, double *x, int n, int m, void (*R)(double*, double*), void (*J)(double*, double**), double nu, double tol, int M);

void quadraticPenaltyMethod(double (*Q)(double*, double), void (*gradQ)(double*, double*, double), void (*g)(double*, double*), void (*h)(double*, double*), int n, int m, int r, double *x0, double tol, double N);

void quadraticPenaltyFunction(double (*f)(double*), int n, void (*g)(double*, double*), int m, void (*h)(double*, double*), int r, double *x, double *gradx);

double line_search_for_QPM(int n, double (*f)(double*, double mu), void (*g)(double*, double*, double), double mu, double * x0, double * xoptimo, int M, double alpha0, double tolg, double tolx, double tolf);

double backtrackingStepForQPM(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*, double), int n, double mu);

#endif /* optimizacion_hpp */
