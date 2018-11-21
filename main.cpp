//
//  main.cpp
//  tarea12 Problema 3
//
//  Created by Jaime Francisco Aguayo on 20/05/18.
//  Copyright © 2018 Jaime Francisco Aguayo. All rights reserved.
//

#include <iostream>
#include "optimizacion.hpp"

#define INTRO "\n"
#define PI 3.141589

int n=2, m=1, r=4;

double f(double *x)
{
    double fx;
    fx = (x[0] - 6.0)*(x[0] - 6.0) + (x[1] - 7.0)*(x[1] - 7.0);
    return fx;
}
void g(double *x, double *gx)
{
    gx[0] = 0.0;
}
void h(double *x, double *hx)
{
    hx[0] = -3.0*x[0] - 2.0*x[1] + 6.0;
    hx[1] = -x[0] + x[1] - 3.0;
    hx[2] = x[0] + x[1] - 7.0;
    hx[3] = (2.0/3.0) * x[0] - x[1] - 4.0/3.0;
}

double Q(double *x, double mu)
{
    double * gx = createVector(m);
    double Qx;
    double * hx = createVector(r);
    
    g(x, gx);
    h(x, hx);
    for (int i=0; i<r; i++) {
        hx[i] = (hx[i]>0.0) ? (-hx[i]) : 0.0;
    }
    
    Qx = f(x) + (mu*0.5) * ( productoPunto(gx, gx, m) + productoPunto(hx, hx, r) );
    
    free(gx);
    free(hx);
    return Qx;
}

void gradQ(double *x, double *gradx, double mu)
{
    double * xp = createVector(n);
    copyVector(x, xp, n);
    
    for (int i=0; i<n; i++) {
        xp[i] += 0.0001;
        gradx[i] = (Q(xp, mu) - Q(x, mu)) / 0.0001;
        xp[i] -= 0.0001;
    }
}

int main(int argc, const char * argv[])
{
    // Punto inicial
    double * x0 = createVector(n);
    double * hx = createVector(r);
    
    x0[0] = -15.0;
    x0[1] = -10.0;
    
    quadraticPenaltyMethod(Q, gradQ, g, h, n, m, r, x0, 0.005, 200);
    
    std::cout<<"Punto solucion: ("<<x0[0]<<", "<<x0[1]<<")\n"
             <<"f(x*) = "<<f(x0)<<"\n";
    for (int i=0; i<r; i++) {
        std::cout<<hx[i]<<"\n";
    }
    
    free(x0);
    free(hx);
    return 0;
}
