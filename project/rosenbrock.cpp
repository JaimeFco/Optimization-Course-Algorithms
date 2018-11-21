//
//  main.cpp
//  Examen 2 práctico
//
//  Created by Jaime Francisco Aguayo on 20/05/18.
//  Copyright © 2018 Jaime Francisco Aguayo. All rights reserved.
//

#include <iostream>
#include "optimizacion.hpp"

#define INTRO "\n"
#define PI 3.141589

int N = 100;
double f(double * x)
{
    double r = 0.0;
    for (int i=0; i<N-1; i++)
        r += 100*pow(x[i+1] - x[i]*x[i], 2) + ((1.0-x[i])*(1.0-x[i]));
    
    return r;
}

void g(double* x, double* gx)
{
    gx[0] = -400.0*x[0]*(x[1] - (x[0]*x[0])) - 2.0*(1 - x[0]);
    for (int i=1; i<N-1; i++) {
        gx[i] = -400.0*x[i]*(x[i+1] - (x[i]*x[i])) - 2.0*(1 - x[i]);
        gx[i] += 200.0*(x[i] - x[i-1]*x[i-1]);
    }
    gx[N-1] = 200.0*(x[N-1] - x[N-2]*x[N-2]);
}


double fwood(double * x)
{
    double r = 100*(x[0]*x[0] - x[1])*(x[0]*x[0] - x[1]) + (x[0] - 1.0)*(x[0] - 1.0) + (x[2] - 1.0)*(x[2] - 1.0) +   90.0*(x[2]*x[2] - x[3])*(x[2]*x[2] - x[3]) + 10.1 *(pow((x[1] - 1.0), 2) + pow((x[3] - 1.0), 2)) + 19.8*(x[1]-1.0)*(x[3] - 1);
    return r;
}

void gwood(double* x, double* gx)
{
    gx[0] = 2.0*(200.0*x[0]*(x[0]*x[0]-x[1]) + x[0] - 1.0);
    gx[1] = 19.8*x[3] - 200.0*x[0]*x[0] + 220.2*x[1] - 40.0;
    gx[2] = 2.0*(180.0*x[2]*(x[2]*x[2]-x[3]) + x[2] - 1.0);
    gx[3] = 200.2*x[3] + 19.8*x[1] - 180.0*x[2]*x[2] - 40.0;
}


int main(int argc, const char * argv[])
{
    double * x = createVector(N);
    double * gx = createVector(N);
    double * lb = createVector(N);
    double * ub = createVector(N);
    for (int i=0; i<N; i++) {
        lb[i] = -2.0;
        ub[i] = 2.0;
    }
    
    cso(18, f, x, lb, ub, N, 10000);
    
    g(x, gx);
    std::cout<<"Norma del gradiente: "<<vectorNorm(gx, N)<<INTRO;
    for (int i=0; i<N; i++) {
        std::cout<<x[i]<<INTRO;
    }
    
    free(x);
    free(gx);
    free(lb);
    free(ub);
    return 0;
}
