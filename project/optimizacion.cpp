//
//  optimizacion.cpp
//  proyecto00
//
//  Created by Jaime Francisco Aguayo on 08/03/18.
//  Copyright © 2018 Jaime Francisco Aguayo. All rights reserved.
//

#include "optimizacion.hpp"



double * lineSearchMethod(int n, double (*f)(double *), void (*g)(double*, double*), void (*H)(double*, double **), double * x0, const char type[], int M, double tolg, double tolx, double tolf, double alpha)
{
    double * xprev = new double[n];
    double * xnext = new double[n];
    double fprev;
    double fnext;
    double norm=0.0;
    double norm2print = 0.0;
    
    // Vector gradiente
    double *vgrad = new double[n];
    double *aux = new double[n];
    // Lee mensaje y decide método
    std::string mensaje(type);
    char method;
    if(mensaje == std::string("StepFijo"))
        method = 'a';
    else if(mensaje == std::string("StepHess"))
        method = 'b';
    else if (mensaje == std::string("StepAprox"))
        method = 'c';
    
    else    // Por defecto tamaño fijo
        method = 'a';
    
    std::cout<<"Se optimizará con un alpha "<<type<<", en el punto inicial (";
    for (int i=0; i<n; i++) {
        std::cout<<x0[i]<<", ";
    }
    std::cout<<"\n";
    
    copyVector(x0, xprev, n);
    // Mientras sean iteraciones menores que M
    for(int k=0; k<M; k++) {
        // Valor de paso fijo
        if(method == 'a') {
            g(xprev, vgrad); // Calculas gradiente
            for(int i=0; i<n; i++) //Creas xk+1 en la direccion del gradiente con paso alpha
                xnext[i] = xprev[i] -alpha*vgrad[i];
        }
        else if (method == 'b') {   // Método del descenso del gradiente de Newton (uso del Hessiano)
            bool ind = NewtonStep(xnext, xprev, n, H, g);
            if(ind)
                std::cout<<"Matriz Hessiana no se pudo invertir"<<"\n";
        }
        else if (method == 'c') {   // Método de aproximacion al hessiano
            bool ind = aproxHessStep(xnext, xprev, n, g, f, alpha);
            if (ind) {
                std::cout<<"No se pudo usar método, bye "<<"\n";
            }
        }
        // Revisa condiciones de paro
        // Primer condición
        for (int j=0; j<n; j++)
            aux[j] = xprev[j] - xnext[j];
        norm = vectorNorm(aux, n);
        if(norm/(std::max(1.0, vectorNorm(xprev, n))) < tolx) {
            std::cout<<"Se cumplió tolerancia en x_k - x_(k+1)\n";
            break;
        }
        norm2print = norm;
        // Segunda condición
        fprev = f(xprev);
        fnext = f(xnext);
        
        norm = abs(fprev-fnext);
        if(norm/(std::max(1.0, abs(fnext))) < tolf) {
            std::cout<<"Se cumplió tolerancia en f_k - f_(k+1)\n";
            break;
        }
        
        // Tercera condición
        g(xnext, vgrad);
        norm = vectorNorm(vgrad, n);
        if(norm < tolg) {
            std::cout<<"Se cumplió tolerancia del valor del gradiente en x\n\n";
            break;
        }
        // Imprimir Iteración y punto
        std::cout<<k+1<<"   Difx = "<< norm2print
        <<", ||∇f (x_k)|| = "<<norm<<", f_k = "<<fnext<<"\n";
        
        // Actualizando el punto
        for (int i=0; i<n; i++)
            xprev[i] = xnext[i];
    }
    
    std::cout<<"Final: "<<"    (";
    for (int i=0; i<n; i++) {
        std::cout<<xnext[i]<<", ";
    }
    
    std::cout<<"), Difx = "<< norm2print
    <<", ||∇f (x_k)|| = "<<norm<<", f_k = "<<fnext<<"\n";
    delete [] aux;
    delete [] xprev;
    delete [] vgrad;
    
    return xnext;
}

bool NewtonStep(double* xnext, double* x, int n, void (*H)(double *, double**), void (*G)(double*, double*))
{
    double ** Hess = createMatrix(n, n);
    H(x, Hess);  // Evaluando el hessinano
    double * g = new double[n];
    G(x, g);                // Evalua el gradiente en x y lo guarga en g
    
    double ** HessI = createMatrix(n, n);
    bool ind = inversa(Hess, HessI, n);    // Obteniendo la inversa del hessiano
    
    if(ind == true) {  // Si hubo problemas al calcular inversa (i.e. Hess no era positiva definida)
        for (int i=0; i<n; i++) {
            xnext[i] = x[i] - 0.0005*g[i];   // Usa 0.5 del gradiente como paso
        }
    }
    else {
        
        prodMatrizVector(HessI, g, xnext, n);   // Crea siguiente valor xnext
        for (int i=0; i<n; i++) {
            xnext[i] = x[i] - xnext[i];
        }
    }
    
    // Libera memoria
    freeMatrix(Hess);
    freeMatrix(HessI);
    delete [] g;
    
    return ind;
}


bool aproxHessStep(double* xnext, double* x, int n, void (*G)(double*, double*), double (*f)(double *), double & alpha)
{
    double * gprev = new double[n];
    double * s = new double[n];
    double fprev = f(x);
    double fnext;
    
    G(x, gprev);        //Calculamos graduente de xprev
    for(int i=0; i<n; i++) //Creas xk+1 en la direccion del gradiente con paso alpha k-1
        s[i] = x[i] - alpha*gprev[i];
    
    fnext = f(s);
    // Actualizando alpha
    double gg = productoPunto(gprev, gprev, n);
    double d = fnext - fprev + alpha*gg;
    
    alpha = (gg*alpha*alpha)/(2*d);
    
    // Actualizando la x
    for (int i=0; i<n; i++) {
        xnext[i] = x[i] - alpha*gprev[i];
    }
    
    // Liberando memoria
    delete [] gprev;
    delete [] s;
    
    return 0;
}


double line_search_steepest_descent(int n, double (*f)(double*), void (*g)(double*, double*), double * x0, double * xoptimo, const char type[], int M, double alpha0, double tolg, double tolx, double tolf) {
    double * xnext = createVector(n);
    double * xprev = createVector(n);
    double * vgrad = createVector(n);
    double * aux = createVector(n);
    double * dir = createVector(n);
    double alpha;
    std::string mensaje(type);
    
    std::cout<<"Se actualizará el paso con método "<<mensaje<<"\n";
    
    copyVector(x0, xprev, n);
    
    // Itera mientras no se cumplan los criterios de paro
    for (int k=0; k<M; k++) {
        std::cout<<k<<": ";
        alpha = alpha0;            // Alpha inicial
        g(xprev, vgrad);        // Obtenemos gradiente para dirección
        // Obtenemos la dirección de máximo descenso
        copyVector(vgrad, dir, n);
        scaleVector(-1.0, dir, n);
        
        // Seleccionamos método para obtener alpha
        if(mensaje == std::string("backtracking")) {             // Si quieren backtraking
            alpha = backtrackingStep(alpha, xprev, dir, vgrad, f, n);// Se busca alpha
        }
        else if (mensaje == std::string("interpolacion2")) {    // Si quieren interpolación cuadrática
            alpha = quadraticInterpolationStep(alpha, xprev, dir, vgrad, f, n);
        }
        else if (mensaje == std::string("interpolacion3")) {    // Si quieren interpolación cúbica
            alpha = cubicInterpolationStep(alpha, xprev, dir, vgrad, f, n);
        }
        // Actualizamos x_k con el alpha encontrado
        for (int j=0; j<n; j++)                             // Se obtiene el x_{k+1}
            xnext[j] = xprev[j] + alpha * dir[j];
        
            /* Condiciones de paro */
        // Primera condición
        for (int j=0; j<n; j++)
            aux[j] = xprev[j] - xnext[j];
        double norm = vectorNorm(aux, n);
//        std::cout<<norm<<"  ";
        if(norm/(std::max(1.0, vectorNorm(xprev, n))) < tolx) {
            std::cout<<"Se cumplió tolerancia en x_k - x_(k+1)\n";
            break;
        }
        
        // Segunda condición
        double fprev = f(xprev);
        double fnext = f(xnext);
        
        norm = abs(fprev-fnext);
//        std::cout<<norm<<"  ";
        if(norm/(std::max(1.0, abs(fnext))) < tolf) {
            std::cout<<"Se cumplió tolerancia en f_k - f_(k+1)\n";
            break;
        }
        
        // Tercera condición
        g(xnext, vgrad);
        norm = vectorNorm(vgrad, n);
        std::cout<<norm<<"  ";
        if(norm < tolg) {
            std::cout<<"Se cumplió tolerancia del valor del gradiente en x\n\n";
            break;
        }
        //std::cout<<f(xnext);
        std::cout<<"\n";
        
        copyVector(xnext, xprev, n);        // Siguiente paso
    }
    copyVector(xnext, xoptimo, n);      // Copia el x*
    
    free(xprev);
    free(xnext);
    free(vgrad);
    free(aux);
    free(dir);
    return f(xoptimo);                  // Regresa evaluación de x* en f
}


double backtrackingStep(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*), int n)
{
    double c = 0.0001;
    double p = 0.27;
    
    double alpha = alpha0;
    
    double fxa;
    double fx;
    double m;
    double * xnext = createVector(n);
    
    m = c * alpha * productoPunto(gk, dk, n);
    fx = f(xk);
    for (int i=0; i<n; i++) {
        xnext[i] = xk[i] + alpha * dk[i];
    }
    fxa = f(xnext);
    
    while (fxa > fx + m) {
        
        alpha = p*alpha;
        m = c * alpha * productoPunto(gk, dk, n);
        for (int i=0; i<n; i++) {
            xnext[i] = xk[i] + alpha * dk[i];
        }
        fxa = f(xnext);
    }
    free(xnext);
    return alpha;
}

double quadraticInterpolationStep(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*), int n)
{
    double c1 = 0.0001;
    double alpha1;
    
    double phi0 = f(xk);
    double phip0 = productoPunto(gk, dk, n);
    
    double * xnext = createVector(n);
    double phi_a0, phi_a1;
    
    for (int i=0; i<n; i++) {
        xnext[i] = xk[i] + alpha0 * dk[i];
    }
    phi_a0 = f(xnext);
    
    // Si cumple con Armijo
    if(phi_a0 <= phi0 + c1 * alpha0 * phip0) {
        free(xnext);
        return alpha0;
    }
    
    // Calcula alpha_1
    alpha1 = (-alpha0*alpha0 * phip0) / (2.0*(phi_a0 - phip0*alpha0 - phi0));
    
    // Para la siguiente alpha
    for (int i=0; i<n; i++) {
        xnext[i] = xk[i] + alpha1 * dk[i];
    }
    phi_a1 = f(xnext);
    // Hasta cumplir Armijo
    while (phi_a1 > phi0 + c1 * alpha1 * phip0) {
        alpha0 = alpha1;  // alpha_0 = alpha_1
        phi_a0 = phi_a1;
        
        // Calcula siguiente alpha
        alpha1 = (-alpha0*alpha0 * phip0) / (2.0*(phi_a0 - phip0*alpha0 - phi0));
        
        // Para la siguiente alpha, calcula phi(alpha)
        for (int i=0; i<n; i++) {
            xnext[i] = xk[i] + alpha1 * dk[i];
        }
        phi_a1 = f(xnext);
    }
    
    free(xnext);
    return alpha1+c1;
}


double cubicInterpolationStep(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*), int n)
{
    double c1 = 0.0001;
    double alpha1, alpha2;
    
    double phi0 = f(xk);
    double phip0 = productoPunto(gk, dk, n);
    
    double * xnext = createVector(n);
    double phi_a0, phi_a1, phi_a2;
    
    for (int i=0; i<n; i++) {
        xnext[i] = xk[i] + alpha0 * dk[i];
    }
    phi_a0 = f(xnext);
    
    // Si cumple con Armijo
    if(phi_a0 <= phi0 + c1 * alpha0 * phip0) {
        free(xnext);
        return alpha0;
    }
    
    // Calcula alpha_1
    alpha1 = (-alpha0*alpha0 * phip0) / (2.0*(phi_a0 - phip0*alpha0 - phi0));
    
    // Para la siguiente alpha
    for (int i=0; i<n; i++) {
        xnext[i] = xk[i] + alpha1 * dk[i];
    }
    phi_a1 = f(xnext);
    
    if(phi_a1 <= phi0 + c1 * alpha1 * phip0) {
        free(xnext);
        return alpha1;
    }
    
    double a, b, den;
    
    den = ( alpha0*alpha0 * alpha1*alpha1 * (alpha1 - alpha0));
    
    a = (alpha0*alpha0 * (phi_a1 - phip0*alpha1 - phi0) - alpha1*alpha1 * (phi_a0 - phip0*alpha0 - phi0) ) / den;
    b = (-alpha0*alpha0*alpha0 * (phi_a1 - phip0*alpha1 - phi0) + alpha1*alpha1*alpha1 * (phi_a0 - phip0*alpha0 - phi0) ) / den;
    
    alpha2 = (-b + sqrt(b*b - 3.0*a*phip0)) / (3.0*a);
    
    for (int i=0; i<n; i++) {
        xnext[i] = xk[i] + alpha2 * dk[i];
    }
    phi_a2 = f(xnext);
    // Hasta cumplir Armijo
    while (phi_a2 > phi0 + c1*alpha2*phip0) {
        alpha0 = alpha1;
        alpha1 = alpha2;
        phi_a0 = phi_a1;
        phi_a1 = phi_a2;
        
        den = (alpha0*alpha0 * alpha1*alpha1 * (alpha1 - alpha0));
        
        a = (alpha0*alpha0 * (phi_a1 - phip0*alpha1 - phi0) - alpha1*alpha1 * (phi_a0 - phip0*alpha0 - phi0) )/ den;
        b = (-alpha0*alpha0*alpha0 * (phi_a1 - phip0*alpha1 - phi0) + alpha1*alpha1*alpha1 * (phi_a0 - phip0*alpha0 - phi0) ) / den;
        
        alpha2 = (-b + sqrt(b*b - 3.0*a*phip0)) / (3.0*a);
        
        for (int i=0; i<n; i++) {
            xnext[i] = xk[i] + alpha2 * dk[i];
        }
        phi_a2 = f(xnext);
    }
    
    free(xnext);
    return alpha2;
}


double dogleg(int n, double (*f)(double*), void (*g)(double*, double*), void (*H)(double*, double **), double * x0, double * xoptimo, double Dmax, int M, double tolg, double told)
{
    double * gprev = createVector(n);
    double * gnext = createVector(n);
    
    double * xprev = createVector(n);
    double * xnext = createVector(n);
    double * Pu = createVector(n);
    double * Pb = createVector(n);
    double * Pc = createVector(n);
    double * P = createVector(n);
    
    
    double ** B = createMatrix(n, n);
    gsl_matrix Bgsl; // Para usarlo con el método de Cholesky de la GSL
    initGslMatrix(&Bgsl, B, n, n);
    gsl_set_error_handler_off ();
    double ** Bo = createMatrix(n, n);
    
    double nu = 0.5;
    double delta = 0.5*Dmax;  // Radio
    copyVector(x0, xprev, n);
    double * vaux = createVector(n);
    
    for (int k=0; k<M; k++) {
        std::cout<<k<<": ";
        // Calula modelo cuadrático
        double f_k = f(xprev);
        
        g(xprev, gprev);
        H(xprev, B);
        copyMatrix(B, Bo, n, n);
        double gnorm = vectorNorm(gprev, n);
        prodMatrizVector(B, gprev, vaux, n);
        double gtBg = productoPunto(gprev, vaux, n);
        
        // Calcula minimo en la dirección del gradiente del modelo cuadrático
        double alpha = - (productoPunto(gprev, gprev, n)) / (gtBg);
        copyVector(gprev, Pu, n);
        scaleVector(alpha, Pu, n);
        
        // Calcula el paso completo, Pb, si B es positiva definida
        int e = gsl_linalg_cholesky_decomp1(&Bgsl);
        
        if(e != GSL_EDOM)   // Si es simetrica definida positiva, usa dogleg
        {
            std::cout<<"Dogleg -> ";
            gsl_linalg_cholesky_invert(&Bgsl);  // Invierte matriz
            
            prodMatrizVector(B, gprev, Pb, n);  // Calcula Pb
            scaleVector(-1, Pb, n);
            
            
            if (vectorNorm(Pb, n) <= delta) {
                copyVector(Pb, P, n);
                std::cout<<" a "<<" | ";
            }
            else    // Calcula tau Dogleg
            {
                // Resuelve  ecuación para obtener tau
                double a, b, c;
                for (int i=0; i<n; i++) {
                    vaux[i] = Pb[i] - Pu[i];
                }
                a = pow(vectorNorm(vaux, n), 2.0);
                b = 2.0 * productoPunto(Pb, vaux, n);
                c = pow(vectorNorm(Pu, n), 2.0) - delta*delta;
                double t1 = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
                double t2 = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a);
                
                alpha = gsl_max(t1, t2) + 1;    // Tau
                
                // Actualiza P
                if (alpha <= 1 && alpha >= 0) {
                    copyVector(Pu, P, n);
                    scaleVector(alpha, P, n);
                }
                else if(alpha>1.0){
                    scaleVector(alpha - 1.0, vaux, n);
                    for (int i=0; i<n; i++) {   // Suma
                        P[i] = Pu[i] + vaux[i];
                    }
                }
                else {
                    if (gtBg <= 0.0)
                        alpha = delta / gnorm;
                    else
                        alpha = gsl_min((gnorm*gnorm) / (gtBg), delta / gnorm);
                    
                    copyVector(gprev, Pc, n);
                    scaleVector(-alpha, Pc, n);
                    copyVector(Pc, P, n); // Actualiza P
                }
                    
                std::cout<<" b "<<" | ";
            }
            
        }
        else // Usamos Cauchy
        {
            std::cout<<"Cauchy -> ";
            if (gtBg <= 0.0)
                alpha = delta / gnorm;
            else
                alpha = gsl_min((gnorm*gnorm) / (gtBg), delta / gnorm);
            
            copyVector(gprev, Pc, n);
            scaleVector(-alpha, Pc, n);
            copyVector(Pc, P, n); // Actualiza P
        }
        
        
        // Calcula ro_k, la medida de ajuste del modelo
        for (int i=0; i<n; i++) {
            xnext[i] = xprev[i] + P[i];
        }
        double num = (f_k - f(xnext));
        prodMatrizVector(Bo, P, vaux, n);
        double lo = productoPunto(P, vaux, n);
        double den = -(productoPunto(gprev, P, n) + 0.5*lo);
        
        double ro = num / den;
        
        std::cout<<ro<<" | ";
        // Ajustamos xnext
        if (ro < nu) {
            copyVector(xprev, xnext, n);
            delta = 0.25 * delta;
        }
        else if (ro > 0.75 && vectorNorm(P, n) == delta) {
            delta = gsl_min(2.0*delta, Dmax);
        }
        
        // Verificar tolerancias e imprimir datos
        g(xnext, gnext);
        
        if (vectorNorm(gnext, n) < tolg) {
            std::cout<<"Se cumplió tolerancia tolg "<<vectorNorm(gnext, n)<<"\n";
            break;
        }
        if(delta < told) {
            std::cout<<"Se cumplió tolerancia región de confianza "<<"\n";
            break;
        }
        std::cout<<vectorNorm(gnext, n);
        std::cout<<"\n";
        copyVector(xnext, xprev, n);
        
    }
    copyVector(xnext, xoptimo, n);
    
    free(gprev);
    free(gnext);
    free(vaux);
    free(Pu);
    free(Pc);
    free(Pb);
    freeMatrix(B);
    freeMatrix(Bo);
    free(xprev);
    free(xnext);
    free(P);
    return f(xoptimo);
}


double conjGrad(int n, double (*f)(double*), void (*g)(double*, double*), double ** Q, double * x0, double * xoptimo, double tolg, double tolx, double tolf)
{
    double * gnext = createVector(n);
    double * gprev = createVector(n);
    double * aux = createVector(n);
    double * d = createVector(n);
    double * xprev = createVector(n);
    double * xnext = createVector(n);
    
    
    double alpha = 0.0;
    double beta = 0.0;
    
    // Inicializa gradiente de x0
    g(x0, gprev);
    copyVector(x0, xprev, n);   // Inicializa xk
    
    // Inicializa dirección
    copyVector(gprev, d, n);
    scaleVector(-1.0, d, n);
    
    // Ciclo de a lo más n
    for (int k=0; k<n; k++) {
        // Actualiza alpha
        prodMatrizVector(Q, d, aux, n);
        double dQd = productoPunto(d, aux, n);
        double gg = productoPunto(gprev, gprev, n);
        alpha = gg / dQd;
        
        // Actualiza x
        copyVector(xprev, xnext, n);
        for (int i=0; i<n; i++) {
            xnext[i] += alpha*d[i];
        }
        
        // Actualiza gradiente
        g(xnext, gnext);
        
        // Actualizando beta
        beta = productoPunto(gnext, gnext, n) / gg;
        
        // Actualizando d
        scaleVector(beta, d, n);
        vectorDif(d, gnext, d, n);
        
        // Verificando convergencia
        double difg = vectorDifNorm(gprev, gnext, n);
        double diff = abs(f(xprev)-f(xnext));
        double difx = vectorDifNorm(xnext, xprev, n);
        
        std::cout<<k<<": "<<difx<<"   "<<diff<<"   "<<difg<<"\n";
        
        if(difx<tolx)
            break;
        if(diff < tolf)
            break;
        if(difg < tolg)
            break;
        
        copyVector(xnext, xprev, n);
        copyVector(gnext, gprev, n);
        
    }
    copyVector(xnext, xoptimo, n);
    free(aux);
    free(d);
    free(gnext);
    free(gprev);
    free(xnext);
    free(xprev);
    return f(xoptimo);
}

double conjGradNL(int n, double (*f)(double*), void (*g)(double*, double*), const char metodo[], double * x0, double * xoptimo, double tolg, double tolx, double tolf, double M)
{
    double * gnext = createVector(n);
    double * gprev = createVector(n);
    double * aux = createVector(n);
    double * d = createVector(n);
    double * xprev = createVector(n);
    double * xnext = createVector(n);
    
    
    double alpha = 0.0;
    double beta = 0.0;
    
    int ind = 0; //Método
    std::string msg(metodo);
    if(metodo == std::string("PR"))
        ind = 1;
    else if (metodo == std::string("HS"))
        ind = 2;
    
    // Inicializa gradiente de x0
    g(x0, gprev);
    copyVector(x0, xprev, n);   // Inicializa xk
    
    // Inicializa dirección
    copyVector(gprev, d, n);
    scaleVector(-1.0, d, n);
    alpha = 1.0;
    // Ciclo de a lo más M
    for (int k=0; k<M; k++) {
        // Actualiza alpha usando backtracking
        alpha = backtrackingStep(alpha, xprev, d, gprev, f, n);
        
        // Actualiza x
        copyVector(xprev, xnext, n);
        for (int i=0; i<n; i++) {
            xnext[i] += alpha*d[i];
        }
        
        // Actualiza gradiente
        g(xnext, gnext);
        
        // Actualizando beta dependiendo del método
        if (ind == 0) { // Fletcher-Reeves
            beta = productoPunto(gnext, gnext, n) / productoPunto(gprev, gprev, n);
        }
        else if (ind == 1) {    // Polak-Ribiere
            vectorDif(gnext, gprev, aux, n);
            beta = productoPunto(gnext, aux, n) / productoPunto(gprev, gprev, n);
        }
        else {  // Hestenes-Stiefel
            vectorDif(gnext, gprev, aux, n);
            beta = productoPunto(gnext, aux, n) / productoPunto(aux, d, n);
        }
        
        // Actualizando d
        scaleVector(beta, d, n);
        vectorDif(d, gnext, d, n);
        
        // Verificando convergencia
        double normg = vectorDifNorm(gprev, gnext, n);
        double diff = abs(f(xprev)-f(xnext));
        double difx = vectorDifNorm(xnext, xprev, n);
        double normd = vectorNorm(d, n);
        
        std::cout<<k<<"  "<<normg<<"   "<<f(xnext)<<"  "<<normd<<"\n";
        
        if(difx<tolx) {
            std::cout<<"Se cumplió tolerancia en x\n";
            break;
        }
        if(diff < tolf) {
            std::cout<<"Se cumplió tolerancia en f\n";
            break;
        }
        if(normg < tolg) {
            std::cout<<"Se cumplió tolerancia en g\n";
            break;
        }
        copyVector(xnext, xprev, n);
        copyVector(gnext, gprev, n);
        
    }
    copyVector(xnext, xoptimo, n);
    free(aux);
    free(d);
    free(gnext);
    free(gprev);
    free(xnext);
    free(xprev);
    return f(xoptimo);
}


void solveSENLNewton(double *x0, double *x, int n, void (*F)(double*, double*), void (*J)(double*, double**), double tol, int M)
{
    double * s;
    double * xnext = createVector(n);
    double * xprev = createVector(n);
    double * Fx = createVector(n);
    double ** Jac = createMatrix(n, n);
    gsl_matrix Jac_gsl;
    initGslMatrix(&Jac_gsl, Jac, n, n);
    
    double ** JacI = createMatrix(n, n);
    gsl_matrix JacI_gsl;
    initGslMatrix(&JacI_gsl, JacI, n, n);
    
    gsl_permutation * permutation = gsl_permutation_alloc(3);
    
    copyVector(x0, xprev, n);
    s = createVector(n);
    // Calcula F(xprev)
    F(xprev, Fx);
    
    for (int k=0; k< M; k++) {
        int sgn;
        double invnorm, norm;
        scaleVector(-1.0, Fx, n);
        // Calcula Jacobian
        J(xprev, Jac);
        norm = matrixInfNorm(Jac, n, n);    // Norma usada para _k_
        
        // Resuelve para s
        gsl_linalg_LU_decomp(&Jac_gsl, permutation, &sgn);
        gsl_linalg_LU_invert(&Jac_gsl, permutation, &JacI_gsl);
//        gsl_linalg_LU_solve(&Jac_gsl, &permutation, &Fx_gsl, &s_gsl);
        prodMatrizVector(JacI, Fx, s, n);
        
        invnorm = matrixInfNorm(JacI, n, n); // Norma usada para _k_
        
        // Calcula xnext
        vectorSum(xprev, s, xnext, n);
        
        // Calcula F(xnext)
        F(xnext, Fx);
        
        // Imprime datos
        std::cout<<k<<": x = ";
        for (int i=0; i<n; i++) {
            std::cout<<xnext[i]<<", ";
        }
        double vn = vectorNorm(Fx, n);
        std::cout<<"Norma: "<<vn<<", Condición: "<<norm*invnorm<<".\n";
        copyVector(xnext, xprev, n);
        
        // Viendo si cumple la tolerancia
        if(vn < tol)
            break;
        
    }
    
    copyVector(xnext, x, n);
    
    freeMatrix(Jac);
    freeMatrix(JacI);
    free(s);
    free(xnext);
    free(xprev);
    free(Fx);
    
}


void solveSENLBroyden2(double *x0, double *x, int n, void (*F)(double*, double*), double ** J0, double tol, int M)
{
    double * s = createVector(n);
    double * y = createVector(n);
    double * xnext = createVector(n);
    double * xprev = createVector(n);
    double * Fxprev = createVector(n);
    double * Fxnext = createVector(n);
    double * aux = createVector(n);
    
    // Inicializa matriz J
    gsl_matrix *Jac_gsl = gsl_matrix_alloc(n, n);
    double ** J = (double**)malloc(n * sizeof(double*));
    J[0] = Jac_gsl->data;
    for(int i=1; i<n; i++)
        J[i] = J[i-1] + n;
    
    gsl_permutation * permutation = gsl_permutation_alloc(3);
    
    // Inicializa matriz Jnext
    double ** Jnext = createMatrix(n, n);
    gsl_matrix *JacI_gsl = gsl_matrix_alloc(n, n);
    
    // Invirtiendo la matriz inicial
    copyMatrix(J0, J, n, n);
    printMatrix(J, n, n);
    int sgn;
    gsl_linalg_LU_decomp(Jac_gsl, permutation, &sgn);
    gsl_linalg_LU_invert(Jac_gsl, permutation, JacI_gsl);
    gsl_matrix_memcpy(Jac_gsl, JacI_gsl);
    
    
    copyVector(x0, xprev, n);
    F(xprev, Fxprev);
    
    for (int k=0; k<M; k++) {
        printMatrixGSL(Jac_gsl);
        std::cout<<"\n";
        // Actualiza s_k
        prodMatrizVector(J, Fxprev, s, n);
        scaleVector(-1.0, s, n);
        
        // Actualiza x_{k+1}
        vectorSum(xprev, s, xnext, n);
        
        // Crea y_k
        F(xnext, Fxnext);
        scaleVector(-1.0, Fxprev, n);
        vectorSum(Fxnext, Fxprev, y, n);
        scaleVector(-1.0, Fxprev, n);
        
        // Actualiza J
        double yy = productoPunto(y, y, n);
        prodMatrizVector(J, y, aux, n);
        scaleVector(-1.0, aux, n);
        
        vectorSum(s, aux, aux, n);
        scaleVector(1.0/yy, aux, n);
        
        prodVectorVector(aux, y, Jnext, n);
        squareMatrixSum(J, Jnext, J, n);
        
        copyVector(xnext, xprev, n);
        
        
        // Imprime datos
        std::cout<<k<<": x = ";
        for (int i=0; i<n; i++) {
            std::cout<<xnext[i]<<", ";
        }
        double vn = vectorNorm(Fxnext, n);
//        double norm = matrixInfNorm(J, n, n);
//        gsl_linalg_LU_decomp(Jac_gsl, permutation, &sgn);
//        gsl_linalg_LU_invert(Jac_gsl, permutation, &JacI_gsl);
//        double normI = matrixInfNorm(Jnext, n, n);
        std::cout<<"Norma: "<<vn<<", Condición: "<<".\n";
        copyVector(xnext, xprev, n);
        
        // Viendo si cumple la tolerancia
        if(vn < tol)
            break;
    }
    
    free(s);
    free(y);
    free(xnext);
    free(xprev);
    free(Fxprev);
    free(Fxnext);
    free(aux);
    free(J);
    gsl_matrix_free(Jac_gsl);
    freeMatrix(Jnext);
}


void solveSENLBroyden(double *x0, double *x, int n, void (*F)(double*, double*), double ** J0, double tol, int M)
{
    double * s;
    double * xnext = createVector(n);
    double * xprev = createVector(n);
    double * Fprev = createVector(n);
    double * Fnext = createVector(n);
    double * y = createVector(n);
    gsl_matrix *Jac_gsl = gsl_matrix_alloc(n, n);
    double ** Jac = (double**)malloc(n * sizeof(double*));
    Jac[0] = Jac_gsl->data;
    for(int i=1; i<n; i++)
        Jac[i] = Jac[i-1] + n;
    
    double * aux = createVector(n);
    
    double ** JacI = createMatrix(n, n);
    gsl_matrix JacI_gsl;
    initGslMatrix(&JacI_gsl, JacI, n, n);
    
    gsl_permutation * permutation = gsl_permutation_alloc(3);
    
    copyVector(x0, xprev, n);
    s = createVector(n);
    // Calcula F(xprev)
    F(xprev, Fprev);
    
    copyMatrix(J0, Jac, n, n);
    
    for (int k=0; k< M; k++) {
        int sgn;
        double invnorm, norm;
        // Calcula
        copyMatrix(Jac, J0, n, n);
        norm = matrixInfNorm(Jac, n, n);    // Norma usada para _k_
        
        // Resuelve para s
        gsl_linalg_LU_decomp(Jac_gsl, permutation, &sgn);
        gsl_linalg_LU_invert(Jac_gsl, permutation, &JacI_gsl);
        //        gsl_linalg_LU_solve(&Jac_gsl, &permutation, &Fx_gsl, &s_gsl);
        scaleVector(-1.0, Fprev, n);
        prodMatrizVector(JacI, Fprev, s, n);
        
        invnorm = matrixInfNorm(JacI, n, n); // Norma usada para _k_
        
        // Calcula xnext
        vectorSum(xprev, s, xnext, n);
        
        // Calcula F(xnext) y y_k
        F(xnext, Fnext);
        vectorSum(Fnext, Fprev, y, n);
        
        // Actualiza J
        prodMatrizVector(Jac, s, aux, n);
        scaleVector(-1.0, aux, n);
        vectorSum(y, aux, aux, n);
        scaleVector(1.0/productoPunto(y, y, n), aux, n);
        prodVectorVector(aux, s, Jac, n);
        squareMatrixSum(Jac, J0, Jac, n);
        
        
        // Imprime datos
        std::cout<<k<<": x = ";
        for (int i=0; i<n; i++) {
            std::cout<<xnext[i]<<", ";
        }
        double vn = vectorNorm(Fprev, n);
        std::cout<<"Norma: "<<vn<<", Condición: "<<norm*invnorm<<".\n";
        copyVector(xnext, xprev, n);
        F(xprev, Fprev);
        // Viendo si cumple la tolerancia
        if(vn < tol)
            break;
        
    }
    
    copyVector(xnext, x, n);
    free(s);
    free(y);
    free(xnext);
    free(xprev);
    free(Fprev);
    free(Fnext);
    free(aux);
    gsl_matrix_free(Jac_gsl);
    freeMatrix(JacI);
}


void BFGS(double *x0, double *x, int n, double (*f)(double*), void (*g)(double*, double*), double ** H0, double tol, int M)
{
    double * xnext = createVector(n);
    double * xprev = createVector(n);
    double * gprev = createVector(n);
    double * gnext = createVector(n);
    double * y = createVector(n);
    double * s = createVector(n);
    double * p = createVector(n);
    double ** Hprev = createMatrix(n, n);
    double ** Hnext = createMatrix(n, n);
    double ** M1 = createMatrix(n, n);
    double ** M2 = createMatrix(n, n);
    
    double ** aux;
    
    // Igualando valores iniciales
    copyVector(x0, xnext, n);
    g(xnext, gnext);
    copyMatrix(H0, Hnext, n, n);
    std::cout<<"BFGS: "<<":  ("<<xnext[0]<<", "<<xnext[1]<<"),  "<<vectorNorm(gnext, n)<<",  "<<f(xnext)<<"\n";
    
    for (int k=0; k<M; k++) {
        // Actualizando
        copyVector(xnext, xprev, n);
        copyVector(gnext, gprev, n);
        aux = Hprev;
        Hprev = Hnext;
        Hnext = aux;
        
        // Verificando condición
        double norma = vectorNorm(gprev, n);
        if(norma < tol)
            break;
        
        // Calculando dirección
        prodMatrizVector(Hprev, gprev, p, n);
        scaleVector(-1.0, p, n);
        
        // Asegurando que la dirección es de descenso
        double pg = productoPunto(p, gprev, n);
        if (pg>0) {
            double factor = pg / productoPunto(gprev, gprev, n);
            for (int i=0; i<n; i++)     // Perturbando matriz
                Hprev[i][i] += 0.00001 + factor;
            prodMatrizVector(Hprev, gprev, p, n);   // Calculando nueva dirección
            scaleVector(-1.0, p, n);
        }
        
        // Calculando tamaño de paso y actualizando xnext y gnext
        double alpha = backtrackingStep(1.0, xprev, p, gprev, f, n);
        scaleVector(alpha, p, n);
        vectorSum(xprev, p, xnext, n);
        g(xnext, gnext);
        // std::cout<<alpha;
        
        // Actualizando s y y
        for (int i=0; i<n; i++) {
            s[i] = xnext[i] - xprev[i];
            y[i] = gnext[i] - gprev[i];
        }
        
        // Actualizando H
        double ys = productoPunto(y, s, n);
        // Falta revisar criterio de curvatura
        if (ys<=0) {
            double lambda = ys/productoPunto(s, s, n);
            copyMatrix(Hprev, Hnext, n, n);
            for (int i=0; i<n; i++)
                Hnext[i][i] += 0.00001 - lambda;
                
        }
        else {
            double ro = 1.0 / ys;
            clearMatrix(M1, n, n);
            prodVectorVector(s, y, M1, n);
            scaleMatrix(-ro, M1, n, n);
            
            clearMatrix(M2, n, n);
            prodVectorVector(y, s, M2, n);
            scaleMatrix(-ro, M2, n, n);
            
            for (int i=0; i<n; i++) {
                M1[i][i] += 1.0;
                M2[i][i] += 1.0;
            }
            productoMatriz(M1, Hprev, Hnext, n);
            productoMatriz(Hnext, M2, M1, n);   // Triple producto en M1
            prodVectorVector(s, s, M2, n);
            scaleMatrix(ro, M2, n, n);          // Segundo término en M2
            
            for (int i=0; i<n; i++) {   // Suma M2 + M1 y guardar en Hnext
                for (int j=0; j<n; j++)
                    Hnext[i][j] = M1[i][j] + M2[i][j];
            }
        }
        
        
        // Imprimiendo k+1, f(xnext), ||gnext|| y xnext
        norma = vectorNorm(gnext, n);
        std::cout<<" "<<k+1<<":  ("<<xnext[0]<<", "<<xnext[1]<<"),  "<<norma<<",  "<<f(xnext)<<"\n";
    }
    
    // Guardando punto óptimo
    copyVector(xnext, x, n);
    
    // Liberando memoria
    free(xnext);
    free(xprev);
    free(gprev);
    free(gnext);
    free(y);
    free(s);
    free(p);
    freeMatrix(Hprev);
    freeMatrix(Hnext);
    freeMatrix(M1);
    freeMatrix(M2);
}


void hessianFiniteDif(double **H, double (*f)(double*), double *x, int n, double h)
{
    double * xaux = createVector(n);
    copyVector(x, xaux, n);
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            H[i][j] = f(x);
            
            xaux[i] += h;
            H[i][j] -= f(xaux);
            xaux[i] -= h;
            
            xaux[j] += h;
            H[i][j] -= f(xaux);
            xaux[j] -= h;
            
            xaux[i] += h;
            xaux[j] += h;
            H[i][j] += f(xaux);
            
            copyVector(x, xaux, n);
        }
    }
    
}


void LM_method(double *x0, double *x, int n, int m, void (*R)(double*, double*), void (*J)(double*, double**), double nu, double tol, int M)
{
    double * xprev = createVector(n);
    double * xnext = createVector(n);
    double * p = createVector(n);
    double * Rx = createVector(m);
    double * JR = createVector(n);
    
    double fnext;
    double fprev;
    double lambda = 1.0;
    
    double ** Jac = createMatrix(m, n);
    double ** JacT = createMatrix(n, m);
    
    gsl_matrix *JJ_gsl = gsl_matrix_alloc(n, n);
    double ** JJ = (double**)malloc(n * sizeof(double*));
    JJ[0] = JJ_gsl->data;
    for(int i=1; i<n; i++)
        JJ[i] = JJ[i-1] + n;
    gsl_permutation * permutation = gsl_permutation_alloc(n);
    
    double ** JJI = createMatrix(n, n);
    gsl_matrix JJI_gsl;
    initGslMatrix(&JJI_gsl, JJI, n, n);
    
    // Inicializando x y R(x)
    copyVector(x0, xprev, n);
    R(xprev, Rx);
    J(xprev, Jac);
    transpose2(Jac, JacT, m, n);
    producto_matricial2(JacT, Jac, JJ, n, m, n);
    for (int i=0; i<n; i++) {
        if(JJ[i][i] > lambda)
            lambda = JJ[i][i];
    }
    
    for (int k=0; k<M; k++) {
        // Actualizando valores
        fprev = 0.5 * productoPunto(Rx, Rx, m);
        J(xprev, Jac);
        transpose2(Jac, JacT, m, n);
        producto_matricial2(JacT, Jac, JJ, n, m, n);
        prodMatrizVector2(JacT, Rx, JR, n, m);  // Obteniendo gradiente de f en xprev
        
        // Imprimiendo datos en la pantalla k, x_1, x_2, lambda, f(prev), ||grad f(xprev)||
        std::cout<<k<<":  "<<xprev[0]<<", "<<xprev[1]<<";  "<<lambda<<", "<<fprev<<", "<<vectorNorm(JR, n)<<"\n";
        
        // Condición de paro
        if (vectorNorm(JR, n)<tol)
            break;
        
        // Obteniendo dirección
        for (int i=0; i<n; i++) {
            JJ[i][i] += lambda;
        }
        int signo;
        gsl_linalg_LU_decomp(JJ_gsl, permutation, &signo);
        gsl_linalg_LU_invert(JJ_gsl, permutation, &JJI_gsl);
        prodMatrizVector(JJI, JR, p, n);
        scaleVector(-1.0, p, n);
        
        // Actualizando x, R(x) y f(x)
        vectorSum(xprev, p, xnext, n);
        R(xnext, Rx);
        fnext = 0.5 * productoPunto(Rx, Rx, m);
        
        // Actualizando lambda
        if (fnext >= fprev) {
            copyVector(xprev, xnext, n);
            R(xnext, Rx);
            lambda = lambda * nu;
        }
        else
            lambda = lambda / nu;
        
        // Actualizando xprev
        copyVector(xnext, xprev, n);
    }
    
    copyVector(xprev, x, n);
    free(xprev);
    free(xnext);
    free(p);
    free(Rx);
    free(JR);
    freeMatrix(Jac);
    freeMatrix(JacT);
    gsl_matrix_free(JJ_gsl);
    free(JJ);
    freeMatrix(JJI);
}

void quadraticPenaltyMethod(double (*Q)(double*, double), void (*gradQ)(double*, double*, double), void (*g)(double*, double*), void (*h)(double*, double*), int n, int m, int r, double *x0, double tol, double N)
{
    double *x = createVector(n);
    double *gx = createVector(m);
    double *hx = createVector(r);
    
    copyVector(x0, x, n);
    double mu = 2.0;
    double tau = 0.001;
    double px;
    
    for (int k=0; k<N; k++) {
        // Optimiza Q
        std::cout<<"Iter: "<<k<<", mi: "<<mu<<", tao:  "<<tau<<" ";
        
        line_search_for_QPM(n, Q, gradQ, mu, x, x, 1000, 0.8, tau, 0.0, 0.0);
        // Calcula p(x)
        g(x, gx);
        h(x, hx);
        for (int i=0; i<r; i++) {
            hx[i] = (hx[i] > 0.0)? (-hx[i]) : 0.0;
        }
        
        px = 0.5*mu*( productoPunto(gx, gx, m) + productoPunto(hx, hx, r));
        
        std::cout<<px<<"\n";
        // Condición de paro
        if (px < tol) {
            break;
        }
        
        mu = 2.0*mu;
        tau = 0.8*tau;
        
    }
    
    copyVector(x, x0, n);
    
    free(x);
    free(gx);
    free(hx);
}

double line_search_for_QPM(int n, double (*f)(double*, double mu), void (*g)(double*, double*, double), double mu, double * x0, double * xoptimo, int M, double alpha0, double tolg, double tolx, double tolf) {
    double * xnext = createVector(n);
    double * xprev = createVector(n);
    double * vgrad = createVector(n);
    double * aux = createVector(n);
    double * dir = createVector(n);
    double alpha;
    
    std::cout<<"Método backtracking ";
    
    copyVector(x0, xprev, n);
    int k;
    // Itera mientras no se cumplan los criterios de paro
    for (k=0; k<M; k++) {
        alpha = alpha0;            // Alpha inicial
        g(xprev, vgrad, mu);        // Obtenemos gradiente para dirección
        // Obtenemos la dirección de máximo descenso
        copyVector(vgrad, dir, n);
        scaleVector(-1.0, dir, n);
        
        // Seleccionamos método para obtener alpha
        alpha = backtrackingStepForQPM(alpha, xprev, dir, vgrad, f, n, mu);// Se busca alpha
        
        // Actualizamos x_k con el alpha encontrado
        for (int j=0; j<n; j++)                             // Se obtiene el x_{k+1}
            xnext[j] = xprev[j] + alpha * dir[j];
        
        /* Condiciones de paro */
        // Primera condición
        for (int j=0; j<n; j++)
            aux[j] = xprev[j] - xnext[j];
        double norm = vectorNorm(aux, n);
        //        std::cout<<norm<<"  ";
        if(norm/(std::max(1.0, vectorNorm(xprev, n))) < tolx) {
            std::cout<<"Se cumplió tolerancia en x_k - x_(k+1)\n";
            break;
        }
        
        // Segunda condición
        double fprev = f(xprev, mu);
        double fnext = f(xnext, mu);
        
        norm = abs(fprev-fnext);
        //        std::cout<<norm<<"  ";
        if(norm/(std::max(1.0, abs(fnext))) < tolf) {
            std::cout<<"Se cumplió tolerancia en f_k - f_(k+1)\n";
            break;
        }
        
        // Tercera condición
        g(xnext, vgrad, mu);
        norm = vectorNorm(vgrad, n);
        //std::cout<<norm<<"  ";
        if(norm < tolg) {
            std::cout<<"Se cumplió tolerancia del valor del gradiente en x\n\n";
            break;
        }
        
        copyVector(xnext, xprev, n);        // Siguiente paso
    }
    copyVector(xnext, xoptimo, n);      // Copia el x*
    
    free(xprev);
    free(xnext);
    free(vgrad);
    free(aux);
    free(dir);
    return f(xoptimo, mu);                  // Regresa evaluación de x* en f
}

double backtrackingStepForQPM(double alpha0, double * xk, double * dk, double * gk, double (*f)(double*, double), int n, double mu)
{
    double c = 0.0001;
    double p = 0.27;
    
    double alpha = alpha0;
    
    double fxa;
    double fx;
    double m;
    double * xnext = createVector(n);
    
    m = c * alpha * productoPunto(gk, dk, n);
    fx = f(xk, mu);
    for (int i=0; i<n; i++) {
        xnext[i] = xk[i] + alpha * dk[i];
    }
    fxa = f(xnext, mu);
    
    while (fxa > fx + m) {
        
        alpha = p*alpha;
        m = c * alpha * productoPunto(gk, dk, n);
        for (int i=0; i<n; i++) {
            xnext[i] = xk[i] + alpha * dk[i];
        }
        fxa = f(xnext, mu);
    }
    free(xnext);
    return alpha;
}



/// ---------   Proyecto  ----------


double cso(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration,
           int G, double FL)
{
    gsl_rng *rng_ptr; // pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    std::srand ( unsigned(std::time(0)) );
    
    const double eps = 0.001;
    
    int rn = ceil(0.15 * n);    // Number of roosters
    int hn = ceil(0.7 * n);     // Number of heins
    int cn = n - rn - hn;       // Number of chicks
    int mn = ceil(0.2 * n);     // Numer of mother heins
    
    // Control of hierarchy
    int roosters[rn];
    int hens[hn];
    int chicks[cn];
    int mothers[mn];
    int chicksMum[cn];  // Save the mums of the chicks
    int hensLider[hn];// Save the hens' rooster lider
    
    // Setting chickens
    double ** agents = createMatrix(n, dimension);
    double ** pbest = createMatrix(n, dimension);
    for (int i=0; i<n; i++) {
        for (int j=0; j<dimension; j++)
            agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
    }
    copyMatrix(agents, pbest, n, dimension);
    
    // Setting fitness
    double * fitness = createVector(n);
    for (int i=0; i<n; i++) {
        fitness[i] = function(agents[i]);
    }
    double * pfit = createVector(n);
    copyVector(fitness, pfit, n);
    
    // To save the best
    double * Pbest = createVector(dimension);
    double* i1; i1 = std::min_element(fitness, fitness+n);
    copyVector(agents[(i1-fitness)], Pbest, dimension);
    double * Gbest = createVector(dimension); copyVector(Pbest, Gbest, dimension);
    
    // Iterations
    for (int k=0; k<iteration; k++) {
        // If is time to update hierarchies
        if (k % G == 0) {
            update_relationship(n, function, agents, rn, hn, cn, mn, roosters, hens, chicks, mothers, chicksMum, hensLider);
        }
        
        // Search for food
        
        // Roosters
        for (int i=0; i<rn; i++) {
            int ri = roosters[i];
            // Pick a different rooster
            int r2 = roosters[(rand() % rn)];
            while (ri == r2) {
                r2 = roosters[(rand() % rn)];
            }
            
            // Compute variance
            double sigma;
            if (pfit[ri] <= pfit[r2])
                sigma = 1.0;
            else
                sigma = sqrt( exp((pfit[r2] - pfit[ri]) / (abs(pfit[ri]) + eps)) );
            
            // Updating position of the i-th rooster
            for (int d = 0; d<dimension; ++d)
                agents[ri][d] = pbest[ri][d]*(1.0 + gsl_ran_gaussian(rng_ptr, sigma));
            
        }
        
        // Hens
        for (int i=0; i<hn; i++) {
            int heni = hens[i];
            
            int henLider = hensLider[i];
            // Choose chicken to steal food
            int aux = (rand() % hn);
            int a = hens[aux]; int ar = hensLider[aux];
            int b = roosters[(rand() % rn)];
            int chicken = (rand()%2)? a : b ;
            ar = ((chicken == a)? ar:b);
            while (henLider == ar) {    // Not in the same group
                aux = (rand() % hn);
                a = hens[aux]; ar = hensLider[aux];
                b = roosters[(rand() % rn)];
                chicken = (rand()%2)? a : b ;
                ar = ((chicken == a)? ar:b);
            }
            
            // Calculate s1 and s2
            double s1 = exp((pfit[heni] - pfit[chicken]) / (abs(pfit[heni]) + eps));
            double s2 = exp(pfit[chicken] - pfit[heni]);
            if(isinf(s2)) s2 = DBL_MAX;
            
            double rand1;
            double rand2;
            // Update position of the ith hen
            for (int d = 0; d<dimension; ++d) {
                rand1 = gsl_ran_flat(rng_ptr, 0.0, 1.0);
                rand2 = gsl_ran_flat(rng_ptr, 0.0, 1.0);
                agents[heni][d] = pbest[heni][d]
                                + s1 * rand1 * (pbest[henLider][d] - pbest[heni][d])
                                + s2 * rand2 * (pbest[chicken][d] - pbest[heni][d]);
            }
        }
        // Chicks
        for (int i=0; i<cn; i++) {
            for (int d = 0; d<dimension; ++d)
                agents[chicks[i]][d] = pbest[chicks[i]][d] * FL * ( pbest[chicksMum[i]][d] - pbest[chicks[i]][d]);
        }
        
        // Checking constrains
        for (int i=0; i<n; i++) {
            for (int j=0; j<dimension; j++) {
                if (agents[i][j] < lb[j])
                    agents[i][j] = lb[j];
                else if (agents[i][j] > ub[j])
                    agents[i][j] = ub[j];
            }
        }
        
        // Update fitness
        for (int i=0; i<n; i++) {
            fitness[i] = function(agents[i]);
        }
        
        // Checking if there were progress
        for (int i=0; i<n; i++) {
            if (fitness[i] < pfit[i]) {
                pfit[i] = fitness[i];
                copyVector(agents[i], pbest[i], dimension);
            }
        }
        
        // Saving the best of the fitness
        i1 = std::min_element(fitness, fitness+n);
        copyVector(agents[(i1-fitness)], Pbest, dimension);
        if (function(Pbest) < function(Gbest))
            copyVector(Pbest, Gbest, dimension);
    }
    
    // Answer
    copyVector(Gbest, x, dimension);
    double bestValue = function(Gbest);
    
    gsl_rng_free(rng_ptr);
    freeMatrix(agents);
    freeMatrix(pbest);
    free(fitness);
    free(pfit);
    free(Pbest);
    free(Gbest);
    
    return bestValue;
}

void update_relationship(int n, double (*function)(double *), double ** agents, int rn, int hn, int cn, int mn, int roosters[], int hines[], int chicks[], int mothers[], int chicksMum[], int hinesLider[])
{
    // Numerate the chickens and sort them by fitness values
    std::vector< std::pair <double,int> > vect;
    for (int i=0; i<n; i++)
        vect.push_back( std::make_pair(function(agents[i]), i) );
    sort(vect.begin(), vect.end());
    
    for (int i=0; i<rn; i++)        // Saving roosters
        roosters[i] = vect[i].second;
    
    for (int i=rn; i<(rn+hn); i++)  // Saving hines
        hines[i-rn] = vect[i].second;
    
    for (int i=rn+hn; i<n; i++)     // Saving chicks
        chicks[i-rn-hn] = vect[i].second;
    
    // Shuffling hines to assign mothers
    std::random_shuffle(hines, hines+hn);
    for (int i=0; i<mn; i++)
        mothers[i] = hines[i];
    
    // Assign every chick to its mummy
    for (int i=0; i<cn; i++) {
        int rnd = rand() % mn;
        chicksMum[i] = mothers[rnd];
    }
    
    // Assign every chicken to a rooster to be its lider
    for (int i=0; i<hn; i++) {
        int rnd = rand() % rn;
        hinesLider[i] = roosters[rnd];
    }
    
}
