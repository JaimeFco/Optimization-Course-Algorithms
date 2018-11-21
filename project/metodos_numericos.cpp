//
//  metodos_numericos.cpp
//  proyecto00
//
//  Created by Jaime Francisco Aguayo on 10/09/17.
//  Copyright © 2017 Jaime Francisco Aguayo. All rights reserved.
//

#include "metodos_numericos.hpp"

double **createMatrix(int nr, int nc)
{
    int i;
    double** mat;
    
    mat = (double**)malloc(nr * sizeof(double*));
    if(mat==NULL) return(NULL);
    mat[0] = (double*)malloc(nr * nc * sizeof(double));
    if(mat[0]==NULL) return(NULL);
    
    for(i=1; i<nr; i++)
        mat[i] = mat[i-1] + nc;
    
    return mat;
}

void freeMatrix(double ** mat)
{
    free(mat[0]);
    free(mat);
}

double * createVector(int n)
{
    double* aux;
    aux = (double*)malloc(sizeof(double)*n);
    return aux;
}

void freeVect(double * vect)
{
    free(vect);
}

int solveDiagonalMatrix(double *D, double *B, double *X, int n)
{
    for (int i = 0; i<n; ++i)
    {
        if(D[i] == 0)
        {
            std::cout<<"ERROR D"<<D[i]<<" = 0\n";
            return 1;
        }
        X[i] = B[i]/D[i];
    }
    return 0;
}

int solveLowerMatrix(double **Mat, int nr, int nc, double *B, double *X)
{
    double suma;
    for(int i = 0; i<nr; ++i)
    {
        suma = 0;
        for(int j = 0; j<i; ++j)
        {
            suma+=Mat[i][j]*X[j];
        }
        if(Mat[i][i] == 0)
            return 1;
        X[i] = (B[i]-suma)/Mat[i][i];
    }
    return 0;
}

int solveUpperMatrix(double **Mat, int nr, int nc, double *B, double *X)
{
    double suma;
    for(int i = nr-1; i>=0; --i)
    {
        suma = 0;
        for(int j = i+1; j<nr; ++j)
        {
            suma+=Mat[i][j]*X[j];
        }
        if(Mat[i][i] == 0)
            return 1;
        X[i] = (B[i]-suma)/Mat[i][i];
    }
    return 0;
}

int fatorizacionLU(double **A, double **L, double **U, int n)
{
    double suma = 0;
    int i, j, k;
    //Diagonal de U = 1
    for (i = 0; i < n; i++)
    {
        U[i][i] = 1;
    }
    
    
    for (j = 0; j < n; j++)
    {
        for (i = j; i < n; i++)
        {
            suma = 0;
            for (k = 0; k < j; k++)
                suma = suma + L[i][k] * U[k][j];
            
            L[i][j] = A[i][j] - suma;
        }
        
        for (i = j; i < n; i++)
        {
            suma = 0;
            for(k = 0; k < j; k++)
                suma = suma + L[j][k] * U[k][i];
            
            if (L[j][j] == 0)   //Evita división por cero
            {
                printf("ERROR: División entre 0, no se pudo encontrar la fatorizaión LU de la matriz\n");
                return 1;
            }
            U[j][i] = (A[j][i] - suma) / L[j][j];
        }
    }
    return 0;
}


double matrixInfNorm(double **A, int nr, int nc)
{
    double max = A[0][0];
    double suma=0.0;
    for(int i = 0; i<nr; i++)
    {
        for(int j = 0; j<nc; j++)
        {
            suma += std::abs(A[i][j]);
        }
        if(suma > max)
            max = suma;
        suma = 0.0;
    }
    return max;
}

void prodMatrizVector(double **M, double *A, double *B, int n)
{
    for(int i=0; i<n; ++i)
    {
        double r = 0.0;
        for(int j = 0; j<n; ++j)
        {
            r += M[i][j]*A[j];
        }
        B[i] = r;
    }
}

void prodMatrizVector2(double **M, double *A, double *B, int n, int m)
{
    for(int i=0; i<n; ++i)
    {
        double r = 0.0;
        for(int j = 0; j<m; ++j)
        {
            r += M[i][j]*A[j];
        }
        B[i] = r;
    }
}

void productoMatriz(double **P, double **Q, double **R, int n)
{
    int i;

    iterate(i, n)
    {
        for(int j = 0; j<i; j++)
        {
            double aux = Q[i][j];
            Q[i][j] = Q[j][i];
            Q[j][i] = aux;
        }
    }
    
    iterate(i, n)
    {
        for(int k = 0; k<n; k++)
        {
            double r = 0.0;
            for(int j = 0; j<n; ++j)
            {
                r += P[i][j]*Q[k][j];
            }
            R[i][k] = r;
        }
    }
}

void producto_matricial(double** P, double** Q, double** R, int n)
{
    //#pragma omp parallel for
    for (int i=0; i<n; i++) {
        //#pragma omp parallel for
        for(int j=0; j<n; j++)  {
            double sum=0;
            //#pragma omp parallel for reduction(+:sum)
            for(int k=0; k<n; k++)
                sum += P[i][k] * Q[k][j];
            R[i][j] = sum;
        }
    }
}

void producto_matricial2(double** P, double** Q, double** R, int n, int m, int p)
{
    //#pragma omp parallel for
    for (int i=0; i<n; i++) {
        //#pragma omp parallel for
        for(int j=0; j<p; j++)  {
            double sum=0;
            //#pragma omp parallel for reduction(+:sum)
            for(int k=0; k<m; k++)
                sum += P[i][k] * Q[k][j];
            R[i][j] = sum;
        }
    }
}

void prodVectorVector(double *A, double *B, double ** C, int n)
{
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++)
            C[i][j] = A[i] * B[j];
    }
}


double errorLU(double **A, double **L, double **U, int size)
{
    double **LU = createMatrix(size, size);
    productoMatriz(L, U, LU, size);
    for(int i = 0; i<size; i++)
        for(int j = 0; j<size; ++j)
            LU[i][j] = A[i][j]-LU[i][j];
    double r = matrixInfNorm(LU, size, size);
    
    freeMatrix(LU);
    return r;
}


int solveTridiagonalMatrix2(double *a, double *b, double *c, int n, double *v, double *x)
{
    double *cprime = new double[n];
    
    cprime[0] = c[0] / b[0];
    x[0] = v[0] / b[0];
    
    for (int i = 1; i < n; i++)
    {
        float r = b[i] - (a[i-1] * cprime[i - 1]);
        if (r==0)
            return 1;
        if(i < (n-1))
            cprime[i] = c[i] / r;
        x[i] = (v[i] - a[i-1] * x[i - 1]) / r;
    }
    
    
    for (int i = n - 2; i-- > 0; )
        x[i] = v[i] - cprime[i] * x[i + 1];
    
    //Libera espacio usado para calcular la solución
    delete [] cprime;
    
    return 0;
}


int solveTridiagonalMatrix(double *a, double *b, double *c, int n, double *d, double *x)
{
    double *b_, *c_, *d_;
    b_ = new double[n];
    c_ = new double[n];
    d_ = new double[n];
    a-=1;
    b_[0] = b[0];
    c_[0] = c[0];
    d_[0] = d[0];
    
    for(int i = 1; i<n-1; ++i)
    {
        b_[i] = b_[i-1]*b[i] - a[i]*c_[i-1];
        c_[i] = b_[i-1]*c[i];
        d_[i] = b_[i-1]*d[i] - a[i]*d_[i-1];
    }
    d_[n-1] = b_[n-2]*d[n-1] - a[n-1]*d_[n-2];
    b_[n-1] = b_[n-2]*b[n-1] - a[n-1]*c_[n-2];
    
    if(b_[n-1] == 0)
        return 2;
    x[n-1] = d_[n-1]/b_[n-1];
    for(int i=n-2; i>=0; --i)
    {
        if(b_[i] == 0)
            return 1;
        x[i] = (d_[i] - c_[i]*x[i+1])/b_[i];
    }
    return 0;
}

double vectorDifNorm(double *u, double *v, int n)
{
    double suma= 0;
    for(int i = 0; i<n; ++i)
    {
        double aux = (u[i]-v[i]);
        suma += (aux*aux);
    }
    return sqrt(suma);
}

double vectorNorm(double *v, int n)
{
    double suma=0.0;
    for(int i=0; i<n; ++i)
        suma+= (v[i]*v[i]);
    return sqrt(suma);
}

double metodoDeLaPotencia(double **A, int N, double *v, double tol, double MAX)
{
    double *y = new double[N];
    double *Aux = new double[N];
    double norma, eigenValue, convergencia;
    int k = 0;
    if(v == NULL)
    {
        v = new double[N];
        for(int i=0; i<N; ++i)
            v[i] = 1.0;
    }
    
    do {
        prodMatrizVector(A, v, y, N);
        //Norma
        norma = vectorNorm(y, N);
        for(int i = 0; i<N; i++)
        {
            v[i] = y[i]/norma;
        }
        eigenValue = vectorMatrizVector(v, A, v, N);
        
        prodMatrizVector(A, y, Aux, N);
        for(int i=0; i<N; ++i)
            Aux[i] -= eigenValue*y[i];
        convergencia = vectorNorm(Aux, N);
        ++k;
        
    } while (k<MAX && convergencia>tol);
    std::cout<<"Se realizaron "<<k<<" iteraciones, para una matriz de "<<N<<"x"<<N<<"\n";
    std::cout<<"Se encontro valor propio "<<eigenValue<<" con error "<<convergencia<<"\n";
    
    delete [] Aux;
    delete [] y;
    return eigenValue;
}

double vectorMatrizVector(double *v1, double **A, double *v2, int n)
{
    double aux;
    double R = 0.0;
    for(int j=0; j<n; ++j)
    {
        aux = 0.0;
        for(int i=0; i<n; ++i)
            aux += v1[i]*A[i][j];
        R += aux*v2[j];
    }
    return R;
}

//REVISAR!!!!!!!!!!!
int solve(double ** A, int n, double *X, double *B)
{
    //Primero intenta Factorizacion LU
    
    double **U = createMatrix(n, n), **L = createMatrix(n, n);
    if(!fatorizacionLU(A, L, U, n))
    {
        double *Aux = new double[n];
        int j = solveLowerMatrix(L, n, n, B, Aux);
        
        if(j)   //Verifica si hubo errores
        {
            freeMatrix(U);
            freeMatrix(L);
            delete [] Aux;
            return 1;
        }
        j = solveUpperMatrix(U, n, n, Aux, X);
        
        if(j)   //Verifica si hubo errores
        {
            freeMatrix(U);
            freeMatrix(L);
            delete [] Aux;
            return 1;
        }
        
        delete [] Aux;
        freeMatrix(U);
        freeMatrix(L);
        return 0;
    }
    freeMatrix(U);
    freeMatrix(L);
    return 1;
    //Intenta Eliminación Gaussiana
}

double productoPunto(double *A, double *B, int n)
{
    double suma=0;
    for(int i=0; i<n; ++i)
    {
        suma += A[i]*B[i];
    }
    return suma;
}


double metodoPotenciaInversa(double **A, int n, double *X, double delta, double tol, double & eigenvalue, int N)
{
    int k = 0; //Maximo de iteraciones
    double r[n];
    double *Y = new double[n];
    double *W = new double[n];
    double ro;
    
    for(int i=0; i<n; ++i)  //Traslada la matriz
        A[i][i] -= delta;
    
    do
    {
        if(solve(A, n, Y, X))
        {
            
            std::cout<<"ERROR RESOLVIENDO MATRIZ Potencia Inversa";
            for(int i=0; i<n; ++i)  //Traslada la matriz
                A[i][i] -= delta;
            exit(3);
        }
        double normy = vectorNorm(Y, n);
        for(int i=0; i<n; i++)
        {
            Y[i] /= normy;
            W[i] = X[i]/normy;
        }
        ro = productoPunto(Y, W, n);   //Hacer esta funcion
        eigenvalue = delta+ro;
        
        for(int i=0; i<n; i++)
            r[i] = W[i] - ro*Y[i];
        std::swap(X, Y);
        ++k;
    } while((vectorNorm(r, n) > tol) && k<N);
    
    if(k%2 == 1)
    {
        for(int i=0; i<n; i++)
        {
            Y[i] = X[i];
        }
        delete [] X;
    }
    else
        delete [] Y;
    
    delete [] W;
    
    for(int i=0; i<n; ++i)  //Traslada la matriz de vuelta
        A[i][i] += delta;
    return vectorNorm(r, n);
}

double ** readMatrix(const char file[], int &nr, int &nc)
{
    std::ifstream myfile(file, std::ios::binary);
    myfile.read((char*)&nr, sizeof(nr));  //Lee numero de renglones de la matriz
    myfile.read((char*)&nc, sizeof(nc));   //Lee numero de columnas de la matriz
    double **M = createMatrix(nr, nc);  //Pide memoria para la matriz
    
    if(M==NULL) //Si no se pudo reservar memoria
    {
        myfile.close();
        return NULL;
    }
    
    myfile.read((char*)M[0], sizeof(M[0])*nr*nc); //Lee matriz del archivo
    myfile.close();
    return M;
}

double * readArray(const char file[], int &n)
{
    std::ifstream myfile(file, std::ios::binary);
    myfile.read((char*)&n, sizeof(n));
    
    double *V = new double[n];
    
    if(V==NULL) //Si no se pudo reservar memoria
    {
        myfile.close();
        return NULL;
    }
    
    myfile.read((char*)V, sizeof(V)*n);
    myfile.close();
    return V;
}

void printMatrix(double **M, int nr, int nc)
{
    for(int i=0; i<nr; ++i)
    {
        for(int j=0; j<nc; j++)
            printf("%4f ", M[i][j]);
        printf("\n");
    }
}

double metodoJacobi(double **A, int n, double ** &eigenvectors, double *eigenvalues, double tol, int N)
{
    int i, j, cont=0;
    double delta, theta, s, c;
    double ** Mat = createMatrix(n, n);
    double **G = createMatrix(n, n);
    double **Aux = createMatrix(n, n);
    clearMatrix(G, n, n);
    for(int k =0; k<n; ++k)
    {
        G[k][k] = 1.0;
        for(int l=0; l<n; ++l)
            Mat[k][l] = A[k][l];
    }
    
    while(1)
    {
        greatValue(Mat, n, i, j);
        if(abs(Mat[i][j])<tol || cont>=N)
            break;
        delta = (Mat[j][j] - Mat[i][i])/(2.0*Mat[i][j]);
        theta = (delta/abs(delta)) / (abs(delta)+sqrt(1+delta*delta));
        c = 1/(sqrt(1+theta*theta));
        s = c*theta;
        rot(Mat, n, i, j, s, c);
        G[i][i] = G[j][j] = c;
        G[i][j] = s;
        G[j][i] = -s;
        productoMatriz(eigenvectors, G, Aux, n);
        std::swap(eigenvectors, Aux);
        G[i][i] = G[j][j] = 1.0;
        G[i][j] = G[j][i] = 0.0;
        cont++;
    }
    for(int k=0; k<n; k++)
        eigenvalues[k] = Mat[k][k];
    std::cout<<"Método iterativo de Jacobi realizó "<<cont<<" iteraciones."<<"\n";
    freeMatrix(G);
    freeMatrix(Aux);
    freeMatrix(Mat);
    return abs(Mat[i][j]);
}


void greatValue(double **mat, int n, int &i, int &j)
{
    double max = abs(mat[0][1]);
    for(int k=0; k<n; k++)
        for(int l=0; l<n; l++)
        {
            if(k != l && abs(mat[k][l]) > max)
            {
                i = k;
                j = l;
                max = abs(mat[k][l]);
            }
        }
}

void rot(double **mat, int n, int i, int j, double s, double c)
{
    double aii = mat[i][i];
    double ajj = mat[j][j];
    double akj, aki;
    mat[i][i] = c*c*aii - 2*s*c*mat[i][j] + s*s*ajj;
    mat[j][j] = c*c*ajj + 2*s*c*mat[i][j] + s*s*aii;
    mat[i][j] = mat[j][i] = (c*c - s*s)*mat[i][j] + s*c*(aii - ajj);
    for(int k=0; k<n; k++)
    {
        if(k != i && k != j)
        {
            aki = mat[k][i];
            akj = mat[k][j];
            mat[k][i] = mat[i][k] = c*aki - s*akj;
            mat[k][j] = mat[j][k] = s*aki + c*akj;
        }
    }
}

void clearMatrix(double **A, int nr, int nc)
{
    for(int i=0; i<nr; ++i)
        for(int j=0; j<nc; j++)
            A[i][j] = 0.0;
}

void clearVector(double *v, int n)
{
    int i;
    iterate(i, n)
        v[i] = 0.0;
}

inline void copyMatrix(double **A, double **B, int nr, int nc)
{
    int i, j;
    iterate(i, nr)
    {
        iterate(j, nc)
            B[i][j] = A[i][j];
    }
}

void copyVector(double *a, double *b, int n)
{
    int i;
    iterate(i, n)
        b[i] = a[i];
}

void writeMatrix(const char name[], double ** M, int nr, int nc)
{
    std::ofstream file (name, std::ios::binary);
    file.write((char*)&nr, sizeof(nr));
    file.write((char*)&nc, sizeof(nc));
    file.write((char*)M[0], sizeof(M[0][0])*nr*nc);
    
    file.close();
}

double metodoGradienteConjugado(double **A, int n, double *x, double *b, double tol, int N)
{
    if(N == 0)
        N = n;
    int k;
    double *p = new double[n];
    double *q = new double[n];
    double *r = new double[n];
    double e, d, alpha, beta;
    int j;
    iterate(j, n)
    {
        r[j] = -b[j];
        p[j] = -b[j];
    }
    clearVector(x, n);
    e = sqrt( productoPunto(r, r, n) / n );
    
    iterate(k, n+1)
    {
        if(e<tol)
            break;
        
        prodMatrizVector(A, p, q, n);
        d = productoPunto(r, r, n);
        alpha = d / productoPunto(p, q, n);
        
        for(int j=0; j<n; j++)
        {
            x[j] = x[j] + alpha*p[j];
            r[j] = r[j] + alpha*q[j];
        }
        
        beta = productoPunto(r, r, n)/d;
        for(int j=0; j<n; j++)
        {
            p[j] = -r[j] + beta*p[j];
        }
        e = sqrt( productoPunto(r, r, n) / n );
    }
    
    std::cout<<"Método del gradiente conjugado realizó "<<k<<" iteraciones\n";
    delete [] p;
    delete [] q;
    delete [] r;
    return e;
}


void vandermondeMatrix(double *x, double **A, int n)
{
    int i, j;
    iterate(i, n)
        iterate(j, n)
        {
            A[i][j] = pow(x[i], j);
        }
}

void polinomioInterpolador(double *x, double *y, double *c, int n)
{
    double ** A = createMatrix(n, n);
    
    vandermondeMatrix(x, A, n);
    if(solve(A, n, c, y))
        std::cout<<"No se pudo resolver sistema\n";
    
    freeMatrix(A);
}

double evaluaPolinomio(double x, double *c, int n)
{
    int i;
    double ans = 0.0;
    iterate(i, n)
    {
        ans += ( c[i] * pow(x, i) );
    }
    
    return ans;
}


double evaluaLagrangePol(double *X, int n, int i, double x)
{
    if(x == X[i])
        return 1.0;
    
    int j;
    n++;
    double prod = 1.0;
    iterate(j, n)
    {
        if(j != i)
            prod *= ( (x - X[j]) / (X[i] - X[j]) );
    }
    return prod;
}


double interpolaLagrange(double x, double * X, double * Y, int n)
{
    int j;
    double suma = 0.0;
    iterate(j, n)
    {
        suma += Y[j] * evaluaLagrangePol(X, n-1, j, x); // como n es el tam del arreglo, tiene elementos de 0 a n-1
    }
    return suma;
}


double NewtonBasisPolynomial(double *X, int i, double x)
{
    if( i == 0 )
        return 1.0;
    
    int j;
    double prod = 1.0;
    
    iterate(j, i)
    {
        prod *= (x-X[j]);
    }
    return prod;
}


void NewtonDividedDifference(double *X, double *Y, double *a, int n)
{
    double ** A = createMatrix(n, n);
    
    int i, j;
    iterate(i, n)
    {
        A[i][0] = Y[i];
    }
    
    for(j=1; j<n; j++)
    {
        for(i=0; i<n-j; i++)
            A[i][j] = (A[i+1][j-1] - A[i][j-1]) / (X[i+j] - X[i]);
    }
    
    iterate(i, n)
    a[i] = A[0][i];
    
    freeMatrix(A);
}

double interpolaNewton(double x, double *X, double *Y, int n)
{
    double *a = new double[n];
    NewtonDividedDifference(X, Y, a, n);
    
    int i;
    double sum = 0.0;
    iterate(i, n)
    {
        sum += a[i]*NewtonBasisPolynomial(X, i, x);
    }
    
    delete [] a;
    return sum;
}

void interpolaNewton(double *x, double *y, int N, double *X, double *Y, int n)
{
    double *a = new double[n];
    NewtonDividedDifference(X, Y, a, n);
    
    int i, j;
    double sum;
    iterate(j, N)
    {
        sum = 0.0;
        iterate(i, n)
        {
            sum += a[i]*NewtonBasisPolynomial(X, i, x[j]);
        }
        y[j] = sum;
    }
    
    delete [] a;
}

void interpolaSplineCubico(double *x, double *y, int n, double *X, double *Y, int m)
{
    // Pide memoria para vector de los coeficientes Mi donde M0 = Mn = 0
    double * M = new double[n];
    // Llama a la función para calcular los coeficientes Mi
    coeficientes_splines_cubicos(x, y, n, M);
    
    
    // Evalua los polinomios con los puntos de X
    int i, j=0;
    double A, B, C;
    
    A = (M[j+1] - M[j]) / ( 6.0*(x[j+1]-x[j]) );
    B = M[j] / 2.0;
    C = - (x[j+1]-x[j]) * M[j+1] / 6.0 - (x[j+1]-x[j]) * M[j] /3.0 + (y[j+1] - y[j]) / (x[j+1]-x[j]);
    
    j++;
    
    iterate(i, m)
    {
        if(j>n)
            break;
        
        if(X[i] <= x[j])
        {
            j--;
            Y[i] = y[j] + (X[i]-x[j]) * ( C + (X[i] - x[j])*(B + (X[i] - x[j]) * A) );
            j++;
        }
        else
        {
            A = (M[j+1] - M[j]) / ( 6.0*(x[j+1]-x[j]) );
            B = M[j] / 2.0;
            C = - (x[j+1]-x[j]) * M[j+1] / 6.0 - (x[j+1]-x[j]) * M[j] /3.0 + (y[j+1] - y[j]) / (x[j+1]-x[j]);
            j++;
            i--;
        }
        
    }
    
    delete [] M;
}

void coeficientes_splines_cubicos(double *x, double *y, int n, double *M)
{
    n--;    // Ahora n no es el tamaño, sino que x tiene indices 0, 1, ..., n.
    double * h = new double[n];
    h--;                            // En cada uno recorremos uno hacia arriba para que corresponda indice 1 a indice 0
    double * mu = new double[n];
    mu--;
    double * lambda = new double[n];
    lambda--;
    double * d = new double[n];
    d--;
    
    double *c = new double[n+1];
    
    for(int i = 1; i<=n; i++)
    {
        h[i] = x[i] - x[i-1];
    }
    
    for(int i = 1; i<=n-1; i++)
    {
        mu[i] = h[i]/(h[i]+h[i+1]);
        lambda[i] = h[i+1]/(h[i]+h[i+1]);
        d[i] = ( 6.0/(h[i]+h[i+1]) ) * ( (y[i+1]-y[i])/h[i+1] - (y[i]-y[i-1])/h[i] );
        c[i] = 2.0;
    }
    c[0] = c[n] = 2.0;
    M[0] = 0;
    M[n] = 0;
    
    M++;
    
    // Regresando los indices
    mu++;
    lambda++;
    d++;
    h++;
    
    // Resolver para M (1, ..., n-1) tridiagonal con diagonales u, c y l.
    solveTridiagonalMatrix(mu, c, lambda, n-1, d, M);
    
    
    delete [] h;
    delete [] mu;
    delete [] lambda;
    delete [] d;
    delete [] c;
}


void calcula_aproximacion(double *phi, int m, double *x, double *y, int n, double lambda)
{
    // Se asume x esta ordenado
    double min = x[0], max = x[n-1];
    double dp = (max - min) / (double)m;
    
    double * l = new double[m];
    double * d = new double[m+1];
    double * u = new double[m];
    double * b = new double[m+1];
    
    // Llena matriz y vector de terminos independientes
    d[0] = 0.0;
    b[0] = 0.0;
    int i = 0;
    for(int k = 0; k<m; k++)
    {
        double z = min + dp * (double)k;
        double zp = min + dp * (double)(k+1);
        
        d[k] += lambda/(zp-z);
        d[k+1] = lambda/(zp-z);
        
        l[k] = -lambda/(zp-z);
        
        b[k+1] = 0.0;
        
        // Suma para z <= xi < zp
        for( ; i<n; i++)
        {
            if(x[i] > zp || x[i] < z)
                break;
            d[k] += pow((zp - x[i])/(zp-z), 2);
            d[k+1] += pow((x[i]-z)/(zp-z), 2);
            
            l[k] += (zp - x[i])*(x[i] - z) / pow((zp-z), 2);
            
            b[k] += y[i] * ( (zp-x[i])/(zp-z) );
            b[k+1] += y[i] * ( (x[i]-z)/(zp-z) );
        }
        
        u[k] = l[k];
        
    }
    // Llena el ultimo renglon de d y b
    solveTridiagonalMatrix(l, d, u, m+1, b, phi);
    
    delete [] l;
    delete [] d;
    delete [] u;
    delete [] b;
}


void interpolacionElementosFinitos(double *X, double *Y, int N, double *x, double *y, int n, int m, double lambda)
{
    double *phi = new double[m+1];
    lambda = 0.5*lambda;    //Usando otra escala de lambda
    // Ordena x y y en orden ascendente en x
    
    // Obten los coeficientes phi
    calcula_aproximacion(phi, m, x, y, n, lambda);
    
    // Evalua los puntos
    double min = x[0], max = x[n-1];
    double dp = (max - min) / (double)m;
    int j=0;
    for(int k=0; k<m; k++)
    {
        double z = min + dp * (double)k;
        double zp = min + dp * (double)(k+1);
        
        while(X[j]<=zp && j<N)
        {
            Y[j] = phi[k] * (zp-X[j])/(zp-z) + phi[k+1] * (X[j]-z)/(zp-z);
            j++;
        }
    }
    
}


double rombergIntegration(double (*func)(double), double a, double b, int n)
{
    n++;
    double ** dp = createMatrix(n, n);
    double ans;
    
    dp[0][0] = 0.5 * (b-a)*(func(a) + func(b));
    
    
    for (int i=1; i<n; i++)
    {
        dp[i][0] = nextTrapeziumRule(dp[i-1][0], i, a, b, func);
    }
    
    
    for(int j = 1; j<n; j++)
        for(int i = j; i<n; i++)
            dp[i][j] = dp[i][j-1] + (dp[i][j-1] - dp[i-1][j-1]) / (pow(4.0, j) - 1.0);
    
    ans = dp[n-1][n-1];
    
    freeMatrix(dp);
    
    return ans;
}


double nextTrapeziumRule(double previous, int n, double a, double b, double (*func)(double))
{
    double ans;
    
    double h = (b-a)/(double)(1<<n);
    double sum = 0.0;
    for(int k = 1; k<=1<<(n-1); k++)
        sum += func(a + (2.0*k-1.0)*h );
    
    ans = 0.5* previous + h*sum;
    
    return ans;
}

void saveTxt(double *a1, double *a2, int n, const char *cfile)
{
    FILE   *fp=fopen(cfile, "wt");
    int     i;
    
    if(!fp)
    {
        printf("No se puede abrir el archivo\n");
        exit(0);
    }
    printf("Generando el archivo %s ...\n", cfile);
    
    for(i=0; i<n; ++i)
    {
        fprintf(fp, "%.8e ", a1[i]);
        fprintf(fp, "%.8e ", a2[i]);
        fprintf(fp, "\n");
    }
    
    fclose(fp);
}


/** GSL */



gsl_matrix *readMatrix(char *cfile) {
    gsl_matrix *mat;
    int         nr, nc;
    FILE       *f1 = fopen(cfile, "rb");
    
    if(!f1) return(NULL);
    fread(&nr, sizeof(int), 1, f1);
    fread(&nc, sizeof(int), 1, f1);
    mat = gsl_matrix_alloc(nr, nc);
    gsl_matrix_fread(f1, mat);
    fclose(f1);
    return(mat);
}

int writeMatrix(gsl_matrix *mat, char *cfile) {
    FILE       *f1 = fopen(cfile, "wb");
    
    if(!f1) return(1);
    fwrite(&(mat->size1), sizeof(int), 1, f1);
    fwrite(&(mat->size2), sizeof(int), 1, f1);
    gsl_matrix_fwrite(f1, mat);
    fclose(f1);
    return(0);
}

void printMatrixGSL(gsl_matrix *mat) {
    unsigned int i, j;
    for(i=0; i<mat->size1; ++i) {
        for(j=0; j<mat->size2; ++j)
            printf("% 6.2f   ", gsl_matrix_get(mat, i, j));
        printf("\n");
    }
}

void scaleMatrixGSL(gsl_matrix *mat, double dft) {
    unsigned int i, j;
    double   dval;
    for(i=0; i<mat->size1; ++i)
        for(j=0; j<mat->size2; ++j) {
            dval = gsl_matrix_get(mat, i, j);
            gsl_matrix_set(mat, i, j, dft*dval);
        }
}

double **gslmat2array2d(gsl_matrix *mat) {
    double **array2d;
    unsigned int i;
    
    array2d = (double **) malloc(mat->size1 * sizeof(double *));
    if(array2d==NULL) return(NULL);
    
    for(i=0; i<mat->size1; ++i)
        array2d[i] = mat->data + mat->size2*i;
    return(array2d);
}

void printVectorGSL(gsl_vector *vec) {
    unsigned int i;
    double *ptr = vec->data;
    for(i=0; i<vec->size; ++i)
        printf("% 6.2f   ", ptr[i]);
    printf("\n");
}


gsl_vector *readVector(char *cfile) {
    gsl_vector *vec;
    int         nr;
    FILE       *f1 = fopen(cfile, "rb");
    
    if(!f1) return(NULL);
    fread(&nr, sizeof(int), 1, f1);
    vec = gsl_vector_alloc(nr);
    gsl_vector_fread(f1, vec);
    fclose(f1);
    return(vec);
}

void initGslVector(gsl_vector   *vec, double *v, int dim) {
    gsl_block *block=vec->block;
    
    vec->data   = v;
    vec->size   = dim;
    vec->owner  = 0;
    vec->stride = 0;
    if(block!=NULL) {
        block->data = v;
        block->size = dim;
    }
}

void initGslMatrix(gsl_matrix  *mat, double **M, int nr, int nc) {
    gsl_block *block=mat->block;
    
    mat->data   = M[0];
    mat->tda    = nc;
    mat->size1  = nr;
    mat->size2  = nc;
    mat->owner  = 0;
    
    if(block != NULL) {
        block->data = M[0];
        block->size = nr*nc;
    }
}



/** Algoritmos de optimización */



double ** cholesky(double **M, int n)
{
    double ** L = createMatrix(n, n);
    
    clearMatrix(L, n, n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            
            if (j == i) // summation for diagnols
            {
                for (int k = 0; k < j; k++)
                    sum += pow(L[j][k], 2);
                if(M[j][j] - sum <= 0)
                {
                    freeMatrix(L);
                    return NULL;    // Si no es definida positiva y no se pudo hacer la factorización
                }
                L[j][j] = sqrt(M[j][j] - sum);
                
            }
            else {
                
                // Evaluating L(i, j) using L(j, j)
                for (int k = 0; k < j; k++)
                    sum += (L[i][k] * L[j][k]);
                L[i][j] = (M[i][j] - sum) / L[j][j];
            }
        }
    }
    
    return L;
}


void traspose(double ** M, int n)
{
    for(int i=0; i<n; i++)
        for(int j=i; j<n; j++) {
            double aux = M[i][j];
            M[i][j] = M[j][i];
            M[j][i] = aux;
        }
    
}

void transpose2(double ** O, double ** D, int n, int m)
{
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++) {
            D[j][i] = O[i][j];
        }
}

bool inversa(double ** M, double ** Mi, int n)
{
    double ** L = cholesky(M, n);
    
    if (L == NULL) {
        std::cout<<" no se pudo cholesky\n";
        return 1;
    }
    for (int i=0; i<n; i++) {
        if (L[i][i] == 0) {
            return 1;
        }
    }
    
    double ** Li = createMatrix(n, n);
    
    for (int j=0; j<n; j++) {
        for (int i=0; i<n; i++) {
            if(i<j)
                Li[i][j] = 0.0;
            else if (i == j)
                Li[i][j] = 1.0/L[i][i];
            else {
                double sum = 0.0;
                for (int k=j; k<i; k++)
                    sum += L[i][k] * Li[k][j];
                Li[i][j] = -sum/L[i][i];
            }
        }
    }
    
//    // Borrar esto
//
//    std::cout<<"\n\n";
//    productoMatriz(Li, L, Mi, n);
//    printMatrix(Mi, n, n);
//    std::cout<<"\n\n";
    
    copyMatrix(Li, L, n, n);
    traspose(L, n); // Ahora L es Li^T
    productoMatriz(L, Li, Mi, n);
    
    return 0;
}



void scaleVector(double s, double * v, int n)
{
    for (int i=0; i<n; i++)
        v[i] = s*v[i];
}


bool inversaCholesky(double ** M, double ** Mi, int n)
{
    double ** L = createMatrix(n, n);
    copyMatrix(M, L, n, n);
    
    for (int i=0; i<n; i++) {
        if (L[i][i] == 0) {
            return 1;
        }
    }
    
    double ** Li = createMatrix(n, n);
    
    for (int j=0; j<n; j++) {
        for (int i=0; i<n; i++) {
            if(i<j)
                Li[i][j] = 0.0;
            else if (i == j)
                Li[i][j] = 1.0/L[i][i];
            else {
                double sum = 0.0;
                for (int k=j; k<i; k++)
                    sum += L[i][k] * Li[k][j];
                Li[i][j] = -sum/L[i][i];
            }
        }
    }
    
    //    // Borrar esto
    //
    //    std::cout<<"\n\n";
    //    productoMatriz(Li, L, Mi, n);
    //    printMatrix(Mi, n, n);
    //    std::cout<<"\n\n";
    
    copyMatrix(Li, L, n, n);
    traspose(L, n); // Ahora L es Li^T
    productoMatriz(L, Li, Mi, n);
    
    return 0;
}


void vectorSum(double *v1, double *v2, double *r, int n)
{
    for (int i=0; i<n; i++) {
        r[i] = v1[i] + v2[i];
    }
}

void vectorDif(double *v1, double *v2, double *r, int n)
{
    for (int i=0; i<n; i++) {
        r[i] = v1[i] - v2[i];
    }
}


void squareMatrixSum(double ** A, double ** B, double ** R, int n)
{
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            R[i][j] = A[i][j] + B[i][j];
}


void scaleMatrix(double s, double **M, int r, int c)
{
    for (int i=0; i<r; i++) {
        for (int j=0; j<c; j++)
            M[i][j] *= s;
    }
}
