#include <cmath>
#include <iostream>

using namespace std;

double function(double x){
    return sin(x);
}

double* inversione(double* h, double* d, double* b, int n);

int main(int argc, char const *argv[])
{
    //spline cubica naturale
    //--------- DEFINIZIONE -----------
    double* x, *f, *g, *h, *b, *d, *z, *alfa, *beta, *gamma, *eta, *p, *x_graf;
    int n, n2;
    double inizio, fine, step, step2;

    //assegnazione
    n = 5;
    n2 = 100;
    inizio = 0.0;
    fine = 2.0*M_PI;
    step = (fine-inizio)/(n);

    x = (double*) new double[n+1];
    f = (double*) new double[n+1];
    g = (double*) new double[n];
    h = (double*) new double[n];
    b = (double*) new double[n];
    d = (double*) new double[n];
    z = (double*) new double[n];
    alfa = (double*) new double[n];
    beta = (double*) new double[n];
    gamma = (double*) new double[n];
    eta = (double*) new double[n];
    p = (double*) new double[(n+1)*(n2+1)];
    x_graf = (double*) new double[(n+1)*(n2+1)];

    //inizializzazione griglia
    for (int i = 0; i < n+1; i++){
        x[i] = inizio + step * i;
        f[i] = function(x[i]);
    }

    //calcolo parametri
    for (int i = 0; i < n; i++){
       g[i] = f[i+1] - f[i];
       h[i] = x[i+1] - x[i];
    }

    //calcolo vettore b
    for (int i = 1; i < n; i++){
        b[i] = (g[i]/h[i] - g[i-1]/h[i-1])*6.0;
        d[i] = 2.0*(h[i-1]+h[i]);
    }
    
    //calcolo vettore derivata seconda
    z = inversione(h,d,b,n);
    
    //Calcolo coefficienti
    for (int i = 0; i < n; i++){
        alfa[i] = z[i+1]/(6.0*h[i]);
        beta[i] = - z[i]/6.0/h[i];
        gamma[i] = -pow(h[i],2) * alfa[i] + f[i+1]/h[i];
        eta[i] = -pow(h[i], 2) * beta[i] - f[i]/h[i];
    }

    //calcolo funzione interpolante
    cout << "x;y;func" << endl;
    for (int i = 0; i < n+1; i++){
        step2 = h[i]/n2;
        for (int j = 0; j < n2; j++){  
            int w;
            double delta, delta2;
            w = n2*i+j;
            x_graf[w] = i*step+j*step2;
            delta = x_graf[w]-x[i];
            delta2 = x_graf[w] - x[i+1];

            p[w] = alfa[i]*pow((delta),3) + beta[i] * pow((delta2),3) + gamma[i] * delta+ eta[i] * delta2;
            cout << x_graf[w] << ";" << p[w]  << ";" << f[i] << endl;
        }
    }
    return 0;
}

double* inversione(double* h, double* d, double* b, int n){
    double *v, *w, *t, *y, *z; 
    v = (double*) new double[n];
    w = (double*) new double[n];
    t = (double*) new double[n];
    y = (double*) new double[n];
    z = (double*) new double[n];
    //primo passo iterativo
    t[0] = 0;
    v[0] = 0;
    v[1] = h[1];
    w[1] = d[1];
    y[1] = b[1]/w[1];

    for (int i = 1; i < n+1; i++){
        v[i] = h[i];
        w[i] = d[i] - t[i-1]*v[i-1];
        t[i] = h[i]/w[i];
    }
    
    for (int i = 2; i < n+1; i++)
        y[i] = (b[i] - v[i-1] * y[i-1])/w[i];
    
    //iterazioni per z
    z[n-1] = y[n-1];
    for (int i = n-2; i >=0; i--)
        z[i] = y[i] - t[i]*z[i+1];

    return z;
}