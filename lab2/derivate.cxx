#include <cmath>
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    //Array di punti
    double* x;
    double* f;
    double* df0;
    double* df1;
    double* df2;
    double* df4; 
    double* diff;
    int n;
    double inizio, fine, step;

    //assegnazione
    n = 100;
    inizio = 0;
    fine = 2 * M_PI;
    step = (fine-inizio)/(n-1);

    x = (double*) new double[n];
    f = (double*) new double[n];
    df0 = (double*) new double[n];
    df1 = (double*) new double[n-1];
    df2 = (double*) new double[n-2];
    df4 = (double*) new double[n-4];

    //riempio l'array di punti
    for (int i = 0; i < n; i++)
    {
        x[i] = inizio + step * i;
        f[i] = sin(x[i]);  
        df0[i] = cos(x[i]);
    }

    cout << "Derivata destra" << endl;
    //calcolo derivata destra
    for (int i = 0; i < n-1; i++)
    {
        df1[i] = (f[i+1]-f[i])/step;
        cout << x[i] << "   " << (df0[i] - df1[i]) << endl;
    }
    cout << endl;
    
    cout << "Derivata O(h^2)" << endl;
    //calcolo con O(h^2)
    for (int i = 1; i < n-1; i++)
    {
        df2[i] = (f[i+1]-f[i-1])/(2*step);
        cout << x[i] << "   " << (df0[i] - df2[i]) << endl;
    }
    cout << endl;

    cout << "Derivata O(h^4)" << endl;
    //calcolo con O(h^4)
    for (int i = 2; i < n-2; i++)
    {
        df4[i] = (-f[i+2] + 8*f[i+1]- 8*f[i-1] + f[i-2])/(12*step);
        cout << x[i] << "   " << (df0[i] - df4[i]) << endl;
    }
    
    return 0;
}