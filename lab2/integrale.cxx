#include <cmath>
#include <iostream>

using namespace std;

void algorithm(int n)
{
        //Array di punti
    long double* x;
    long double* f;
    long double* f2;
    long double* diff;
    long double inizio, fine, step, res, rett_naif, trap, simp;

    //assegnazione
    inizio = 0;
    fine = 1;
    res = exp(1) - 1;
    rett_naif = 0;
    trap = 0;
    simp = 0;
    step = (fine-inizio)/(n);

    x = (long double*) new long double[n+1];
    f = (long double*) new long double[n+1];

    //riempio l'array di punti
    for (int i = 0; i < n+1; i++)
    {
        x[i] = inizio + step * i;
        f[i] = exp(x[i]);
        
    }

    //metodo rettangoli naif
    for (int i = 0; i < n; i++)
    {
        rett_naif+= f[i];
    }
    rett_naif*=step;
    cout << res - rett_naif;

    //metodo trapezi
    trap = (f[0]+f[n])/2;
    for (int i = 1; i < n; i++){
        trap+=f[i];
    }
    trap*=step;
    cout << " " << res - trap;

    //metodo simpson
    simp = f[0] + f[n];
    for (int i = 1; i < n; i++)
    {
        simp+=2*(i%2 + 1)*f[i];
    }
    
    simp*=step/3;

    cout << " " << res - simp << endl;
    
}

int main(int argc, char const *argv[])
{
    cout << "Rettangoli_naif Trapezi Simpson" << endl;
    for (int i = 10; i < 110; i+=10)
    {
        cout << i << " ";
        algorithm(i);
    }
    for (int i = 100; i < 1001; i+=100)
    {
        cout << i << " ";
        algorithm(i);
    }
    for (int i = 1000; i != 100000; i+=1000)
    {
        cout << i << " ";
        algorithm(i);
    }
    for (int i = 100000; i != 1000000; i+=10000)
    {
        cout << i << " ";
        algorithm(i);
    }
    return 0;
}

