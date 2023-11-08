#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    string in_file;
    string out_file_name;
    ofstream out_file;
    double inizio, fine, step;
    //numero punti
    int n;
    //Array di punti
    double* x;
    double* f;
    
    //lettura dati  
    cin >> n;
    cin >> inizio;
    cin >> fine;
    cin >> out_file_name;  

    out_file.open(out_file_name);
    x = (double*) new double[n];
    f = (double*) new double[n];
    step = (fine-inizio)/n;

    //riempio l'array di punti
    for (int i = 0; i < n; i++)
    {
        x[i] = inizio + step * i;
        //calcolo il seno
        f[i] = sin(x[i]);   
        //scrivo il risultato
        out_file << x[i] << " " << f[i] << endl; 
    }
    
    delete []x;
    delete []f;
    
    return 0;
}
