#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

double const L = 10;
int const N = 1000;
double const STEP = L/N;
double const A = 3;
double const B = 7;
double const V0 = 1;
double const K = 1;


double numerov(double E, double y_1, double y0, double c_1, double c_0, double c1, double step=STEP);
double bisection(double a, double b, double s, double (*potential)(double x));
double f(double x, double (*potential)(double x), double step=STEP);
double *f_stamp(double x, double (*potential)(double x), ofstream &out_file, double step=STEP);
double step_potential(double x);
double flatHole_potential(double x);
double harmonic_potential(double x);
double simpson_integration(double* f, double step, int n);

int main(int argc, char const *argv[]){
    cout.precision(10);
    double const Emin = 0;
    double const Emax = 10;
    double const dE = 0.001;
    double const s = 1e-8;
    short int choice;
    double (*potential)(double);
    ofstream out_file;

    out_file.open("eigenEnergies.txt");

    int nEGrid = (Emax - Emin)/dE + 1;

    double* energyGrid = new double[nEGrid];
    vector<double> eigenEnergy;

    cout << "Scegli il potenziale: " << endl;
    cout << "   1: buca piatta" <<endl;
    cout << "   2: buca a gradino" << endl;
    cout << "   3: armonico" << endl;
    cin >> choice;

    switch (choice)
    {
    case 1:
        potential = &flatHole_potential;
        break;
    case 2:
        potential = &step_potential;
        break;
    case 3: 
        potential = &harmonic_potential;
        break;
    
    default:
        return 1;
    }

    for (int i = 0; i < nEGrid; i++){
        double E = dE * i;//candidato autovalore energia
        energyGrid[i] = f(E, potential);
    }   
    //ora col metodo della bisezione trovo tutti gli autovalori 
    for (int i = 1; i < nEGrid; i++){
        //se cambia segno -> avvio bisezione
        if (energyGrid[i]*energyGrid[i+1]<=0){
            eigenEnergy.push_back(bisection(dE*i, dE*(i+1), s, potential));
            out_file << setprecision(10) << fixed  <<  eigenEnergy.back() << endl;        
        }    
    }
    delete[] energyGrid;
    //calcolo autofunzioni
    for (int i = 0; i < eigenEnergy.size(); i++){
        out_file.close();
        out_file.open("armoniche/armonica"+to_string(i)+".txt");
        f_stamp(eigenEnergy.at(i), potential, out_file);
    }
    
    return 0;
}

double f(double x, double (*potential)(double x), double step){
    double y0 = 1;
    double y_1 = 0;
    double y1;
    for (int j = 1; j <= N; j++){
        y1 = numerov(-2.0*x, y_1, y0, -2.0*potential((j-1)*step),-2.0*potential(j*step),-2.0*potential((j+1)*step), step); 
        y_1 = y0;
        y0 = y1;
    }
    return y0;
}

double* f_stamp(double x, double (*potential)(double x), ofstream &out_file, double step){
    double *function = new double[N+1];
    double *squared_funtion = new double[N+1];
    double norm = 0;
    function[0] = 0;
    function[1] = 1;
    squared_funtion[0] = 0;
    squared_funtion[1] = 1;

    for (int j = 2; j <= N; j++){
        function[j] = numerov(-2.0*x, function[j-2], function[j-1], -2.0*potential((j-1)*step),-2.0*potential(j*step),-2.0*potential((j+1)*step), step); 
        squared_funtion[j]=pow(function[j],2);
    }
    norm = sqrt(simpson_integration(squared_funtion, step, N));
    for (int i = 0; i <=N; i++){
        function[i]/=norm;
        out_file << i*step << ";" << function[i] << endl;
    }
    return function;   
}

double bisection(double a, double b, double s, double (*potential)(double x)){
    double c;
    do{
        c = (b+a)*0.5;
        if((f(c, potential)*f(a, potential)) < 0)
            b = c;
        else
            a = c;
    }while((b-a) > s);
    return c;
}

double numerov(double E, double y_1, double y0, double c_1, double c_0, double c1, double step){
    double y1 = ((E-c_1)*step*step - 12.0) * y_1;
    y1 += ((E-c_0) * 10.0 * step * step + 24.0)*y0;
    y1 /= (12.0 - (E-c1) * step * step);
    return y1;
}

double step_potential(double x){
    if(x>A &&  x<B)
        return V0;
    return 0;
}

double flatHole_potential(double x){
    return 0;
}

double harmonic_potential(double x){
    return K * pow((x-L/2.0), 2);
}

double simpson_integration(double* f, double step, int n){
    //metodo simpson
    double simp = f[0] + f[n];
    for (int i = 1; i < n; i++)
        simp+=2*(i%2 + 1)*f[i];   
    simp*=step/3;
    return simp;
}