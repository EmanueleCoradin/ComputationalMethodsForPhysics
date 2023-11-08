#include <cmath>
#include <iostream>

using namespace std;

double f(double x){
    return cos(x)-x;
}

double df(double x){
    return -sin(x)-1;
}

int main(int argc, char const *argv[])
{
    double a, b, c, e;

    cout.precision(10);
    cout << "Metodo della bisezione" << endl;
    cout << "Inserire precisione e: " << endl;
    cin >> e;

    a = 0.0;
    b = 1.0;

    cout << endl << "Elaborazione" << endl;
    do{
        c = (b+a)*0.5;
        if((f(c)*f(a)) < 0)
            b = c;
        else
            a = c;
        cout << c << endl;
    }while((b-a) > e);
    
    cout << "La soluzione è: " << c << endl; 

    //metodo di Newton Raphson
    
    double x_trial, x_sol;
    
    cout << "Metodo di Newton" << endl;
    x_sol = 1.0;
    do
    {
        x_trial = x_sol;
        x_sol = x_trial - f(x_trial)/df(x_trial); 
        cout << x_sol << endl;
    } while (abs(x_trial-x_sol) > e);
    
    cout << endl << "La soluzione è: " << x_sol;
    
    return 0;
}
