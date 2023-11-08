#include <cmath>
#include <iostream>

using namespace std;

double omega, l;
const double G = 9.80665;
int N,M;//dimensione spazio fasi e numero di dati da salvare

double newton_law(double **y, double t=0);
void Eulero_esplicito(double *y0, double*y1, double timestep);
void Runge_kutta(double *y0, double*y1, double timestep);

int main(int argc, char const *argv[]){
    //----------------- DEFINIZIONE -------------------
    double theta_0, v_0, periodo_app, periodo_eulero, periodo_runge, timestep, timetot;
    double *mat, *f;
    double **y;
    
    //----------------- ASSEGNAZIONE -----------------
    N = 2;
    M = 2;  
    l = 1;//lunghezza del pendolo
    theta_0 = -M_PI/6;//angolo iniziale
    v_0 = 0;//velocità angolare iniziale
    timestep = 0.000001;
    timetot = 0;

    omega = sqrt(G/l);

    mat = new double[N*M];
    y = new double*[N];
    
    for (int i = 0; i < N; i++)
        y[i] = & mat[M*i];

    f = new double[N];    

    //----------------- ANALISI -----------------
    periodo_app = 2*M_PI/omega;
    cout << "Periodo di piccola oscillazione: " << periodo_app << endl;
    
    //-------- METODO EULERO ESPLICITO --------
    //dato iniziale
    y[1][0] = theta_0;
    y[1][1] = v_0;

    do{
        //resetto dato iniziale
        y[0][0] = y[1][0];
        y[0][1] = y[1][1];

        Eulero_esplicito(y[0], y[1], timestep);
        //aggiungo timestep
        timetot+=timestep;
    }while (y[1][0]<0);
    //interpolo intercetta
    periodo_eulero = timetot - timestep - y[0][0] * timestep / (y[1][0]-y[0][0]);
    periodo_eulero*=4;
    
    cout << "Periodo Eulero: " << periodo_eulero << endl;

    //-------- METODO RUNGE KUTTA --------
    timetot = 0;
    y[1][0] = theta_0;
    y[1][1] = v_0;
    
    do{
        //resetto dato iniziale
        y[0][0] = y[1][0];
        y[0][1] = y[1][1];
        
        Runge_kutta(y[0], y[1], timestep);

        timetot+=timestep;
    }while (y[1][0]<0);
    
    //interpolo intercetta
    periodo_runge = timetot - timestep - y[0][0] * timestep / (y[1][0]-y[0][0]);
    periodo_runge*=4;
    
    cout << "Periodo Runge-Kutta: " << periodo_runge << endl;

    return 0;
}

double newton_law(double *y, double t=0){
    return - omega * omega* sin(y[0]);
}

void Eulero_esplicito(double *y0, double*y1, double timestep){
    y1[0] = y0[0] + y0[1] * timestep;
    y1[1] = y0[1] + newton_law(y0) * timestep;
}

void Runge_kutta(double *y0, double *y1, double timestep){
    double *mat;
    double **Y, **F;

    Y = new double*[4];
    F = new double*[4];

    mat = new double[4*N];

    for (int i = 0; i < 4*N; i++)
        Y[i] = &mat[N*i];

    mat = new double[4*N];
    for (int i = 0; i < 4*N; i++)
        F[i] = &mat[N*i];
 
    Y[0][0] = y0[0];
    Y[0][1] = y0[1];

    F[0][0] = Y[0][1];
    F[0][1] = newton_law(Y[0]);

    Y[1][0] = y0[0] + F[0][0]*timestep/2;
    Y[1][1] = y0[1] + F[0][1]*timestep/2;

    F[1][0] = Y[1][1];
    F[1][1] = newton_law(Y[1]);

    Y[2][0] = y0[0] + F[1][0]*timestep/2;
    Y[2][1] = y0[1] + F[1][1]*timestep/2;

    F[2][0] = Y[2][1];
    F[2][1] = newton_law(Y[2]);

    Y[3][0] = y0[0] + F[2][0]*timestep;
    Y[3][1] = y0[1] + F[2][1]*timestep;

    F[3][0] = Y[3][1];
    F[3][1] = newton_law(Y[3]);

    for (int i = 0; i < N; i++)
        y1[i] = y0[i] + (F[0][i] + 2*F[1][i] + 2*F[2][i] + F[3][i])*timestep/6;
}