#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//Prototipi
double eulero(int N, double l, double dt, double theta0);
double runge(int N, double l, double dt, double theta0);

//Main
int main(){

  float l = 1.0;

  double theta0 = 0.523;

  double dt = 0.001;

  int N = 1000;
  
  int select;

  cin >> select;

  double T;

  switch(select){
    case (1): T = eulero(N,l,dt,theta0); 
    break;
    case (2): T = runge(N,l,dt,theta0);
    break;
  };

  cout << endl << "Periodo: " << T << endl;
  
  return 0;
}


//Metodo di Runge-Kutta
double runge(int N, double l, double dt, double theta0){

  //Contatore di cicli
  int time = 0;

  //creo la matrice y_t
  double *ymat;
  double **y;
  ymat = new double[2*2];
  y = new double*[2];

  for(int i=0; i<2; i++)
    y[i] = &ymat[2*i];

  y[1][0] = theta0;
  y[1][1] = 0.0;

  //Inizializzo gli Y
  double *Y1;
  double *Y2;
  double *Y3;
  double *Y4;
  Y1 = new double[2];
  Y2 = new double[2];
  Y3 = new double[2];
  Y4 = new double[2];

  //tempi di cambio segno
  float tpiu = 0;
  float tmeno = 0;
  
  do{
    y[0][0] = y[1][0];
    y[0][1] = y[1][1];

    Y1[0] = y[0][0];
    Y1[1] = y[0][1];

    Y2[0] = y[0][0]  +  Y1[1]*dt*0.5;
    Y2[1] = y[0][1]  - (9.81/l)*sin(Y1[0])*dt*0.5;

    Y3[0] = y[0][0]  +  Y2[1]*dt*0.5;
    Y3[1] = y[0][1]  - (9.81/l)*sin(Y2[0])*dt*0.5;

    Y4[0] = y[0][0]  +  Y3[1]*dt;
    Y4[1] = y[0][1]  - (9.81/l)*sin(Y3[0])*dt;

    y[1][0] = y[0][0] + (Y1[1]+2.0*Y2[1]+2.0*Y3[1]+Y4[1])*dt/6.0;
    y[1][1] = y[0][1] - (9.81/l)*(sin(Y1[0])+2.0*sin(Y2[0])+2.0*sin(Y3[0])+sin(Y4[0]))*dt/6.0;

    cout << "Angolo: " << y[1][0] << " Velocita': " << y[1][1] << endl;  
    
    time += 1;

  }while (y[1][0]>0);
  
  tpiu = time*dt;
  tmeno = (time-1)*dt;

  double periodo = 2*(tpiu+tmeno);
  return periodo;

};



//Metodo di eulero
double eulero(int N, double l, double dt, double theta0){

  //Contatore di cicli
  int time = 0;

  //creo la matrice y_t
  double *ymat;
  double **y;
  ymat = new double[2*2];
  y = new double*[2];

  for(int i=0; i<2; i++)
    y[i] = &ymat[2*i];

  y[1][0] = theta0;
  y[1][1] = 0.0;

  //tempi di cambio segno
  float tpiu = 0;
  float tmeno = 0;
  
  do{
    y[0][0] = y[1][0];
    y[0][1] = y[1][1];

    y[1][0] = y[0][0] + y[0][1]*dt;
    y[1][1] = y[0][1] - (9.81/l)*sin(y[0][0])*dt;

    cout << "Angolo: " << y[1][0] << " Velocita': " << y[1][1] << endl;  
    
    time += 1;

  }while(y[0][0]*y[1][0]>0);
  
  tpiu = time*dt;
  tmeno = (time-1)*dt;

  double periodo = 2*(tpiu+tmeno);
  return periodo;

};