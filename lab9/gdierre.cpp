//calcola la pair correlation function
//usa unità di Angstrom


#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>

double gaussiana(double x, double x0, double sigma);

int main(int argc, const char * argv[])
{
    int N;//numero atomi
    double L;//lato cella di simulazione
    double *r0mat;
    double **r0;
    double Rmax, sigma;
    double *gr,*xr;
    double dist,dist_min;
    int ns;
    
    std::ifstream flast;
    std::ofstream fgr;
    
    flast.open("argonlast.txt", std::ios::in);
    flast >> N;
    flast >> L;
    std::cout << "Atomi: " << N << "  Lato (A): " << L << '\n';
    r0mat= new double[N*3];
    r0= new double*[N];
    for (int i=0; i<N; i++) {
        r0[i]=&r0mat[i*3];
    }
    for (int i=0; i<N; i++) {
        for(int k=0;k<3;k++) flast >> r0[i][k];
    }
    flast.close();
    
    
    std::cout << "Rmax (A): \n";
    std::cin >> Rmax;
    std::cout << "Broadening gaussiano (A): \n";
    std::cin >> sigma;
    std::cout << "Numero intervalli: \n";
    std::cin >> ns;
    
    gr = new double[ns];
    xr = new double[ns];
    for(int l=0;l<ns;l++) gr[l]=0.;
    for(int l=0;l<ns;l++) xr[l]=Rmax/((double) ns)*((double) l);
    
    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if(j!=i){
              dist_min=1e10;
              for(int kx=-1;kx<2;kx++){
                  for(int ky=-1;ky<2;ky++){
                      for(int kz=-1;kz<2;kz++){
                          dist=pow(r0[i][0]-r0[j][0]+L*kx,2.)+pow(r0[i][1]-r0[j][1]+L*ky,2.)+pow(r0[i][2]-r0[j][2]+L*kz,2.);
                          dist=sqrt(dist);
                          if(dist<dist_min) dist_min=dist;
                      }
                  }
              }
              for(int l=0;l<ns;l++) gr[l]+=gaussiana( xr[l], dist_min,  sigma);
            }
   
        }
    
    }
    
    for(int l=0;l<ns;l++) gr[l]*=(1./(4.*M_PI*pow(xr[l]*N,2.)))*pow(L,3.);
    
    fgr.open("gidierre.txt", std::ios::out);
    for(int l=0;l<ns;l++){
        fgr << xr[l] << "  "  << gr[l] << '\n';
    }
    fgr.close();
    
    //dealloca
    delete [] r0mat;
    delete [] r0;
    delete [] gr;
    delete [] xr;
    
}


double gaussiana(double x, double x0, double sigma){
    double g;
    g=(1./(sigma*sqrt(2.*M_PI)))*exp(-pow((x-x0)/sigma,2.)/2.);
    return g;
}
