// we use SI units
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>

double const kB = 1.38064852e-23; // Boltzmann constant
double const A = 1e-10;           // Angstrom

void confiniziale(int N, double l, double m, double T, double **r, double **v);
double maxwell_boltzmann(double m, double T, double v);
void energia_temp(double eps, double sigma, double m, int n, double l, double **r0, double **v0, double &ene, double &t);
void energia_potenziale(double eps, double sigma, double m, int n, double l, double **r0, double &ene);
void delta_energia_potenziale(double eps, double sigma, double m, int n, double l, int ia, double **r0, double **r1, double &dene);
void forza(double eps, double sigma, double m, int n, double l, double **r0, double **f);
double drand(void);

using namespace std;

int main(int argc, const char *argv[])
{
    double rho, l, eps, sigma; // densità, lato cella di simulazione, parametri LJ

    double const A = 1e-10;        // Angstrom
    double m = 39.95 * 1.6747e-27; // mass atomo Argon
    int N;                         // numero di atomi
    long int Ntr;                  // numero time steps di rilassamento
    long int Nt;                   // numero time steps di simualazione
    double dt;                     // time step
    double *r0mat, *r1mat, *v0mat, *v1mat, *f0mat, *f1mat, *rDmat, *rD0mat;
    double **r0, **r1, **v0, **v1, **f0, **f1, **rD, **rD0; // array bidimensionali posizioni e velocità
    double T;                                               // temperatura
    int Nstampa;                                            // intervallo di passi per stampa e analisi
    double smax;                                            // massimo spostamento lungo una direzione Cartesiana
    int ia;                                                 // indice random atomo da spostare
    double sr;                                              // spostamento random
    double delta_e;                                         // differenza energetica

    double ene, temp;
    char cstampa;
    string namexyz;
    ofstream fxyz, flast, fene;

    double ene0, ene1;

    eps = 120 * kB;
    sigma = 3.4 * A;

    cout << "Programma avviato : \n";
    /*
    cout << "Numero di atomi: \n";
    cin >> N;
    cout << "Densità: (kg/m^3) \n";
    cin >> rho;
    cout << "Tempertura (K): \n";
    cin >> T;
    cout << "Spostamento massimo (Angstrom): \n";
    cin >> smax;
    smax*=A;
    cout << "Numero di passi rilassamento: \n";
    cin >> Ntr;
    cout << "Numero di passi simulazione: \n";
    cin >> Nt;
    cout << "Stampa ogni : \n";
    cin >> Nstampa;
    cout << "Scrivi xyz file (t/f) : \n";
    cin >> cstampa;
    */

    N = 864;
    rho = 1374;
    T = 94.4;
    smax = 5.;
    smax*=A;
    Ntr = 2E05;
    Nt = 1E06;
    Nstampa = 100;
    cstampa = 'f';

    if (cstampa == 't')
    {
        cout << "Nome file xyz : \n";
        cin >> namexyz;
        fxyz.open(namexyz, ios::out);
    }

    // trova lato scatola di simulazione
    l = pow(m * N / rho, 1. / 3.);

    cout << "Lato cella Simulazione (m): " << l << '\n';
    // alloca
    r0mat = new double[N * 3];
    r1mat = new double[N * 3];
    v0mat = new double[N * 3];
    v1mat = new double[N * 3];
    f0mat = new double[N * 3];
    f1mat = new double[N * 3];
    rDmat = new double[N * 3];
    rD0mat = new double[N * 3];
    r0 = new double *[N];
    r1 = new double *[N];
    v0 = new double *[N];
    v1 = new double *[N];
    f0 = new double *[N];
    f1 = new double *[N];
    rD = new double *[N];
    rD0 = new double *[N];
    for (int i = 0; i < N; i++)
    {
        r0[i] = &r0mat[i * 3];
        r1[i] = &r1mat[i * 3];
        v0[i] = &v0mat[i * 3];
        v1[i] = &v1mat[i * 3];
        f0[i] = &f0mat[i * 3];
        f1[i] = &f1mat[i * 3];
        rD[i] = &rDmat[i * 3];
        rD0[i] = &rD0mat[i * 3];
    }

    fene.open("energia.dat", ios::out);
    
    // starting coordinartes
    cout << "Genera configurazione di partenza \n :";
    confiniziale(N, l, m, T, r0, v0);

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            rD[i][k] = r0[i][k];
            rD0[i][k] = r0[i][k];
            r1[i][k] = r0[i][k];
        }
    }
    cout << "Calcolo energia potenziale \n";
    energia_potenziale(eps, sigma, m, N, l, r0, ene0);
    fene << scientific << ene0 << '\n';
    fene.flush();
    cout <<  "Inizio Monte Carlo\n";
    // Monte Carlo
    for (long it = 0; it < Ntr + Nt; it++)
    {
        // scelgo atomo da spostare
        ia = round(drand() * (N - 1));
        // sposto random stando attento a PBC
        for (int j = 0; j < 3; j++)
        {
            rD[ia][j] = rD[ia][j] + (drand() * 2.0 * -1.0) * smax;
            if (rD[ia][j] > l)
                rD[ia][j] -= l;
            else if (rD[ia][j] < 0)
                rD[ia][j] += l;
        }
        // calcolo differenza energia
        delta_energia_potenziale(eps, sigma, m, N, l, ia, rD0, rD, delta_e);
        // vedo se accettare o no
        if (delta_e < 0)
        {
            // accetto xtrial
            for (int j = 0; j < 3; j++)
            {
                rD0[ia][j] = rD[ia][j];
                r1[ia][j] = rD[ia][j];
            }
        }
        else
        {
            double rand = drand();
            if (exp(-delta_e / (kB * T)) > rand)
            {
                // accetto xtrial
                for (int j = 0; j < 3; j++)
                {
                    rD0[ia][j] = rD[ia][j];
                    r1[ia][j] = rD[ia][j];
                }
            }
        }
        // aggiorno configurazione ed energia, messa in ene0
        energia_potenziale(eps, sigma, m, N, l, r1, ene0);
        fene << it << "  " << ene0 << '\n';

        if (it > Ntr && it % Nstampa == 0)
        {
            if (cstampa == 't')
            {
                fxyz << N << '\n';
                fxyz << "Argon\n";
                for (int i = 0; i < N; i++)
                {
                    fxyz << "Ar " << r0[i][0] / A << "  " << r0[i][1] / A << "  " << r0[i][2] / A << '\n';
                }
            }
        }
    }

    // sciave su disco l'ultimo snapshot in unità di Angstrom
    flast.open("argonlast.txt", ios::out);
    flast << N << "  " << l / A << '\n';
    for (int i = 0; i < N; i++)
    {
        flast << r0[i][0] / A << "  " << r0[i][1] / A << "  " << r0[i][2] / A << '\n';
    }
    flast.close();

    // dealloca
    delete[] r0mat;
    delete[] r1mat;
    delete[] v0mat;
    delete[] v1mat;
    delete[] rDmat;
    delete[] rD0mat;
    delete[] r0;
    delete[] r1;
    delete[] v0;
    delete[] v1;
    delete[] f0;
    delete[] f1;
    delete[] rD;
    delete[] rD0;

    if (cstampa == 't')
        fxyz.close();
}

void confiniziale(int N, double l, double m, double T, double **r0, double **v0)
{
    int i;
    double vp;   // velocità più probabile
    double fmax; // massimo della funzione f
    double f, c, v, d;
    int ni;
    ni = (int)pow((double)N, 1. / 3.);
    if (ni * ni * ni < N)
        ni++;

    vp = sqrt(2. * kB * T / m);
    fmax = maxwell_boltzmann(m, T, vp);
    i = 0;
    for (int ix = 0; ix < ni; ix++)
    {
        for (int iy = 0; iy < ni; iy++)
        {
            for (int iz = 0; iz < ni; iz++)
            {
                if (i < N)
                {
                    r0[i][0] = l / ((double)ni) * ix;
                    r0[i][1] = l / ((double)ni) * iy;
                    r0[i][2] = l / ((double)ni) * iz;
                    do
                    {
                        v = 4. * vp * drand();
                        c = fmax * drand();
                        f = maxwell_boltzmann(m, T, v);

                    } while (c > f);
                    // prima la direzione
                    d = 0.;
                    for (int j = 0; j < 3; j++)
                    {
                        v0[i][j] = drand() - 0.5;
                        d += pow(v0[i][j], 2.);
                    }
                    d = sqrt(d);
                    for (int j = 0; j < 3; j++)
                        v0[i][j] *= v / d;
                    i++;
                }
            }
        }
    }
    // metto  a zero la velocità del centro di massa
    // tutte le masse sono uguali
    double vcm[3] = {0., 0., 0.};
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            vcm[k] += v0[i][k];
        }
    }
    for (int k = 0; k < 3; k++)
        vcm[k] /= ((double)N);
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            v0[i][k] -= vcm[k];
        }
    }
}

double maxwell_boltzmann(double m, double T, double v)
{
    double f;
    f = sqrt(2. / M_PI) * pow(m / (kB * T), 3. / 2.) * pow(v, 2.) * exp(-m * pow(v, 2.) / (2. * kB * T));
    return f;
}

void energia_temp(double eps, double sigma, double m, int n, double l, double **r0, double **v0, double &ene, double &t)
{
    double ekin = 0;

    ene = 0;
    ekin = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 3; j++)
            ekin += 0.5 * m * pow(v0[i][j], 2.);
    }
    t = 2. * ekin / (double(3 * n - 3)) / kB;

    double epot = 0.;
    double dist_min, dist;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {

            dist_min = 1e10;
            for (int kx = -1; kx < 2; kx++)
            {
                for (int ky = -1; ky < 2; ky++)
                {
                    for (int kz = -1; kz < 2; kz++)
                    {
                        dist = pow(r0[i][0] - r0[j][0] + l * kx, 2.) + pow(r0[i][1] - r0[j][1] + l * ky, 2.) + pow(r0[i][2] - r0[j][2] + l * kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min)
                            dist_min = dist;
                    }
                }
            }
            //       cout << "dist " << dist_min << "  l " << l << '\n';
            epot += 4. * eps * (pow(sigma / dist_min, 12.) - pow(sigma / dist_min, 6.));
        }
    }
    ene = ekin + epot;
}

void energia_potenziale(double eps, double sigma, double m, int n, double l, double **r0, double &ene)
{

    ene = 0;

    double epot = 0.;
    double dist_min, dist;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {

            dist_min = 1e10;
            for (int kx = -1; kx < 2; kx++)
            {
                for (int ky = -1; ky < 2; ky++)
                {
                    for (int kz = -1; kz < 2; kz++)
                    {
                        dist = pow(r0[i][0] - r0[j][0] + l * kx, 2.) + pow(r0[i][1] - r0[j][1] + l * ky, 2.) + pow(r0[i][2] - r0[j][2] + l * kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min)
                            dist_min = dist;
                    }
                }
            }
            //cout << "dist " << dist_min << "  l " << l << '\n';
            epot += 4. * eps * (pow(sigma / dist_min, 12.) - pow(sigma / dist_min, 6.));
        }
    }
}

void delta_energia_potenziale(double eps, double sigma, double m, int n, double l, int ia, double **r0, double **r1, double &dene)
{

    dene = 0;

    double depot = 0.;
    double dist_min, dist;

    for (int j = 0; j < n; j++)
    {
        if (j != ia)
        {
            dist_min = 1e10;
            for (int kx = -1; kx < 2; kx++)
            {
                for (int ky = -1; ky < 2; ky++)
                {
                    for (int kz = -1; kz < 2; kz++)
                    {
                        dist = pow(r0[ia][0] - r0[j][0] + l * kx, 2.) + pow(r0[ia][1] - r0[j][1] + l * ky, 2.) + pow(r0[ia][2] - r0[j][2] + l * kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min)
                            dist_min = dist;
                    }
                }
            }
            depot -= 2. * eps * (pow(sigma / dist_min, 12.) - pow(sigma / dist_min, 6.));
        }
    }
    for (int j = 0; j < n; j++)
    {
        if (j != ia)
        {
            dist_min = 1e10;
            for (int kx = -1; kx < 2; kx++)
            {
                for (int ky = -1; ky < 2; ky++)
                {
                    for (int kz = -1; kz < 2; kz++)
                    {
                        dist = pow(r1[ia][0] - r1[j][0] + l * kx, 2.) + pow(r1[ia][1] - r1[j][1] + l * ky, 2.) + pow(r1[ia][2] - r1[j][2] + l * kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min)
                            dist_min = dist;
                    }
                }
            }
            depot += 2. * eps * (pow(sigma / dist_min, 12.) - pow(sigma / dist_min, 6.));
        }
    }

    dene = depot;
}

double drand(void)
{
    double c;
    c = ((double)rand()) / ((double)RAND_MAX);
    return c;
}

void forza(double eps, double sigma, double m, int n, double l, double **r0, double **f)
{
    int kk[3];
    double fact;
    double dist, dist_min;
    double ff;
    double *p;
    p = &f[0][0];
    for (int k = 0; k < n * 3; k++)
    {
        *p = 0.;
        p++;
    }

    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < i; j++)
        {

            dist_min = 1e10;
            for (int kx = -1; kx < 2; kx++)
            {
                for (int ky = -1; ky < 2; ky++)
                {
                    for (int kz = -1; kz < 2; kz++)
                    {
                        dist = pow(r0[i][0] - r0[j][0] + l * kx, 2.) + pow(r0[i][1] - r0[j][1] + l * ky, 2.) + pow(r0[i][2] - r0[j][2] + l * kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min)
                        {
                            dist_min = dist;
                            kk[0] = kx;
                            kk[1] = ky;
                            kk[2] = kz;
                        }
                    }
                }
            }
            fact = 24. * eps / pow(dist_min, 2.) * (2. * pow(sigma / dist_min, 12.) - pow(sigma / dist_min, 6.));
            for (int k = 0; k < 3; k++)
            {
                ff = fact * (r0[j][k] - (r0[i][k] + kk[k] * l));
                f[i][k] -= ff;
                f[j][k] += ff;
            }
        }
    }
}
