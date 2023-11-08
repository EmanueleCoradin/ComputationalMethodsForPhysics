#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

long double *difference(long double *v1, long double *v2);
long double magnitude(long double *v);
long double *force(long double **positions, long double *masses, int w);
long double *mul(long double *v, long double s);
long double *sum(long double *v1, long double *v2);

long double const G = 6.67*pow(10,-17);
int n; // numero di corpi

struct planet
{
    string name;
    long double mass;
    long double *positions[3];
    long double *velocity[3];
};

int main(int argc, char const *argv[])
{
    struct planet
    {
        string name;
        long double mass;
        long double *position;
        long double *velocity;
    };

    long timestep = 24 * 60 * 60;
    long unsigned tempo = 0;

    long double *f, *f1, *util;          // vettore forza
    long double **positions, *cm, *v_cm; // vettore delle posizioni
    long double **velocities;            // vettore delle velocità
    long double **r1;                    // vettore delle posizioni
    long double **v1;                    // vettore delle velocità
    long double *masses, *mat, mtot;

    mtot = 0;

    ofstream *out_file = new ofstream[10];
    ifstream fin("sistema.dat");
    for (int i = 0; i < 10; i++)
    {
        out_file[i].open("object" + to_string(i) + ".dat");
    }

    fin >> n;

    positions = new long double *[n];
    velocities = new long double *[n];
    masses = new long double[n];
    f = new long double[3];
    f1 = new long double[3];
    util = new long double[3];
    cm = new long double[3];
    v_cm = new long double[3];

    mat = new long double[3 * n];
    for (int i = 0; i < n; i++)
        positions[i] = &mat[3 * i];
    mat = new long double[3 * n];
    for (int i = 0; i < n; i++)
        velocities[i] = &mat[3 * i];

    struct planet *system = new struct planet[n];

    // importazione condizioni iniziali
    for (int i = 0; i < n; i++)
    {
        fin >> system[i].mass;
        system[i].position = new long double[3];
        system[i].velocity = new long double[3];
        for (int j = 0; j < 3; j++)
            fin >> system[i].position[j];
        for (int j = 0; j < 3; j++)
            fin >> system[i].velocity[j];
        positions[i] = system[i].position;
        velocities[i] = system[i].velocity;
        masses[i] = system[i].mass;
    }

    // calcolo ed elimino velocità cm
    for (int i = 0; i < 3; i++)
    {
        cm[i] = 0;
        v_cm[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        cm = sum(cm, mul(positions[i], masses[i]));
        v_cm = sum(v_cm, mul(velocities[i], masses[i]));
        mtot += masses[i];
    }

    cm = mul(cm, 1.0 / mtot);
    v_cm = mul(v_cm, 1.0 / mtot);

/*
    for (int i = 0; i < n; i++)
    {
        positions[i] = difference(positions[i], cm);
        out_file[i] << 0 << " " << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;
        velocities[i] = difference(velocities[i], v_cm);
    }
  */  

    for (int j = 0; j < 1000; j++)
    {
        tempo += timestep*j;
        cout << j << endl;
        r1 = new long double *[n];
        v1 = new long double *[n];

        mat = new long double[3 * n];
        for (int i = 0; i < n; i++)
            r1[i] = &mat[3 * i];
        mat = new long double[3 * n];
        for (int i = 0; i < n; i++)
            v1[i] = &mat[3 * i];
        // propago con velocity Verlet
        // calcolo posizioni
        for (int i = 0; i < n; i++)
        {
            f = force(positions, masses, i);
            r1[i] = mul(velocities[i], timestep);
            long double numero = timestep * timestep / 2.0 / masses[i];
            util = mul(f, numero);
            r1[i] = sum(r1[i], util);
            r1[i] = sum(r1[i], positions[i]);
            out_file[i] << tempo << " " << r1[i][0] << " " << r1[i][1] << " " << r1[i][2] << endl;
        }
        // calcolo velocità
        for (int i = 0; i < n; i++)
        {
            f1 = force(r1, masses, i);
            v1[i] = sum(f, f1);
            v1[i] = mul(v1[i], timestep / 2 / masses[i]);
            v1[i] = sum(v1[i], velocities[i]);
        }
        positions = r1;
        velocities = v1;
    }
    return 0;
}

long double *force(long double **positions, long double *masses, int w)
{
    long double *f = new long double[3];
    long double *d = new long double[3];
    long double m = 0;

    for (int i = 0; i < 3; i++)
    {
        f[i] = 0;
        d[i] = 0;
    }

    for (int i = 0; i < w; i++)
    {
        d = difference(positions[i], positions[w]); // versore della forza
        m = 1.0 / pow(magnitude(d), 3);
        m *= G * masses[i] * masses[w]; // modulo della forza
        f = sum(f, mul(d, m));
    }

    for (int i = w + 1; i < n; i++)
    {
        d = difference(positions[i], positions[w]);
        m = 1.0 / pow(magnitude(d), 3);
        m *= G * masses[i] * masses[w];
        f = sum(f, mul(d, m));
    }

    return f;
}

long double *sum(long double *v1, long double *v2)
{
    long double *sum = new long double[3];
    for (int i = 0; i < 3; i++)
        sum[i] = v1[i] + v2[i];
    return sum;
}

long double *difference(long double *v1, long double *v2)
{
    long double *diff = new long double[3];
    for (int i = 0; i < 3; i++)
        diff[i] = v1[i] - v2[i];

    return diff;
}

long double magnitude(long double *v)
{
    long double m = 0;
    for (int i = 0; i < 3; i++)
        m += v[i] * v[i];
    return sqrt(m);
}

long double *mul(long double *v, long double s)
{
    long double *v2 = new long double[3];
    for (int i = 0; i < 3; i++)
        v2[i] = v[i] * s;
    return v2;
}
