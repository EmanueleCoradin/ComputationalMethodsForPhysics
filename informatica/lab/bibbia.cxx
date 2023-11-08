#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

struct studente
    {
        string nome;
        int matricola;
        int scritto_inf;
        int scritto_stat;
        int rel;
        int orale;
        int votoFinale;
};

double mediaPesata(studente &s){
    double media = s.scritto_inf / 4.0 + s.scritto_stat / 4.0 + s.rel / 5.0 + s.orale * 3.0 / 10.0;
    return media;
}

int main(int argc, char const *argv[])
{
    /*
    string nomeFile = "";
    cout << "Inserire nome del file ";
    cin >> nomeFile;
    ifstream inputFile(nomeFile);
    ofstream outputFile("out.txt");

    string linea= "";
    string letto;
    while(getline(inputFile, letto))
        linea+=letto;
    outputFile << "Il file " << nomeFile << " contiene " << linea.size() << " caratteri";

    studente matte;
    matte.nome = "Matteo";
    matte.matricola = 12;
    matte.orale =  25;
    matte.rel = 30;
    matte.scritto_inf = 26;
    matte.scritto_stat = 27;

    cout << "La media di " << matte.nome << " è " << mediaPesata(matte);
    
    double media = 27.0, devstd = 1.5;

    uniform_int_distribution<> randI (0,9);
    uniform_real_distribution<> randR (0,1);
    normal_distribution<> norm (media, devstd);

    default_random_engine e;
    e.seed(25);
    int intero = randI(e);
    double reale = randR(e);
    double gauss = norm(e);

    cout << intero << " " << endl << reale << " " << endl << " "  << gauss;

    vector<int> casuale;
    for (int i = 0; i < 50; i++)
    {
        casuale.push_back(randI(e));
        cout << casuale.at(i) << endl;
    }
    string s;
    cin >> s;
    sort(casuale.begin(), casuale.end());
    vector<int>::iterator walk = casuale.begin();
    while(walk!= casuale.end()){
        cout << * walk <<endl;
        walk++;
    }
    */
    int numero = 18;
    int &ref = numero;
    int *ptr = &numero;
    cout << numero << ptr << ref << endl;
    cout << &numero << &ptr << &ref << endl;

    return 0;
}
