#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;

double media (vector<double> dati, int indice_iniziale, int indice_finale){
    double media = 0;

    for (int i = indice_iniziale; i <= indice_finale; i++)    {
        media+=dati.at(i);
    }
    media/=indice_finale - indice_iniziale;
    return media;
}

double devstd (vector<double> dati, int indice_iniziale, int indice_finale, double media){
    double dev = 0;

    for (int i = indice_iniziale; i <= indice_finale; i++)
    {
        dev+= pow((media - dati.at(i)), 2);
    }
    dev/=indice_finale - indice_iniziale; 
    return sqrt(dev);
}

int main(int argc, char const *argv[])
{   
    string nome_file;
    double valore;
    vector<double> dati;

    cout << "Inserire il nome del file da cui leggere i valori" << endl;
    cin >> nome_file;
    ifstream inputFile(nome_file);
    
    if(!inputFile){
        cout << "errore, file di lettura non esistente";
        return (-1);
    }
    while(inputFile >> valore){
        dati.push_back(valore);
    }

    random_shuffle(dati.begin(), dati.end());
    
    //SET INTERO
    double media_set = media(dati, 0, dati.size()-1);
    double dev_set = devstd(dati, 0, dati.size()-1, media_set);
    //PRIMO QUARTO
    double media_primo = media(dati, 0, dati.size()/4);
    double dev_primo = devstd(dati, 0, dati.size()/4, media_primo);
    //SECONDO QUARTO
    double media_sec = media(dati, dati.size()/4, dati.size()/2-1);
    double dev_sec = devstd(dati, dati.size()/4, dati.size()/2-1, media_sec);
    //TERZO QUARTO
    double media_terzo = media(dati, dati.size()/2, (3 * dati.size()/4)-1);
    double dev_terzo = devstd(dati, dati.size()/2, (3 * dati.size()/4)-1, media_terzo);
    //ULTIMO QUARTO
    double media_quarto = media(dati, 3 * dati.size()/4, dati.size()-1);
    double dev_quarto = devstd(dati, 3 * dati.size()/4, dati.size()-1, media_quarto);
    
    //RESTITUZIONE OUTPUT
    cout << "Set: (" << media_set << " ± " << dev_set << ")" << endl;
    cout << "Primo quarto: (" << media_primo << " ± " << dev_primo << ")" << endl;
    cout << "Secondo quarto: (" << media_sec << " ± " << dev_sec << ")" << endl;
    cout << "Primo quarto: (" << media_terzo << " ± " << dev_terzo << ")" << endl;
    cout << "Primo quarto: (" << media_quarto << " ± " << dev_quarto << ")" << endl;
    
    return 0;
}