#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>

using namespace std;

int main(int argc, char const *argv[])
{
    default_random_engine e;
    int count;
    double media_distribuzione, dev_std;
    string nome_file = "";

    cout << "Quanti numeri casuali desidera generare?";
    cin >> count;
    cout << "Quale desidera sia la media?";
    cin >> media_distribuzione;
    cout << "Quale desidera sia la deviazione standard?";
    cin >> dev_std;
    cout << "Che nome desidera assegnare al file di output?";
    cin >> nome_file;
    
    ofstream outputFile(nome_file);
   
    normal_distribution<> norm(media_distribuzione, dev_std);
    for (int i = 0; i < count; i++)
    {
        outputFile << norm(e) << endl;
    }
    

    return 0;
}
