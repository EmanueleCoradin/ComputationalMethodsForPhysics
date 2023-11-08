#include <iostream>
#include <vector>
#include <random>

using namespace std;

int main(int argc, char const *argv[])
{
    default_random_engine e;

    vector <int> frequenzeVoti(13, 0);
    int n, voto;
    double media_distribuzione, dev_std;
    cout << "Quanti voti desidera generare? ";
    cin >> n;
    cout << "Inserire la media della distribuzione:  ";
    cin >> media_distribuzione;
    cout << "Inserire la deviazione standard: ";
    cin >> dev_std;

    normal_distribution<> norm(media_distribuzione, dev_std);

    for (int i = 0; i < n; i++)
    {
        voto = lround(norm(e));
        if(voto > 17 && voto < 31)
            frequenzeVoti.at(30-voto)++;
        else i--;
    }
    for (int i = 0; i < frequenzeVoti.size(); i++)
    {
        cout << 30-i << ": ";
        for (int j = 0; j < frequenzeVoti.at(i); j++)
        {
            cout <<  "*";
        }
        cout << endl;
    }
    return 0;
}
