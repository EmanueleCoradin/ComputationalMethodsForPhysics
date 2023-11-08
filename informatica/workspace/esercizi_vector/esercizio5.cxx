#include <iostream>
#include <vector>
#include <random>

using namespace std;

int main(int argc, char const *argv[])
{
    uniform_int_distribution<unsigned> u(18,30);
    default_random_engine e;

    vector <int> frequenzeVoti(13, 0);
    int n, voto;
    cout << "Quanti voti desidera generare? ";
    cin >> n;
    for (int i = 0; i < n; i++)
    {
        voto = u(e);
        frequenzeVoti.at(30-voto)++;
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
