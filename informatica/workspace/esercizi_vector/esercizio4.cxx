#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char const *argv[])
{
    vector <int> voti;
    int singoloVoto;
    cout << "Inserire i voti degli esami superati: " << endl;
    while (cin >> singoloVoto){
        if(singoloVoto > 17 && singoloVoto < 31)        
            voti.push_back(singoloVoto);
    }
    cout << "Sono stati inseriti " << voti.size() << " voti" << endl;
    
    for (int i = 18; i <= 30; i++){
        cout << i << ": ";
        for(auto &v : voti){
            if(v == i){
                cout << "*";
            }
        }
        cout << endl;
    }
    return 0;
}
