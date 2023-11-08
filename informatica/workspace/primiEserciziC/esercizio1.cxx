using namespace std;
#include <iostream>

int main(int argc, char const *argv[])
{
    //due numeri in input, stampa somma e differenza
    
    int a, b, somma, diff;
    
    cout << "Inserire un numero a" << endl;
    cin >> a;
    cout << "Inserire un numero b" << endl;
    cin >> b;

    somma = a + b;
    diff = a - b;

    cout << "La loro somma: " << somma << endl;
    cout << "La loro differenza: " << diff << endl;

    //input base e altezza triangolo, calcolare area
    float base, altezza, area;
    
    cout << "Inserire base" << endl;
    cin >> base;
    cout << "Inserire altezza" << endl;
    cin >> altezza;

    area = base * altezza / 2;
    
    cout << "Area del triangolo: " << area << endl;
    
    //input x y z, calcolare area laterale e volume
    float x, y, z, areaL, volume;
    
    cout << "Inserire larghezza" << endl;
    cin >> x;
    cout << "Inserire lungezza" << endl;
    cin >> y;
    cout << "Inserire profondità" << endl;
    cin >> z;

    volume = x * y * z;
    //areaL = 2 * x * y + 
    
    cout << "Volume: " << volume << endl;
    
    return 0;
}
