#include <iostream>
#include <math.h>

using namespace std;

void numeriImmaginari(){
 double re, im, modulo, fase;
    struct complex {double re;  double im;} z1, z2, somma, prodotto;
    
    cout << "dato numero immmaginario, calcolo modulo e fase"<<endl<<endl;

    cout << "Inserire parte reale: " << endl;
    cin >> re;

    cout << "Inserire parte immmaginaria " << endl;
    cin >> im; 

    z1.re = re;
    z1.im = im;
    //calcolare modulo e fase numero complesso
    modulo = sqrt(re*re + im*im);
    fase = atan(im/re);

    cout << "Modulo: " << modulo << endl;
    cout << "Fase: " << fase << " radianti" << endl;

    //calcolare somma e prodotto di due numeri complessi
    cout << "inserendo un secondo complesso, calcolo la somma e il prodotto con il primo"<<endl<<endl;
    
    cout << "Inserire parte reale: " << endl;
    cin >> re;

    cout << "Inserire parte immmaginaria " << endl;
    cin >> im; 

    z2.re = re;
    z2.im = im;

    somma.re = z1.re + z2.re;
    somma.im = z1.im + z2.im;
    prodotto.re = z1.re * z2.re - z1.im * z2.im;
    prodotto.im = z1.im * z2.re + z1.re * z2.im;
    
    cout << "La somma vale: " << somma.re << " + i" << somma.im << ". " << endl;
    cout << "Il prodotto vale: " << prodotto.re << " + i" << prodotto.im << ". " << endl;

}
void polinomiSecondo(){
    double a,b,c, delta, x1, x2;
    //inserimento parametri da terminale
    cout << "Inserire parametri equazione secondo grado, in output gli zeri" << endl;
    cout << "Inserire a ";
    cin >> a;
    cout << "Inserire b ";
    cin >> b;
    cout << "Inserire c ";
    cin >> c;
    //calcolo 
    delta = b * b - 4.0 * a * c;
    if(delta > 0){
        x1 = (-b + sqrt(delta)) / (2 * a);
        x2 = (-b - sqrt(delta)) / (2 * a);
        cout << "Gli zeri sono : x1 = " << x1 << ", x2 = " << x2 << endl;
    }
    else if(delta == 0){
        x1 = (-b) / (2 * a);
        cout << "Lo zero e' : x1 = " << x1 << endl;
    }
    else{
        cout << "Il polinomio non presenta zeri (delta negativo)"<<endl;
    }

}

void mediaNumeri(){
    double media, i, a;
    i = -1;
    media = 0;
    a = 0;
    
    do
    {
        media+=a;
        i++;
        cout << "Inserire un numero ";
    } while (cin >> a);
    
    media/=i;
    cout << endl << "La media dei valori inseriti: " << media << endl; 
    
}


void mediaNumeri2(){
    double media, i, a;
    i = 0;
    media = 0;
    a = 0;
    
    cout << "Inserire un numero ";
    while (cin >> a)
    {
        media+=a;
        i++;
        cout << "Inserire un numero ";
    }
    
    media/=i;
    cout << endl << "La media dei valori inseriti: " << media << endl; 
    
}
int main(int argc, char const *argv[])
{   
    //numeriImmaginari();
    //polinomiSecondo();
    mediaNumeri2();
    return 0;
}
