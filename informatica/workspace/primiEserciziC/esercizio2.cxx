#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    int a = 2048;
    int &ref = a;
    int *puntatore = &a;

    cout << ref << endl;
    cout << &ref << endl;
    cout << puntatore << endl;
    cout << &a << endl;
    return 0;
}

