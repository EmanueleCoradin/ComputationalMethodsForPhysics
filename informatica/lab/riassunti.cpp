#include <string>
using namespace std;

struct pippo{
	int a;
	int b;
	int c;
	float d;
	bool e;
	bool f;
	string g;
	int h[10];
}

//prototipo
double media(int*);

double media(int* vettore){
	double m = 0;
	int i = 0;
	for(int n: vettore){
		m+=n;
		i++;	
	}
	m/=i;
	return m;
}
