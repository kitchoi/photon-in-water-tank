#include <math.h>
#include <iostream>
using namespace std;

int main(){
	for (int i=1;i<=1000;i++){
		double wavelength=210.+540*i/1000.;
		cout << wavelength << '\t' << exp(19.648)/pow(wavelength,4.1538) <<endl;
	}
}
