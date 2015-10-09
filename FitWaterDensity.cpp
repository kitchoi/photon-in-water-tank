//This program calculate the residue of the fit of water density as a function of temperature and pressure
//The fit (from mathematica): waterdensity=999.825+0.0511962*P-9.14698e-6*P*P+0.0510779*T-0.00681379*T*T-0.000261832*P*T
//in 99% confidence level

#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main(){
	ifstream datafile("./data/waterdensity.dat");
	double p; //pressure
	double t; //temperature
	double rho; //water density
	
	double a[6]={999.825,0.0511962,-9.14698e-6,0.0510779,-0.00681379,-0.000261832};

	while(true){
		datafile >> p >> t >> rho;
		if (datafile.eof()) break;
		cout << a[0]+a[1]*p+a[2]*p*p+a[3]*t+a[4]*t*t+a[5]*p*t-rho <<endl;
	}
return 0;
}
