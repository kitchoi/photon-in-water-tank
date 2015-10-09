// This program applies the formula given by Schiebener et al. to calculate the index of refraction of water at different temperature and pressure
// Variables with underscores replace those with stars in the formula

//Unit of lambda (wavelength) : nanometer
//Unit of pressure: Pa
//Unit of temperature: Kelvin

// The following line includes the FreeSteam software for calculating the water density
//#include "/home/user/freesteam-0.8.1/steamcalculator.h"

#include <math.h>
#include <iostream>
#include "RefIndex.h"
using namespace std;

RefIndex::RefIndex(double inputp, double inputt, double inputl){
  this->p=inputp;
  this->t=inputt;
  this->lambda=inputl;
}

double RefIndex::nfind() {
	double a[8]={0.243905091,9.53518094e-3,-3.64358110e-3,2.65666426e-4,1.59189325e-3,2.45733798e-3,0.897478251,-1.63066183e-2};

// Calculate the water density

// The followings are codes for applying the FreeSteam software
/*	SteamCalculator SC;	
	Temperature T1 = t * Kelvin;
	Pressure P1 = p * Pascal;
	SC.set_pT(P1,T1);
	rho_=SC.dens().getRawValue()/1000.;
*/

#define atm 101325.
#define celcius(t) ((t)-273.15)

// The following formula is fit from values generated from FreeSteam
// C.I.=99% for pressure with 1-7atm and temperature at 1-18 degree Celcius
	rho_=999.825+0.0511962*p/atm-9.14698e-6*p*p/atm/atm+0.05100779*celcius(t) \
		-0.00681379*celcius(t)*celcius(t)-0.000261832*p/atm*celcius(t);
	rho_=rho_/1000.;

// for debugging
//	cout << "Temperature = " << t << " K, Pressure = "<< p << " Pa, Density= " << rho_*1000. << "kg/m3" << endl;

// Calculate the refractive index
	double lamb_=lambda/589.;
	double t_=t/273.15;
	double lamb_2=pow(lamb_,2);
	double lambuv_2=0.2292020*0.2292020;
	double lambir_2=5.432937*5.432937;

	// solve the quadratic equation
	double F=rho_*(a[0]+a[1]*rho_+a[2]*t_+a[3]*lamb_2*t_+a[4]/lamb_2+a[5]/(lamb_2-lambuv_2)+a[6]/(lamb_2-lambir_2)+a[7]*rho_*rho_);
	double n=sqrt((2*F+1)/(1-F));
	return n;
}
