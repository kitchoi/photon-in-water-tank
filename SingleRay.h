#include <math.h>
#include <string>
#include "RefIndex.h"
#include <fstream>
#include <iostream>
#include <cstdlib> //for random number generation, and terminating program from functions
using namespace std;

class SingleRay{
 public:
	SingleRay(double step,double t,double x,double z,double lambda,double theta,double TopTemp,double BottomTemp,double TopPressure,double Diameter,double Height,bool absorption_turn_on,bool rayleigh_turn_on,bool tempGradient, bool presGradient, string AbsorptionModel);
	bool Absorption(double wavelength, double distance, double ran);
	bool Rayleigh(double wavelength, double distance, double ran);
	double CurrentTemp(double z);
	double CurrentPressure(double z);
	double RayleighAngle();
	void WithDiffract();
	void WithoutDiffract();
	double step; //step of numerical calculation
	double t; // time for the photon after emission, in nanosecond
	double x; // radial position of photon, zero at tank center
	double z; // vertical position of photon, zero at tank center
	double lambda; //wavelength of photon, in micrometer
	double theta; // angle of elevation, in radian
	double TopTemp; //Temperature of tank top
	double BottomTemp; //Temperature of tank  bottom, in Kelvin
	double TopPressure; //Pressure of tank top, in Pascal
	double Diameter; //Diameter of tank, in meter
	double Height; //Height of the tank, in meter
	bool absorbed; //if the photon is absorbed
	bool absorption_turn_on;
	bool rayleigh_turn_on;
	bool tempGradient;
	bool presGradient;
	string AbsorptionModel; //"Fry" or "SuperK"
	int scattered; //0: not scattered, 1:scattered at least once
};

