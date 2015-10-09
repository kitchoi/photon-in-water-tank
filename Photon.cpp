#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cstdlib> // for random number generation, and main program termination
#include <string>
#include "SingleRay.h"

#define PI 3.14159265358979323846

double Gaussian(){
	// This is a Gaussian Distribution with center at x=380
	double sigma=17.4847;
//	double sigma=8.;
	double x,y,r2,result;
	do{
		x=-1+2*rand()%30000/30000.;			//random number in [0,1)
		y=-1*2*rand()%30000/30000.;
		r2=x*x+y*y;
	}while(r2>1.0 || r2==0.);

	result=380.+sigma * y * sqrt (-2.0 * log (r2) / r2);
	if (result < 300.) return Gaussian();
	if (result > 600.) return Gaussian();
	return result;
}

double Spectrum(){
	ifstream inFile("../reference/spectrum.txt");
	if (!inFile.good()){
		cerr<<endl<<":::::::Failed to open ../reference/spectrum.txt for wavelength generation";
		exit(1);
	}
	double ran=rand()%30000/30000.;			//random number in [0,1)
	double xLower; 					//lower x value of the bin
	double xUpper; 					//upper x value of the bin
	double yLower; 					//y value at x=xLower
	double yUpper; 					//y value at x=xUpper
	double Prob=0; 					//cumulative probablility of getting x<xUpper
	double ProbLower=0; 				//cumulative probability of getting x<xLower
	double slope;					//useful at the end of spectrum
	double width;					//useful at the end of spectrum
	
	inFile >> xLower >> yLower;

	cout << "Random number =" << ran << endl;
	
	while(true){
		inFile >> xUpper >> yUpper;

		if (inFile.eof()) break;

		ProbLower=Prob;
		Prob+=(yUpper+yLower)/2.*(xUpper-xLower);
		
		if (Prob<0.) {
			cout<< "The spectrum is unphysical."<<endl;
			return -1.;
		}		

		if (Prob>1.) {
			cout<< "The spectrum maybe not normalized. It has to be normalized."<<endl;
			return -1; 	
		}


		if(ran<Prob){
			double delta=(ran-ProbLower)/(Prob-ProbLower);
			return xLower+delta*(xUpper-xLower);
		}
                if (xUpper!=xLower) {
			slope=(yUpper-yLower)/(xUpper-xLower);
			width=xUpper-xLower;
		}
		xLower=xUpper;
		yLower=yUpper;
	}

	inFile.close();

	//End of spectrum
	// Do extrapolation at the end of spectrum
	//	double xmax=xUpper-yUpper/slope;
	//	cout<<"End of Spectrum"<<endl;
	//	cout<<"slope= "<<slope<<endl;
	//	cout<<"xmax= "<<xmax<<endl;
	//	cout<<"Prob= "<<Prob<<endl;
	//	return xmax+yUpper/slope*sqrt((1-ran)/(1-Prob));
	
	cout<<"End of Spectrum"<<endl;
	double x=xUpper;
	while(true){
		x+=width/100.;
		double area=(x-xUpper)*yUpper+0.5*slope*(x-xUpper)*(x-xUpper);
		if (ran-Prob < area) return x;
	}
}

int main(){

/////// Setting begins

	string AbsorptionModel="Fry";	//"Fry" for fry.dat or "SuperK" for superk.dat
	bool include_undiffract=true;	//include un undiffracted photon as reference or not

	double photonStep=1.;		//in meter
	double initialT=0.0;		//in ns
	double initialX=0.0;		//in meter
	double initialZ=0.0;		//in meter
	double TankDiameter=50.;	//in meter
	double TankHeight=60.;		//in meter

	// For the test photon
	bool absorption=true;		//apply absorption of photon or not
	bool rayleigh=true;		//apply rayleigh scattering or not
	bool tempGradient=true;		//apply temperature gradient or not, if not Temp=TempTop
	bool presGradient=true;		//apply pressure gradient or not, if not Pressure=PressureTop
	double wavelength;		//please define in the for-loop
	double TempTop=273.15+15.;	//in Kelvin
	double TempBottom=273.15+5.;	//in Kelvin
	double PressureTop=101325.;	//in Pa

	// For the reference photon only
	bool ref_absorption=false;
	bool ref_rayleigh=false;
	bool ref_tempGradient=false;
	bool ref_presGradient=false;
	double ref_wavelength=380.;	//fixed wavelength
	double ref_TempTop=273.15+10.;	//in Kelvin
	double ref_TempBottom=275.15+10.;	//in Kelvin
	double ref_PressureTop=101325.+1000.*9.8*30.;	//in Pa, though named Top, the reference photon is living at this pressure

/////// Setting Ends

	int number_absorbed=0;


        // Open output file

        FILE * outputfile;
        outputfile = fopen("output.dat","w");
        if (include_undiffract) {
		fprintf(outputfile,"# theta wavelength deltaX deltaZ deltaT deltaD X1(no diffract) Z1 T1 X2(has diffract) Z2 T2 lambda scatter\n");
	}
	else{
		fprintf(outputfile,"# theta wavelength X2(has diffract) Z2 T2 lambda\n");
	}


	// Do the iteration of theta, ranging from -pi/2 to pi/2

	for (int i=500;i<=500;i++)	//iteration of angle
	{
		for (int j=0;j<1000;j++) //iteration of photon
		{
		double initialTheta=-PI/2.+PI*i/1000.;


////////////////Select how wavelength is generated///////////////
		wavelength=Spectrum();				//
//		wavelength=380.;				//
//		wavelength=Gaussian();				//
/////////////////////////////////////////////////////////////////

		cout<<endl;
		cout<<"*******************************************"<<endl;
		cout<<"InitialTheta="<<initialTheta/PI<<" PI"<<endl;
		cout<<"InitialX="<<initialX<<"m"<<endl;
		cout<<"InitialZ="<<initialZ<<"m"<<endl;
		cout<<"Test Photon Wavelength="<<wavelength<<" nanometer"<<endl;
		cout<<"Step="<<photonStep<<"m"<<endl;
		cout<<"*******************************************"<<endl;
	
	//
	//Run the test photon
	//
		SingleRay photon2(photonStep,initialT,initialX,initialZ,wavelength,initialTheta,\
				TempTop,TempBottom,PressureTop,TankDiameter,TankHeight,\
				absorption,rayleigh,tempGradient,presGradient,AbsorptionModel);

		//photon2 stores the data of a diffracted ray
	        photon2.WithDiffract();

		if (photon2.absorbed) {		//photon absorbed
			number_absorbed++;	//record number of absorbed photons
			continue;		//run the next photon
			}

		cout << "WithDiffract: x="<<photon2.x<<", z="<<photon2.z<<", t="<<photon2.t<<endl;

		if (!include_undiffract) {
			fprintf (outputfile, "%3.3f %3.3f %3.5f %3.5f %3.5f %3.5f \n",\
			initialTheta/PI,wavelength,\
			photon2.x,photon2.z,photon2.t,photon2.lambda);
			}


	//
	//Run an undiffracted photon for reference
	//
		double deltaX,deltaZ,deltaT,deltaD;
		if(include_undiffract){

			SingleRay photon1(photonStep,initialT,initialX,initialZ,ref_wavelength,initialTheta, \
				ref_TempTop,ref_TempBottom,ref_PressureTop,TankDiameter,TankHeight, \
				ref_absorption,ref_rayleigh,ref_tempGradient,ref_presGradient,AbsorptionModel);

			// photon1 stores the data of an undiffracted ray
			photon1.WithoutDiffract();
			
			cout << "WithoutDiffract: x="<<photon1.x<<", z="<<photon1.z<<", t="<<photon1.t<<endl;
			// Calculate the difference in x,z and t
			deltaX=photon2.x-photon1.x;
			deltaZ=photon2.z-photon1.z;
			deltaT=photon2.t-photon1.t;
			cout << "deltaX="<<deltaX<<"m, deltaZ="<<deltaZ<<"m, deltaT="<<deltaT<<"ns"<<endl;

			// Calculate sqrt(x^2+z^2)
			deltaD=sqrt(deltaX*deltaX+deltaZ*deltaZ);

			//Output to data to file
			fprintf (outputfile, "%3.3f %3.3f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.3f %d\n",\
			initialTheta/PI,wavelength,deltaX,deltaZ,\
			deltaT,deltaD,photon1.x,photon1.z,\
			photon1.t,photon2.x,photon2.z,photon2.t,
			photon2.lambda,photon2.scattered);
			}
		}
	}
	fclose(outputfile);
	cout << "Number of absorbed photons: "<< number_absorbed <<endl;
}
