#include "SingleRay.h"

#define PI 3.14159265358979323846
#define lightSpeed 0.299792458

//This is a function for linear interpolation
//date file should contain two columns (x y)
double Interpolate(string filename_str, double x0){
	const char* filename=filename_str.c_str();
	ifstream datafile(filename,ifstream::in);
	if (!datafile.is_open()) {
		cerr << ":::::::Error opening "<<filename<<" for interpolation"<<endl;
		exit(1);		//terminate the main program
	}

	double x,xprev,y,yprev,slope;
	datafile >> x >> y;
	xprev=x; yprev=y;
	
	//The input x0 is smaller than the smallest value in the data file
	if(x0<x){
		cerr<<":::::::Out of range error! Cannot interpolate with file "<<filename<<" and x0=" <<x0<<endl;
		exit(1); 
	}
	
	while(true){
		datafile >> x >> y;
		if (datafile.eof()) break;
		if (x>x0) return (y-yprev)*(x0-xprev)/(x-xprev)+yprev;
		if (x!=xprev) slope=(y-yprev)/(x-xprev);
		xprev=x;
		yprev=y;
	}

	datafile.close();

	//End of data
	cerr<<":::::::Out of range error! cannot interpolate with file "<<filename<<" and x0= "<<x0<<endl;
	exit(1);	

	//Extrapolation
	//return y+slope*(x0-x);
}

double SingleRay::CurrentTemp(double z){
	if (!tempGradient) return TopTemp;
	if (tempGradient) return (TopTemp+BottomTemp)/2.+(TopTemp-BottomTemp)*z/Height;
}

double SingleRay::CurrentPressure(double z){
	if (!presGradient) return TopPressure;
	if (presGradient) return TopPressure+1000.*9.8*(Height/2.-z);
}

SingleRay::SingleRay(double step,double t,double x,double z,double lambda,double theta, \
		double TopTemp,double BottomTemp,double TopPressure,double Diameter,double Height, \
		bool absorption_turn_on, bool rayleigh_turn_on, bool tempGradient, bool presGradient, string AbsorptionModel){
	this->step=step;
	this->t=t;
	this->x=x;
	this->z=z;
	this->lambda=lambda;
	this->theta=theta;
	this->TopTemp=TopTemp;
	this->BottomTemp=BottomTemp;
	this->TopPressure=TopPressure;
	this->Diameter=Diameter;
	this->Height=Height;
	this->absorbed=false;
	this->absorption_turn_on=absorption_turn_on;
	this->rayleigh_turn_on=rayleigh_turn_on;
	this->tempGradient=tempGradient;
	this->presGradient=presGradient;
	this->AbsorptionModel=AbsorptionModel;
	this->scattered=0;
}

bool SingleRay::Absorption(double wavelength, double distance, double ran){

	// obtain logged absorption coefficient
	double log_a;
	if (AbsorptionModel.compare("Fry")==0) log_a=Interpolate("fry.dat",wavelength);
	if (AbsorptionModel.compare("SuperK")==0) log_a=Interpolate("superk.dat",wavelength);

	double a=pow(10,log_a);					// restore to absorption coefficient
	if (ran <= step*a) {						// the photon is absorbed
		absorbed=true;
		return true;
	}
	return false;
}

bool SingleRay::Rayleigh(double wavelength, double distance, double ran){
	double b=exp(19.648)/pow(wavelength,4.1538);	//fit from data
	if (ran >= 1-step*b) {
		scattered=1; //change status of photon
		return true;
	}
	return false;
}

#define phasef(u) ((u)*(u)*(u)+3*(u)-4+8*(r))
#define phasef_d(u) (3*(u)*(u)+3)

double SingleRay::RayleighAngle(){

//Scattering phase function is 3(1+cos*cos)/16/pi
//random number = integrating the solid angle with phase from 0 to x
//Equivalent to solve for root of (cube of u)+3u-4+8r=0
//where u = cos(phase)

//Save the phase angle in 3D space
ofstream rayleigh_angle_file("rayleigh_angle.dat",ios_base::app);

	double r=rand()%30000/30000.;
	if (r==0.5) {
		rayleigh_angle_file<<0.5*PI<<endl;
		rayleigh_angle_file.close();
		return 0.5*PI;
	}

	double x;
	//Apply Newton's method
	x=-1.;
	do{
		x=x-phasef(x)/phasef_d(x);
	}while (fabs(phasef(x))>1.e-4);	

	rayleigh_angle_file<<acos(x)<<endl;
	rayleigh_angle_file.close();
	return acos(x);

}

#undef phasef
#undef phasef_d

void SingleRay::WithDiffract(){


	// Calculate the temperature and pressure at the current position
	// Assumed linear variation in temperature
	// Asssumed pressure changes with depth as density*g*h

	double temp=CurrentTemp(z);
	double pressure=CurrentPressure(z);

	// Calculate the refractive index at the current position
	RefIndex refIndex(pressure,temp,lambda);
	double n=refIndex.nfind();
	double nbefore=n;

	while (true)
	{
		if (t<1.e-5) {  
			//first step
			
			x+=step*cos(theta);
			z+=step*sin(theta);
			t+=n*step/lightSpeed;
		}
		else {
			// Not the first step

			double old_theta=theta;
			
			//Calculate the temperature and pressure at the current position
        		temp=CurrentTemp(z);
		        pressure=CurrentPressure(z);
			RefIndex refIndex2(pressure,temp,lambda);
		        n=refIndex2.nfind();
			lambda*=nbefore/n;	//wavelength changed in new medium			

//			cout<<"x="<<x<<", z="<<z<<", t="<<t<<"ns, theta="<<theta<<", nbefore="<<nbefore<<", n="<<n<<endl;
			
			if (cos(theta) > n/nbefore) 
			{
				cout << "***********************************" << endl;
				cout << "Total Internal Reflection occurs" << endl;
				cout << "***********************************" << endl;
				theta=-1.*theta;	 //Total Internal Reflection
			}
			else
			{
				// No total internal reflection
				// From snell's law, nb*cos(old_theta)=nf*cos(theta)
				theta=acos(nbefore*cos(old_theta)/n); 	//angle of refraction
				if (old_theta<0.) theta*=-1.;
			}

			double xnew=step*cos(theta);
			double znew=step*sin(theta);
			double tnew=n*step/lightSpeed;
			
		//
		// Absorption process turned on
		//
			if (absorption_turn_on){

				//photon absorbed
				if (Absorption(lambda,step,rand()%30000/30000.)) return;

				//if not absorbed, continue
			}

		//
		// Rayleigh Scattering turned on
		//
			if (rayleigh_turn_on){
				if (Rayleigh(lambda,step,rand()%30000/30000.)) {
					//ran is a random number determine at which point the particle is scattered
					//phi is the azimuthal angle runs from 0 to 2pi
					//psi is the scattered angle according to the distribution
					double ran=rand()%30000/30000.;
					double phi=rand()%30000/30000.*2.*PI;
					double psi=RayleighAngle();
				
					xnew=step*ran*cos(theta)+step*(1-ran)*cos(psi)*cos(theta) \
						+step*(1-ran)*sin(psi)*cos(phi)*sin(theta);
					znew=step*ran*sin(theta)+step*(1-ran)*cos(psi)*sin(theta) \
						-step*(1-ran)*sin(psi)*cos(phi)*cos(theta);
					tnew=n*step/lightSpeed;	
					theta=theta-atan(tan(psi)*cos(phi));

					//Save PROJECTED scattered angle
					ofstream scatterangle("scatterangle.dat",ios_base::app);
					double projected_angle=atan(tan(psi)*fabs(cos(phi)));
					if (projected_angle<0.0) projected_angle+=PI;
					scatterangle << projected_angle <<endl;
					scatterangle.close();
				}

			}


		// At last step
			if (fabs(x+xnew) >= Diameter/2. || fabs(z+znew) >= Height/2.)
			{
				// The photon will pass through the tank wall in this step
				// Use the WithoutDiffract() to confine the photon within the tank
				WithoutDiffract();

				// exit the while-loop
				break;
			}
			x+=xnew;
			z+=znew;
			t+=tnew; 
			nbefore=n;
		}

	} 
}


double reduced_angle(double theta){

//reduce theta to the range of [-pi,pi] in radian

	if (theta>=0.){
		double within2pi=theta-2*PI*int(theta/2./PI);
		if (within2pi < PI) return within2pi;
		if (within2pi >= PI) return within2pi-2.*PI;
	}

	if (theta<0.){
		double within2pi=theta-2*PI*int(theta/2./PI);
		if (within2pi < -1.*PI) return 2.*PI+within2pi;
		if (within2pi >= -1.*PI) return within2pi;
	}
}


void SingleRay::WithoutDiffract(){

//The trajectory made by this function is a straight line
//The photon's wavelength is constant, and it is stopped when it hits any of the walls
	
	//Calculate the current temperature and pressure
	double temp=CurrentTemp(z);
        double pressure=CurrentPressure(z);

	//Calculate the refractive index
        RefIndex refIndex(pressure,temp,lambda);
        double n=refIndex.nfind();

	cout << "WithoutDiffract():Refractive Index = "<< n << endl;

	// The photon hits the left vertical wall
	if (reduced_angle(theta)>atan2(Height/2.-z,-Diameter/2.-x)) {
		cout <<"The photon hits the left vertical wall"<< endl;
                z+=(-1.*Diameter/2.-x)*tan(theta);
                t+=n*fabs(-Diameter/2.-x)/lightSpeed/fabs(cos(theta));
                x=-1.*Diameter/2.;
	}

	// The photon hits the upper horizontal wall
	else if (reduced_angle(theta)>atan2(Height/2.-z,Diameter/2.-x)) {
                cout <<"The photon hits the upper horizontal wall"<<endl;
                x+=(Height/2.-z)/tan(theta);
                t+=n*fabs(Height/2.-z)/lightSpeed/fabs(sin(theta));
                z=Height/2.;
	}

	// The photon hits the right vertical wall 
	else if (reduced_angle(theta)>atan2(-Height/2.-z,Diameter/2.-x)) {
                cout <<"The photon hits the right vertical wall"<<endl;
                z+=(Diameter/2.-x)*tan(theta);
                t+=n*fabs(Diameter/2.-x)/lightSpeed/fabs(cos(theta));
                x=Diameter/2.;
	}

	// The photon hits the lower horizontal wall
	else if (reduced_angle(theta)>atan2(-Height/2.-z,-Diameter/2.-x)) {
                cout <<"The photon hits the lower horizontal wall"<<endl;
                x+=(-Height/2.-z)/tan(theta);
                t+=n*fabs(-Height/2.-z)/lightSpeed/fabs(sin(theta));
                z=-1.*Height/2.;		
	}

	// The photon hits the left vertical wall
	else {
                cout <<"The photon hits the left vertical wall"<< endl;
                z+=(-1.*Diameter/2.-x)*tan(theta);
                t+=n*fabs(-Diameter/2.-x)/lightSpeed/fabs(cos(theta));
                x=-1.*Diameter/2.;
	}

}
