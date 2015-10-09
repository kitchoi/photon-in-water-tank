#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

double Interpolate(string filename_str, double x0){
        const char* filename=filename_str.c_str();
        ifstream datafile(filename,ifstream::in);
        if (!datafile.is_open()) {
                cerr << ":::::::Error opening "<<filename<<" for interpolation"<<endl;
                exit(1);                //terminate the main program
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
        cerr<<":::::::Out of range error! cannot interpolate with file "<<filename<<" and x= "<<x0<<endl;
        exit(1);

        //Extrapolation
        //return y+slope*(x0-x);
}

int main(){

	double log_a;
	for (int i=1;i<=1000;i++){
		double wavelength=210.+400*i/1000.;
		log_a=Interpolate("fry.dat",wavelength);
		cout << wavelength << '\t' << pow(10,log_a) <<endl;
	}
}
