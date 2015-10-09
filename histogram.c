#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

int main(){
	const char* filename="hisInput.dat";	

	ifstream input_dump(filename,ifstream::in);
	double xmax=0.;
	double xmin=0.;
	double mean=0.;
	int counter=0;			// number of data
	int bin;			// number of bins
	double x;
	double x2sum=0;			// sum of (x-mean)*(x-mean)
	double width;			// bin width
	double scale=1.;		// scale up raw data input
	//
        // determine the maximum, minimum and the number of x values
	//
	input_dump >> x;
	x*=scale;
	xmax=x;
	xmin=x;
	mean+=x;
	counter++;

	while (true){
		input_dump >> x;
		x*=scale;
		if (input_dump.eof()) break;
		if (x>xmax) xmax=x;
		if (x<xmin) xmin=x;
		mean+=x;
		counter++;
	}

	bin = int(sqrt(counter));
	mean=mean/counter;			// calculate the mean of x
        width=(xmax-xmin)/bin;			// calculate bin width
        int *value = new int[bin];
	for (int j=0; j<bin; j++){
		value[j]=0;			// initialize histogram
	}	

	input_dump.close();


	ifstream input(filename,ifstream::in);
	while (true){
		input >> x;
		x*=scale;
		if (input.eof()) break;
		int i=int((x-xmin)/width);

		value[i]+=1;

		x2sum+=(x-mean)*(x-mean);
	}
	
	input.close();
	
	ofstream outfile ("hisOutput.dat",ofstream::out);
	
	//
	// Output histogram as probability density
	//
	outfile << "# Number of data =" << counter << endl;
	outfile << "# xmax= " << xmax/scale << ", xmin= " << xmin/scale << endl;
	outfile << "# Variance = " << x2sum/counter/scale/scale << endl;
	outfile << "# Standard Deviation =  "<< sqrt(x2sum/counter)/scale<< endl;
	outfile << "# Mean = " << mean/scale << endl;
	outfile << "# SD/Mean ="<< sqrt(x2sum/counter)/mean <<endl;

	for (int j=0;j<bin;j++){
		outfile << xmin+width*(j+0.5)<<'\t'<<value[j]/width/counter<<endl;
	}
	outfile.close();

        cout << "Number of data =" << counter << endl;
        cout << "xmax= " << xmax/scale << ", xmin= " << xmin/scale << endl;
        cout << "Variance = " << x2sum/counter/scale/scale << endl;
        cout << "Standard Deviation =  "<< sqrt(x2sum/counter)/scale << endl;
        cout << "Mean = " << mean/scale << endl;
	cout << "SD/Mean ="<< sqrt(x2sum/counter)/mean <<endl;

}
