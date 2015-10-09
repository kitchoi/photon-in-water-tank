//Upon passing parameter to the calculator "RefIndex", units are as followed
//Unit of pressure: Pa
//Unit of temperature: Kelvin
//Unit of wavelength: nanometer

//This program requires file "GenerateIndex.txt" for iteration limits, the units in the file as followed
//Unit of pressure: atm <--- NOTE!!!
//Unit of temperature: Celcius <--- NOTE!!!
//Unit of wavelength: nanometer


#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "RefIndex.h"

#define atm 101325.	//1 atm = 101325 Pa

int main(){
double rindex;
double t;
double p;
double lambda;
double InputPara[9];
	for(int N=0;N<9;N++) InputPara[N]=-100.;
string ParaName[9]={"MinTemp","MaxTemp","TempStep","MinPressure","MaxPressure","PressureStep","MinWavelength","MaxWavelength","WavelengthStep"};
string Unit[3]={"C","atm","nm"};

//
// Read-in iteration limits
//
ifstream parameter("GenerateIndex.txt", ios::in);	// File contains parameter for run

if(!parameter.is_open()) {
	cout << "Error opening file GenerateIndex.txt" <<endl;
	return 1;
}

while(parameter.is_open()){
	string str1;
	string str2;
	double value;
	parameter >> str1 >> str2 >> value;	//Read in parameter name, unit and value
	
	if(parameter.eof()) break;
	
	for (int k=0; k<9; k++){
		if(str1.compare(ParaName[k])==0){
			if(InputPara[k]!=-100.) {
				cout << "Error! Duplicate input of parameter " << ParaName[k] << endl;
				return 1;
			}
			if(str2.compare(Unit[k/3])!=0) {
				cout << "Unit Error! Please use unit "<< Unit[k/3] <<" for "<< ParaName[k] <<endl;
				return 1;
			}
			InputPara[k]=value;
		}
	}
}
parameter.close();
//
// End read-in iteration limits
//

//Check if all parameters are initialized
{
	bool initialize=true;
	for (int k=0; k<9; k++){
		cout << ParaName[k] <<" is " <<InputPara[k] << Unit[k/3] <<endl;
		if (InputPara[k]==-100.){
			initialize=false;
			cout <<"Error!" << ParaName[k] << " is not initialized" <<endl;
		}
	}
	if (!initialize) return 1;
}



ofstream output("n.dat", ios::out);

cout << "Pressure(atm)" << '\t' << "Temperature(C)" << '\t'<< "Water Density(kg/m3)"<<'\t'<< "Wavelength(nm)" << '\t' << "RefIndex" << endl;
output << "#Pressure(atm)" << '\t' << "Temperature(C)" <<'\t'<<"Water Density(kg/m3)" <<'\t'<<"Wavelength(nm)" << '\t' << "RefIndex" << endl;

//Iterate Temperature
for (t=InputPara[0];t<=InputPara[1];t=t+InputPara[2]){
	
	//Iterate Pressure
	for (p=InputPara[3];p<=InputPara[4];p=p+InputPara[5]){

		//Iterate Wavelength
		for (lambda=InputPara[6];lambda<=InputPara[7];lambda=lambda+InputPara[8]){

			RefIndex refn(p*atm,t+273.15,lambda);
			rindex=refn.nfind();

			cout << p << '\t' << t << '\t' << (refn.rho_)*1000.<< '\t'<< lambda << '\t' << rindex << endl;
			output << p << '\t' << t << '\t' << (refn.rho_)*1000.<< '\t'<<lambda << '\t' << rindex << endl;

		} //end lambda iteration

	} // end pressure iteration

} //end temperature iteration

output.close();
} //end main program
