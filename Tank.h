#include <iostream>
#include <math.h>
using namespace std;

class Tank{
 public:
  Tank(double diameter,double height,double TopTemp,double BottomTemp, double TopPressure);
  //  void SpacetimeDiff(double result[],double InitialPos[],double InitialTheta,double InitialWLength,double step);
  //  void SingleRay(double pos[],double InitialPos[],double InitialTheta,double InitialWLength,double step,int DiffOnOff); // no diffraction:0, with diffraction:1
  double getDiameter();
  double getHeight();
  double getTopTemp();
  double getBottomTemp();
  double getTopPressure();

 private:
  double diamater;
  double height;
  double TopTemp;
  double BottomTemp;
  double TopPressure;
  

};
