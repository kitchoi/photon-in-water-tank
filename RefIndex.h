#include <math.h>
#include <iostream>
using namespace std;

// #include "WaterDensity.h"

class RefIndex{
 public:
  RefIndex(double p, double t, double lambda);
  double nfind();

  double p; // pressure
  double t; // temperature
  double lambda; // wavelength
  double rho_; //water density
};
