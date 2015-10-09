#include "Tank.h"
using namespace std;

Tank::Tank(double inputden,double inputh,double inputTopT, double inputBottomT, double inputTopP){
  this->diameter=inputden;
  this->height=inputh;
  this->TopTemp=inputTopT;
  this->BottomTemp=inputBottomT;
  this->TopPressure=inputTopP;
}

double Tank::getDiameter(){
  return diameter;
}

double Tank::getHeight(){
  return height;
}

double Tank::getTopTemp(){
  return TopTemp;
}

double Tank::getBottomTemp(){
  return BottomTemp;
}

double Tank::getTopPressure(){
  return TopPressure;
}
