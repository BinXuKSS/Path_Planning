#include "vehicle.h"
#include <iostream>

Vehicle::Vehicle() {};
Vehicle::Vehicle(double s, double s_d, double s_dd, double d, double d_d, double d_dd)
{
	this->s = s;
	this->s_d = s_d;
	this->s_dd = s_dd;
	this->d = d;
	this->d_d = d_d;
	this->d_dd = d_dd;
}

