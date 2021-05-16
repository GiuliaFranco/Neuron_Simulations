#ifndef __Integrator_H__
#define __Integrator_H__
#include "HH.h"
#include "Stoch.h"

class Solver{
public:
	double I0;
	double Period;
	const double dt = 0.01;  // Time steps (ms)
	const double  T = 100000.;    // Total time of observations
	const double V_0= -60;
	const int num=(T/dt);
	float * V =  new float[num]; // Membrane potential list
	float * n = new float[num]; // n potential list
	float * m = new float[num];  // m potential list
	float * h = new float[num]; // h potential list
	float * s1 = new float[num];  // Slow sodium activation
	float * s2 = new float[num];  // Slow potassium activation
	float * t = new float[num];
	HH * model;
	Stoch_HH * new_model;
	Solver(double I,double Periodicity);
	void Update_Chi(double dt,double V,int i);
	void Update_Psi(double dt,double V,int i);
	double * Calculator_aux(double y_1,double y_2,double y_3,double y_4,double y_5,double y_6,double t);
	void Euler();
	void Runge_Kutta();
	void Temporal_Integral(std :: string funz,double s);
};

#endif