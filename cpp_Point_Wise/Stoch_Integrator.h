#ifndef __Stoch_Integrator_H__
#define __Stoch_Integrator_H__
#include "HH.h"
#include "Stoch.h"

class Solver_Stoch{
public:
	double I0;
	double Period;
	const double dt = 0.01;  //0.01 Time steps (ms)
	const double  T = 600000.;    // Total time of observations
	const double V_0= -60;
	const int num=(T/dt);
	float * V =  new float[num]; // Membrane potential list
	float * n = new float[1]; // n potential list
	float * m = new float[1];  // m potential list
	float * h = new float[1]; // h potential list
	float * s1 = new float[1];  // Slow sodium activation
	float * s2 = new float[1];  // Slow potassium activation
	float * q_Na = new float[1];  // Slow sodium activation
	float * q_K = new float[1];  // Slow potassium activation
	float * p_Na = new float[1];  // Slow sodium activation
	float * p_K = new float[1];  // Slow potassium activation
	float * t = new float[num];
	/*float** Chi=new float*[num+1];
	float** Psi=new float*[num+1];
	float** Chi_Orio=new float*[num+1];
	float** Psi_Orio=new float*[num+1];*/
	Stoch_HH * model;
	Solver_Stoch(double I,double Periodicity);
	void InitTime(double start);
	double * Calculator_aux(double y_1,double y_2,double y_3,double y_4,double y_5,double y_6,double t);
	double * Calculator_Guler(double y_1,double y_2,double y_3,double y_4,double y_5,double y_6,double y_7,double y_8,double y_9,double t);
	void Euler();
	void Euler_Orio();
	void Euler_Guller();
};

#endif
