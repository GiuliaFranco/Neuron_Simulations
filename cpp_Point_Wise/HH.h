#ifndef __HH_H__
#define __HH_H__
#include "Alpha_Beta.h"
#include "In_current.h"
#include <string>

class HH{

public:
	const double Cm = 0.5;  // Membrane Capacitance (uF/cm^2)
	const double Phi= 2.;  // Temperture factor
	const double V_0=-60;  // mV Initial value
	const double ENa=50; // mV Na potential
	const double  EK=-77; // mV K potential
	const double El=-54; // mV Leakage potential
	const double gbarNa=120; // (cm^2/kΩ) Na conductance evolving m
	//const double gbarNa=80; // (cm^2/kΩ) Na conductance not evolving m
	const double gbarK=36; //(cm^2/kΩ) K conductance
	const double gbarl=0.3; // (cm^2/kΩ) l conductance
	const double N= 1000000; //number of channels
	const int N_NA=240000;
	const int N_K=72000;
	const double gbarm=0.01*gbarK; 
	//const double gbarm=0.05*gbarK;
	int  mode;  // model selection
	int eps;   // evolution of slow dynamics  (0 -> not evolving,1 -> evolving)
	double s1_const;
	double s1_init;
	double s2_const; 
	double I_0;
	double Periodicity;
	double T;
	double Na_noise;
	double K_noise;
	Alpha_Beta * constants=new Alpha_Beta();
	In_current * current;
	HH(double I,double P,double Tot,int mod);
	double * model_call(double *y,double *s,double t);
	double * model_core(double t,double gNa,double gK,double * y,double s_1,double s_2);
	double * HHSAP(double *y,double *s,double t);
	double * HHS(double *y,double *s,double t);
	double * HHSIP(double *y,double *s,double t);
	double * HH_basic(double *y,double t);
	/*double * HHS_Meir_stoc(double *y,double *s,double t);
	double * HHS_OUstoc(double *y,double *s,double t);
	float Fox_stoc_gating(double n,double alpha,double beta,int Num);
	float Fox_stoc_gatingWB(double V);
	double * HHS_Fox_stoc(double *y,double *s,double t);*/
	double updateV(double * y,double t);
	double updaten(double V,double n);
	double updatem(double V,double m);
	double updateh(double V,double h);
	/*float  Update_Chi_0(double dt,float V,float olds);
	float  Update_Chi_1(double dt,float V,float olds);
	float  Update_Chi_2(double dt,float V,float olds);
	float  Update_Chi_3(double dt,float V,float olds);
	float  Update_Chi_4(double dt,float V,float olds);
	float  Update_Chi_5(double dt,float V,float olds);
	float  Update_Chi_6(double dt,float V,float olds);
	float Update_Psi_0(double dt,float V,float olds);
	float Update_Psi_1(double dt,float V,float olds);
	float Update_Psi_2(double dt,float V,float olds);
	float Update_Psi_3(double dt,float V,float olds);*/
};
#endif