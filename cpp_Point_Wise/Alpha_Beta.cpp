#include "Alpha_Beta.h"
#include <cmath>
#include <iostream>
using namespace std;

Alpha_Beta :: Alpha_Beta(){};

float Alpha_Beta :: vtrap( float x,float y){
	float res;
	if (fabs(x/y) < 1e-6){
		res = y*(1 - x/y/2);
	}
	else{
		res = x/(exp(x/y) - 1);
	}
	return res;
};

//float Alpha_Beta ::  G_s1(float V){return 0.51/(exp(-0.3*(V+17.))+1);};	//0.00051
float Alpha_Beta :: G_s1(float V){return 0.00051/(exp(-0.3*(V+17.))+1);}; // high 0.001
//float Alpha_Beta :: D_s1(float V){return 0.05*exp(-(V+85.)/30.);};	//0.00005
float Alpha_Beta :: D_s1(float V){return 0.00005*exp(-(V+85.)/30.);}; // high 0.0001

float Alpha_Beta ::  G_s2(float V){
	float	sigma=1./28.;  //mV
	return (3.3*exp((V+35)*sigma) + exp(-(V+35)*0.05))/(1+ exp((V+35)*0.1));
};

float Alpha_Beta :: D_s2(float V){
	//HZ
	double	sigma=1./28.;  //mV
	return (3.3*exp((V+35)*sigma) + exp(-(V+35)*0.05))/(1+ exp(-(V+35)*0.1));
};

/*
float Alpha_Beta :: A_n( float V){ return .01*vtrap(-(V+34.),10.); };  //kHz

float Alpha_Beta :: A_m(float V){ return .1 * vtrap(-(V+35.),10.);  };  //kHz

float Alpha_Beta :: A_h(float V){ return .07 * exp(-(V+58.)/20.); };   //kHz

float Alpha_Beta :: B_n(float V){ return .125*exp(-(V+44.)/80.);  }; 
	
float Alpha_Beta :: B_m(float V){ return 4. * exp(-(V+60.)/18.);  };

float Alpha_Beta :: B_h(float V){ return 1. / (exp(-(V+28.)/10.) + 1);  };
*/


//WB
float Alpha_Beta :: A_n( float V){ return .01*vtrap(-(V+55.),10.); };  //kHz
float Alpha_Beta :: A_m(float V){ return .1 * vtrap(-(V+40.),10.);  };  //kHz


//HH
float Alpha_Beta :: A_n_HH( float V){ return .01*(V+55.)/(1-exp(-(V+55.)/10.)); };  //kHz
float Alpha_Beta :: A_m_HH(float V){ return .1 * (V+40.)/(1-exp(-(V+40.)/10.));  };  //kHz


float Alpha_Beta :: A_h(float V){ return .07 * exp(-(V+65.)/20.); };   //kHz

float Alpha_Beta :: B_n(float V){ return .125*exp(-(V+65.)/80.);  }; 
	
float Alpha_Beta :: B_m(float V){ return 4. * exp(-(V+65.)/18.);  };

float Alpha_Beta :: B_h(float V){ return 1. / (exp(-(V+35.)/10.) + 1);  };

float Alpha_Beta :: gamma_integral(double Periodicity,float * V){
	double gamma_aux= (V[(int)Periodicity*100+5*100]-V[(int)Periodicity*100-1*100])+(log(1+exp(-0.3*(V[(int)Periodicity*100+5*100]+17)))-log(1+exp(-0.3*(V[(int)Periodicity*100-1*100]+17))))/0.3;
	double gamma_aux1= (V[(int)Periodicity*100+(int)Periodicity*100]-V[(int)Periodicity*100+5*100])+(log(1+exp(-0.3*(V[(int)Periodicity*100+(int)Periodicity*100]+17)))-log(1+exp(-0.3*(V[(int)Periodicity*100+5*100]+17))))/0.3;
	double gamma_H=0.51*gamma_aux/0.005;
	double gamma_L=0.51*gamma_aux1/(Periodicity/1000-0.005);
	return (gamma_H-gamma_L)*0.005*1000/Periodicity+gamma_L;
};

float Alpha_Beta :: delta_integral(double Periodicity,float * V){
	double delta_aux= exp(-V[(int)Periodicity*100+5*100]/30.)-exp(-V[(int)Periodicity*100-1*100]/30.);
	double delta_aux1= exp(-V[(int)Periodicity*100+(int)Periodicity*100]/30.)-exp(-V[(int)Periodicity*100+5*100]/30.);
	double delta_H=-0.05*30.*exp(-85./30.)*delta_aux/0.01;
	double delta_L=-0.05*30.*exp(-85./30.)*delta_aux1/(Periodicity/1000-0.01);
	return (delta_H-delta_L)*0.01*1000/Periodicity+delta_L;
};

float Alpha_Beta :: SimplifiedCBM_s1(double t,double init){
	double gamma=0.00225;  //0.0025
	double delta=0.027;   // 0.027
	double s_inf=delta/(delta+gamma);
	double tau_inf=1/(delta+gamma);
	return s_inf+(init-s_inf)*exp(-t/(tau_inf*100));
}
