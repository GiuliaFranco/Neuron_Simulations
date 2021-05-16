#include "HH.h"
#include <cmath>
#include <string>
#include <iostream>
#include <random>
#include <thread> 
#include <future>

using namespace std;

HH :: HH(double I,double P,double Tot,int mod){
	I_0=I;
	Periodicity=P;
	T=Tot;
	current=new In_current(I_0,T,Periodicity);
	eps=1;   // evolution of slow dynamics  (0 -> not evolving,1 -> evolving)
	s1_const=1;
	s2_const=1;
	s1_init=0.;
	Na_noise=0;
	K_noise=0;
	mode=mod;
};



double * HH :: model_call(double *y,double *s,double t){
		switch(mode){
			case 1: 
				return this->HHS(y,s,t);
			case 2: 
				return this->HHSAP(y,s,t);
			case 3: 
				return this->HHSIP(y,s,t);
			/*case 4: 
				return this->HHS_Meir_stoc(y,s,t);
			case 5:
				return this->HHS_OUstoc(y,s,t);
			case 6:
				return this->HHS_Fox_stoc(y,s,t);*/
			default: 
				return this->HH_basic(y,t);

		};
};
		


double * HH :: model_core(double t,double gNa,double gK,double * y,double s_1,double s_2){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double gl=gbarl;
		double INa=gNa*(V-ENa);
		double IK=gK*(V-EK);
		double Il=gl*(V-El);
		double	I_in=current->Get_current(t);
		//cout<<constants->A_m(V)<<" "<<constants->B_m(V)<<endl;
		double res[6]={((1/Cm)*(I_in-(INa+IK+Il))),Phi*(constants->A_n(V)*(1-n)-constants->B_n(V)*n),
		Phi*(constants->A_m(V)*(1-m)-constants->B_m(V)*m),Phi*(constants->A_h(V)*(1-h)-constants->B_h(V)*h),
		s_1,s_2};
		return res;
	};


double HH :: updateV(double * y,double t){
	double V=y[0];
	double n=y[1];
	double m=y[2];
	double h=y[3];
	double gNa=gbarNa*pow(m,3)*h;
	double gK=gbarK*pow(n,4);
	double gl=gbarl;
	double INa=gNa*(V-ENa);
	double IK=gK*(V-EK);
	double Il=gl*(V-El);
	double I_in=current->Get_current(t);
	return ((1/Cm)*(I_in-(INa+IK+Il)));
		
};
	
double HH :: updaten(double V,double n){ return Phi*(constants->A_n(V)*(1-n)-constants->B_n(V)*n);};

double HH :: updatem(double V,double m){ return Phi*(constants->A_m(V)*(1-m)-constants->B_m(V)*m);};

double HH :: updateh(double V,double h){ return Phi*(constants->A_h(V)*(1-h)-constants->B_h(V)*h);};

double * HH :: HH_basic(double * y,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double gNa=gbarNa*pow(m,3)*h;
		double gK=gbarK*pow(n,4);
		return this->model_core(t,gNa,gK,y,0,0);
/*	std::future<double> th1 = std::async(std::launch::async,&HH::updateV,this,y,t);
	double V =  th1.get();
	std::future<double> th2 = std::async(std::launch::async,&HH::updaten,this,y[0],y[1]);
	double n = th2.get();
	std::future<double> th3 = std::async(std::launch::async,&HH::updatem,this,y[0],y[2]);
	double m = th3.get();
	std::future<double> th4 = std::async(std::launch::async,&HH::updateh,this,y[0],y[3]);
	double h = th4.get();
	double res[6]={V,n,m,h,0,0};
	return res;*/
		
};
		

double * HH :: HHS(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double gNa=gbarNa*pow(m,3)*h*s_1;
		//double gNa=gbarNa*pow(m,3)*h*constants->SimplifiedCBM_s1(t,s1_init);
		double gK=gbarK*pow(n,4);

		return this->model_core(t,gNa,gK,y,(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1),0);
	
};
			
double * HH :: HHSAP(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double s_2=s[1];
		double gNa=gbarNa*pow(m,3)*h*s_1;
		double gK=(gbarK+gbarm*s_2)*pow(n,4);
		return this->model_core(t,gNa,gK,y,constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1,
			constants->D_s2(V)*(1-s_2)-constants->G_s2(V)*s_2);

};


double * HH :: HHSIP(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double s_2=s[1];
		double gNa=gbarNa*pow(m,3)*h*s_1;
		double gK=(gbarK+gbarm*s_2)*pow(n,4);
		return this->model_core(t,gNa,gK,y,constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1,
			constants->G_s2(V)*(1-s_2)-constants->D_s2(V)*s_2);
};
	/*	
double * HH :: HHS_Meir_stoc(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double gNa=gbarNa*pow(m,3)*h*s_1;
		//double gNa=gbarNa*pow(m,3)*h*constants->SimplifiedCBM_s1(t,s1_init);
		double gK=gbarK*pow(n,4);
		std::default_random_engine generator;
  		std::normal_distribution<double> distribution (0.0,1.5);

		return this->model_core(t,gNa,gK,y,
			(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1)+sqrt((constants->D_s1(V)*(1-s_1)+constants->G_s1(V)*s_1)/N)*distribution(generator),0);
};

float HH :: Fox_stoc_gating(double n,double alpha,double beta,int Num){
  	std::default_random_engine generator1;
  	double sigma = (alpha*(1-n)+beta*n)/Num;
  	std::normal_distribution<double> distribution1 (0.0,sigma);
  	float res=Phi*(alpha*(1-n)-beta*n+distribution1(generator1));
  	//float res=Phi*(alpha*(1-n)-beta*n+sqrt(2*0.001*(alpha*(1-n)+beta*n)*((16.11809-log(distribution(generator)))*10000000)/Num)*cos(2*3.14*distribution(generator)));
  	return res;
}

float HH :: Fox_stoc_gatingWB(double V){
  	std::default_random_engine generator1;
  	double alpha=constants->A_m(V);
  	double beta=constants->B_m(V);
  	double sigma = (alpha*beta/(alpha+beta))/N_NA;
  	std::normal_distribution<double> distribution1 (0.0,sigma);
  	return distribution1(generator1);
}

double * HH :: HHS_Fox_stoc(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double gl=gbarl;
		double gNa=gbarNa*pow(m,3)*h*s_1;
		double gK=gbarK*pow(n,4);
		double INa=gNa*(V-ENa);
		double IK=gK*(V-EK);
		double Il=gl*(V-El);
		double	I_in=current->Get_current(t);
		std::default_random_engine generator;
  		std::normal_distribution<double> distribution (0.0,1.5);
		double res[6]={((1/Cm)*(I_in-(INa+IK+Il))),this->Fox_stoc_gating(n,constants->A_n(V),constants->B_n(V),4*N_K),
		this->Fox_stoc_gating(m,constants->A_m(V),constants->B_m(V),3*N_NA),this->Fox_stoc_gating(h,constants->A_h(V),constants->B_h(V),N_NA),
		this->Fox_stoc_gating(s_1,constants->D_s1(V),constants->G_s1(V),N_NA+N_K),0};
		return res;
};

double * HH :: HHS_OUstoc(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double gNa=gbarNa*(pow(m,3)*h+this->Na_noise)*s_1;
		double gK=gbarK*(pow(n,4)+this->K_noise);
		std::default_random_engine generator;
  		std::normal_distribution<double> distribution (0.0,1.5);
		return this->model_core(t,gNa,gK,y,
				(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1)+sqrt((constants->D_s1(V)*(1-s_1)+constants->G_s1(V)*s_1)/N)*distribution(generator),0);
			//(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1),0);

};

float HH :: Update_Chi_0(double dt,float V,float olds){ 
	double tau=1/(constants->A_h(V)+constants->B_h(V));
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(pow(mbar,6)*hbar*(1-hbar)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Chi_1(double dt,float V,float olds){ 
	double tau=1/(constants->A_m(V)+constants->B_m(V));
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(3*pow(mbar,5)*pow(hbar,2)*(1-mbar)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Chi_2(double dt,float V,float olds){ 
	double tau=(1/(constants->A_m(V)+constants->B_m(V)))/2;
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(3*pow(mbar,4)*pow(hbar,2)*pow((1-mbar),2)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Chi_3(double dt,float V,float olds){ 
	double tau=(1/(constants->A_m(V)+constants->B_m(V)))/3;
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(pow(mbar,3)*pow(hbar,2)*pow((1-mbar),3)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Chi_4(double dt,float V,float olds){ 
	double taum=1/(constants->A_m(V)+constants->B_m(V));
	double tauh=1/(constants->A_h(V)+constants->B_h(V));
	double tau= tauh*taum/(taum+tauh);
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(3*pow(mbar,5)*pow(hbar,1)*pow((1-mbar),1)*pow((1-hbar),1)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Chi_5(double dt,float V,float olds){ 
	double taum=1/(constants->A_m(V)+constants->B_m(V));
	double tauh=1/(constants->A_h(V)+constants->B_h(V));
	double tau= tauh*taum/(taum+2*tauh);
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(3*pow(mbar,4)*pow(hbar,1)*pow((1-mbar),2)*pow((1-hbar),1)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Chi_6(double dt,float V,float olds){ 
	double taum=1/(constants->A_m(V)+constants->B_m(V));
	double tauh=1/(constants->A_h(V)+constants->B_h(V));
	double tau= tauh*taum/(taum+3*tauh);
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(1*pow(mbar,3)*pow(hbar,1)*pow((1-mbar),3)*pow((1-hbar),1)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Psi_0(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(4*pow(nbar,7)*pow((1-nbar),1)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);

};

float HH :: Update_Psi_1(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun/2;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(6*pow(nbar,6)*pow((1-nbar),2)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};


float HH :: Update_Psi_2(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun/3;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(4*pow(nbar,5)*pow((1-nbar),3)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float HH :: Update_Psi_3(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun/4;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(1*pow(nbar,4)*pow((1-nbar),4)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};


*/

