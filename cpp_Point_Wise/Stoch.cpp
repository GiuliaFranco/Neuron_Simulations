#include "Stoch.h"
#include <cmath>
#include <string>
#include <iostream>
#include <random>
#include <thread> 
#include <future>

using namespace std;

Stoch_HH :: Stoch_HH(double I,double P,double Tot,int mod){
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



double * Stoch_HH :: model_call(double *y,double *s,double t){
		switch(mode){
			case 1: 
				return this->HHS_Meir_stoc(y,s,t);
			case 2:
				return this->HHS_OUstoc(y,s,t);
			case 3:
				return this->HHS_Fox_stoc(y,s,t);
			case 4:
				return this->HHS_Guller(y,s,t);
			default: 
				return this->HHS_Meir_stoc(y,s,t);

		};
};
		


double * Stoch_HH :: model_core(double t,double gNa,double gK,double * y,double s_1,double s_2){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double gl=gbarl;
		double INa=gNa*(V-ENa);
		double IK=gK*(V-EK);
		double Il=gl*(V-El);
		double	I_in=current->Get_current(t);
		double res[6]={((1/Cm)*(I_in-(INa+IK+Il))),Phi*(constants->A_n(V)*(1-n)-constants->B_n(V)*n),
		Phi*(constants->A_m(V)*(1-m)-constants->B_m(V)*m),Phi*(constants->A_h(V)*(1-h)-constants->B_h(V)*h),
		s_1,s_2};
		return res;
	};

		
double * Stoch_HH :: HHS_Meir_stoc(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double gNa=gbarNa*pow(m,3)*h*s_1;
		double gK=gbarK*pow(n,4);
		std::default_random_engine generator;
  		std::normal_distribution<double> distribution (0.0,1.5);

		return this->model_core(t,gNa,gK,y,
			(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1)+sqrt((constants->D_s1(V)*(1-s_1)+constants->G_s1(V)*s_1)/N)*distribution(generator),0);
};

float Stoch_HH :: Fox_stoc_gating(double n,double alpha,double beta,int Num){
  	std::default_random_engine generator1;
  	std::random_device rd;
  	generator1.seed( rd() );
  	double sigma = (alpha*(1-n)+beta*n)/Num;
  	std::normal_distribution<double> distribution1 (0.0,sigma);
  	float res=Phi*(alpha*(1-n)-beta*n+distribution1(generator1));
  	/*if(flag==1){
  		float res=Phi*(alpha/(alpha+beta)+distribution1(generator1));

  	}*/
  	//float res=Phi*(alpha*(1-n)-beta*n+sqrt(2*0.001*(alpha*(1-n)+beta*n)*((16.11809-log(distribution(generator)))*10000000)/Num)*cos(2*3.14*distribution(generator)));
  	return res;
}

float Stoch_HH :: M_stoc_gating(double n,double alpha,double beta,int Num){
  	std::default_random_engine generator1;
  	std::random_device rd;
  	generator1.seed( rd() );
  	double sigma = (alpha/(alpha+beta))/Num;
  	std::normal_distribution<double> distribution1 (0.0,sigma);
  	float res=Phi*(alpha/(alpha+beta)+distribution1(generator1));
  	//float res=Phi*(alpha*(1-n)-beta*n+sqrt(2*0.001*(alpha*(1-n)+beta*n)*((16.11809-log(distribution(generator)))*10000000)/Num)*cos(2*3.14*distribution(generator)));
  	return res;
}

float Stoch_HH :: Fox_stoc_gatingWB(double V){
  	std::default_random_engine generator1;
  	double alpha=constants->A_m(V);
  	double beta=constants->B_m(V);
  	double sigma = (alpha*beta/(alpha+beta))/N_NA;
  	std::normal_distribution<double> distribution1 (0.0,sigma);
  	return distribution1(generator1);
}

double * Stoch_HH :: HHS_Fox_stoc(double *y,double *s,double t){
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
  		std::normal_distribution<double> distribution (0.0,1);
		double res[6]={((1/Cm)*(I_in-(INa+IK+Il))),this->Fox_stoc_gating(n,constants->A_n(V),constants->B_n(V),4*N_K),
		this->Fox_stoc_gating(m,constants->A_m(V),constants->B_m(V),3*N_NA),this->Fox_stoc_gating(h,constants->A_h(V),constants->B_h(V),N_NA),
		this->Fox_stoc_gating(s_1,constants->D_s1(V),constants->G_s1(V),N)/Phi,0};
		//constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1,0};
		return res;
};

double * Stoch_HH :: HHS_OUstoc(double *y,double *s,double t){
	/*	double V=y[0];
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
			//(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1),0);*/
};

double * Stoch_HH :: HHS_Guller(double *y,double *s,double t){
		double V=y[0];
		double n=y[1];
		double m=y[2];
		double h=y[3];
		double s_1=s[0];
		double q_Na=s[3];
		double q_K=s[1];
		double p_Na=s[4];
		double p_K=s[2];
		double gammak=10;
		double gammaNa=10;
		double tau=10000000;
		double Wk=150;
		double WNa=200;
		double Tk=400;
		double TNa=800;
		double gNa=gbarNa*(pow(m,3)*h+sqrt((pow(m,3)*(1-(pow(m,3))))/N_NA)*h*q_Na)*s_1;
		double gK=gbarK*(pow(n,4)+sqrt((pow(n,4)*(1-(pow(n,4))))/N_K)*q_K);
		/*
		this->generator;
		*/
		/*
		if((t==3600.) || (t==7200.) || (t==10800.)){
			this->generator.seed( rd() );
		}
		*/
		std::random_device rd;
		std::default_random_engine generator;
		generator.seed( rd() );
		std::normal_distribution<double> distribution (0.0,1.5);
		std::normal_distribution<double> distributionK (0.0,gammak*Tk*(constants->A_n(V)*(1-n)+constants->B_n(V)*n));
		std::normal_distribution<double> distributionNa (0.0,gammaNa*TNa*(constants->A_m(V)*(1-m)+constants->B_m(V)*m));

  		double epsK=distributionK(generator);
  		double epsNa=distributionNa(generator);
  		//return this->model_core(t,gNa,gK,y,
		//		(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1),0);
			
  		double gl=gbarl;
		double INa=gNa*(V-ENa);
		double IK=gK*(V-EK);
		double Il=gl*(V-El);
		double I_in=current->Get_current(t);
		
		double res[9]={((1/Cm)*(I_in-(INa+IK+Il))),this->Fox_stoc_gating(n,constants->A_n(V),constants->B_n(V),4*N_K),
		this->M_stoc_gating(m,constants->A_m(V),constants->B_m(V),3*N_NA),this->Fox_stoc_gating(h,constants->A_h(V),constants->B_h(V),N_NA),
		(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1)+sqrt((constants->D_s1(V)*(1-s_1)+constants->G_s1(V)*s_1)/(N))*distribution(generator),
		p_K/tau,(-gammak*p_K-pow(Wk,2)*(constants->A_n(V)*(1-n)-constants->B_n(V)*n)*q_K+epsK)/tau,p_Na/tau,
		(-gammaNa*p_Na-pow(WNa,2)*(constants->A_m(V)*(1-n)-constants->B_m(V)*m)*q_Na+epsNa)/tau};
		
		/*double res[9]={((1/Cm)*(I_in-(INa+IK+Il))),Phi*(constants->A_n(V)*(1-n)-constants->B_n(V)*n),
		Phi*(constants->A_m(V)*(1-m)-constants->B_m(V)*m),Phi*(constants->A_h(V)*(1-h)-constants->B_h(V)*h),
		(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1)+sqrt((constants->D_s1(V)*(1-s_1)+constants->G_s1(V)*s_1)/(N))*distribution(generator),
		p_K/tau,(-gammak*p_K-pow(Wk,2)*(constants->A_n(V)*(1-n)-constants->B_n(V)*n)*q_K+epsK)/tau,p_Na/tau,
		(-gammaNa*p_Na-pow(WNa,2)*(constants->A_m(V)*(1-n)-constants->B_m(V)*m)*q_Na+epsNa)/tau};*/
		return res;
		//substitute q_k con p_k, q_Na con p_Na
		//TODO, add noise to m gate.
};


double * Stoch_HH :: HHS_Orio(double y_1,double y_2,double t,float  n,float  mh){
	double V=y_1;
	double s_1=y_2;
	double gNa=gbarNa*mh*s_1;
	double gK=gbarK*n;
	double gl=gbarl;
	double INa=gNa*(V-ENa);
	double IK=gK*(V-EK);
	double Il=gl*(V-El);
	double	I_in=current->Get_current(t);
	double res[2]={((1/Cm)*(I_in-(INa+IK+Il))),(constants->D_s1(V)*(1-s_1)-constants->G_s1(V)*s_1)};
	return res;
};

double * Stoch_HH :: Update_Orio_Psi(double * Psi_up,double V){
	double alpha=constants->A_n(V);
	double beta=constants->B_n(V);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
  	double eps1=distribution(generator);
  	double eps2=distribution(generator);
  	double eps3=distribution(generator);
  	double eps4=distribution(generator);
	double n0=(-4*alpha*Psi_up[0]+beta*Psi_up[1])+sqrt(4*alpha*Psi_up[0]+beta*Psi_up[1])*eps1/sqrt(N_K);
	double n1=(4*alpha*Psi_up[0]-beta*Psi_up[1])+(2*beta*Psi_up[2]-3*alpha*Psi_up[1])+
	sqrt(2*beta*Psi_up[2]+3*alpha*Psi_up[1])*eps2/sqrt(N_K)-sqrt(4*alpha*Psi_up[0]+beta*Psi_up[1])*eps1/sqrt(N_K);
	double n2=(3*alpha*Psi_up[1]-2*beta*Psi_up[2])+(3*beta*Psi_up[3]-2*alpha*Psi_up[2])-
	sqrt(3*alpha*Psi_up[1]+2*beta*Psi_up[2])*eps2/sqrt(N_K)+sqrt(3*beta*Psi_up[3]+2*alpha*Psi_up[2])*eps3/sqrt(N_K);
	double n3=(2*alpha*Psi_up[2]-3*beta*Psi_up[3])+(4*beta*Psi_up[4]-alpha*Psi_up[3])-
	sqrt(2*alpha*Psi_up[2]+3*beta*Psi_up[3])*eps3/sqrt(N_K)+sqrt(4*beta*Psi_up[4]+alpha*Psi_up[3])*eps4/sqrt(N_K);
	double n4=(alpha*Psi_up[3]-4*beta*Psi_up[4])-sqrt(alpha*Psi_up[3]+4*beta*Psi_up[4])*eps4/sqrt(N_K);
	double res[5]={n0,n1,n2,n3,n4};
	return res;
};

double * Stoch_HH :: Update_Orio_Chi(double * Chi_up,double V){
	double alpham=constants->A_m(V);
	double betam=constants->B_m(V);
	double alphah=constants->A_h(V);
	double betah=constants->B_h(V);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
  	double eps10=distribution(generator);
  	double eps20=distribution(generator);
  	double eps31=distribution(generator);
  	double eps42=distribution(generator);
	double eps23=distribution(generator);
  	double eps53=distribution(generator);
  	double eps64=distribution(generator);
	double eps54=distribution(generator);
	double eps35=distribution(generator);
	double eps75=distribution(generator);
	double eps76=distribution(generator);

  	double s0=(betam*Chi_up[2]-3*alpham*Chi_up[0])+(betah*Chi_up[1]-alphah*Chi_up[0])+
  	sqrt(betam*Chi_up[2]+3*alpham*Chi_up[0])*eps20/sqrt(N_NA)+sqrt(betah*Chi_up[1]+alphah*Chi_up[0])*eps10/sqrt(N_NA);

  	double s1=(alphah*Chi_up[0]-betah*Chi_up[1])+(betam*Chi_up[3]-3*alpham*Chi_up[1])-
  	sqrt(betah*Chi_up[1]+alphah*Chi_up[0])*eps10/sqrt(N_NA)+sqrt(betam*Chi_up[3]+3*alpham*Chi_up[1])*eps31/sqrt(N_NA);

  	double s2=(3*alpham*Chi_up[0]-betam*Chi_up[2])+(2*betam*Chi_up[4]-2*alpham*Chi_up[2])+(betah*Chi_up[3]-alphah*Chi_up[2])-
  	sqrt(3*alpham*Chi_up[0]+betam*Chi_up[2])*eps20/sqrt(N_NA)+sqrt(2*betam*Chi_up[4]+2*alpham*Chi_up[2])*eps42/sqrt(N_NA)+
  	sqrt(betah*Chi_up[3]+alphah*Chi_up[2])*eps23/sqrt(N_NA);


  	double s3=(3*alpham*Chi_up[1]-betam*Chi_up[3])+(2*betam*Chi_up[5]-2*alpham*Chi_up[3])-(alphah*Chi_up[2]-betah*Chi_up[3])-
  	sqrt(3*alpham*Chi_up[1]+betam*Chi_up[3])*eps31/sqrt(N_NA)+sqrt(2*betam*Chi_up[5]+2*alpham*Chi_up[3])*eps53/sqrt(N_NA)-
  	sqrt(alphah*Chi_up[2]+betah*Chi_up[3])*eps23/sqrt(N_NA);

  	double s4=(2*alpham*Chi_up[2]-2*betam*Chi_up[4])+(3*betam*Chi_up[6]-alpham*Chi_up[4])+(betah*Chi_up[5]-alphah*Chi_up[4])-
  	sqrt(2*alpham*Chi_up[2]+2*betam*Chi_up[4])*eps42/sqrt(N_NA)+sqrt(3*betam*Chi_up[6]+alpham*Chi_up[4])*eps64/sqrt(N_NA)+
  	sqrt(betah*Chi_up[5]+alphah*Chi_up[4])*eps54/sqrt(N_NA);

  	double s5=(2*alpham*Chi_up[3]-2*betam*Chi_up[5])+(3*betam*Chi_up[7]-alpham*Chi_up[5])+(betah*Chi_up[4]-alphah*Chi_up[5])-
  	sqrt(2*alpham*Chi_up[3]+2*betam*Chi_up[5])*eps35/sqrt(N_NA)+sqrt(3*betam*Chi_up[7]+alpham*Chi_up[5])*eps75/sqrt(N_NA)-
  	sqrt(betah*Chi_up[5]+alphah*Chi_up[4])*eps54/sqrt(N_NA);


  	double s6=(alpham*Chi_up[4]-3*betam*Chi_up[6])+(betah*Chi_up[7]-alphah*Chi_up[6])-
  	sqrt(alpham*Chi_up[4]+3*betam*Chi_up[6])*eps64/sqrt(N_NA)+sqrt(betah*Chi_up[7]+alphah*Chi_up[6])*eps76/sqrt(N_NA);

  	double s7=(alpham*Chi_up[5]-3*betam*Chi_up[7])-(-betah*Chi_up[7]+alphah*Chi_up[6])-
  	sqrt(alpham*Chi_up[5]+3*betam*Chi_up[7])*eps75/sqrt(N_NA)-sqrt(betah*Chi_up[7]+alphah*Chi_up[6])*eps76/sqrt(N_NA);

  	double res[8]={s0,s1,s2,s3,s4,s5,s6,s7};
	return res;

};

float Stoch_HH :: Update_Chi_0(double dt,float V,float olds){ 
	double tau=1/(constants->A_h(V)+constants->B_h(V));
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(pow(mbar,6)*hbar*(1-hbar)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float Stoch_HH :: Update_Chi_1(double dt,float V,float olds){ 
	double tau=1/(constants->A_m(V)+constants->B_m(V));
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(3*pow(mbar,5)*pow(hbar,2)*(1-mbar)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float Stoch_HH :: Update_Chi_2(double dt,float V,float olds){ 
	double tau=(1/(constants->A_m(V)+constants->B_m(V)))/2;
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(3*pow(mbar,4)*pow(hbar,2)*pow((1-mbar),2)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float Stoch_HH :: Update_Chi_3(double dt,float V,float olds){ 
	double tau=(1/(constants->A_m(V)+constants->B_m(V)))/3;
	double mbar=constants->A_m(V)/(constants->A_m(V)+constants->B_m(V));
	double hbar=constants->A_h(V)/(constants->A_h(V)+constants->B_h(V));
	double var=sqrt(pow(mbar,3)*pow(hbar,2)*pow((1-mbar),3)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float Stoch_HH :: Update_Chi_4(double dt,float V,float olds){ 
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

float Stoch_HH :: Update_Chi_5(double dt,float V,float olds){ 
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

float Stoch_HH :: Update_Chi_6(double dt,float V,float olds){ 
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

float Stoch_HH :: Update_Psi_0(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(4*pow(nbar,7)*pow((1-nbar),1)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);

};

float Stoch_HH :: Update_Psi_1(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun/2;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(6*pow(nbar,6)*pow((1-nbar),2)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};


float Stoch_HH :: Update_Psi_2(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun/3;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(4*pow(nbar,5)*pow((1-nbar),3)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

float Stoch_HH :: Update_Psi_3(double dt,float V,float olds){
	double taun=1/(constants->A_n(V)+constants->B_n(V));
	double tau= taun/4;
	double nbar=constants->A_n(V)/(constants->A_n(V)+constants->B_n(V));
	double var=sqrt(1*pow(nbar,4)*pow((1-nbar),4)/N);
	std::default_random_engine generator;
  	std::normal_distribution<double> distribution (0.0,1);
	return exp(-dt/tau)*olds+var*sqrt(1-exp(-2*dt/tau))*distribution(generator);
};

