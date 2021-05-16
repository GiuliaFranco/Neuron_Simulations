#include "Stoch_Integrator.h"
#include <iostream>
#include <thread>  
#include <future>
using namespace std;

Solver_Stoch :: Solver_Stoch(double I,double Periodicity){
	I0=I;
	Period=Periodicity;
	model=new Stoch_HH(I0,Period,T,3);
	t[0]=0;
	double start = 0;
	for(int i=0;i<=num;i++){ 
		start+=dt;
		t[i]=start; 
	};
	
	V[0]=V_0;
	m[0]=model->constants->A_m(V[0])/(model->constants->A_m(V[0])+model->constants->B_m(V[0]));
	n[0]=model->constants->A_n(V[0])/(model->constants->A_n(V[0])+model->constants->B_n(V[0]));
	h[0]=model->constants->A_h(V[0])/(model->constants->A_h(V[0])+model->constants->B_h(V[0]));
	q_Na[0] = 0;  
	q_K [0] = 0; 
	p_Na[0] = 0;  
	p_K [0] = 0; 
	if(model->eps!=0 && model->mode!=3){
	//	s2[0]=model->constants->D_s2(V[0])/(model->constants->D_s2(V[0])+model->constants->G_s2(V[0]));
		s1[0]=model->constants->D_s1(V[0])/(model->constants->D_s1(V[0])+model->constants->G_s1(V[0]));

		model->s1_init=s1[0];
	}
	if(model->eps!=0 && model->mode==3){
	//	s2[0]=model->constants->G_s2(V[0])/(model->constants->G_s2(V[0])+model->constants->D_s2(V[0]));
		s1[0]=model->constants->D_s1(V[0])/(model->constants->D_s1(V[0])+model->constants->G_s1(V[0]));
	}
	if(model->eps==0){
		s1[0]=model->s1_const;
	//	s2[0]=model->s2_const;
	}

};

void Solver_Stoch :: InitTime(double start){
	for(int i=0;i<=num;i++){ 
		start+=dt;
		t[i]=start; 
	};
	//model->current->T_tot+=start;
	model->current->start_i=start;
};

double * Solver_Stoch :: Calculator_aux(double y_1,double y_2,double y_3,double y_4,double y_5,double y_6,double t){
	double y[4] = {y_1,y_2,y_3,y_4};
	double s[2] = {y_5,y_6};
	return model->model_call(y,s,t);

};

double * Solver_Stoch :: Calculator_Guler(double y_1,double y_2,double y_3,double y_4,double y_5,double y_6,double y_7,double y_8,double y_9,double t){
	double y[4] = {y_1,y_2,y_3,y_4};
	double s[5] = {y_5,y_6,y_7,y_8,y_9};
	return model->model_call(y,s,t);
};

void Solver_Stoch :: Euler(){
	//cout<<V[0]<<" "<<t[0]<<endl;
	for(int i=1;i<=num;i++){
		//this->Update_Chi(dt,this->V[i-1],i);
		//this->Update_Psi(dt,this->V[i-1],i);
		double * K_1=Calculator_aux(V[i-1],n[0],m[0],h[0],s1[0],s2[0],t[i-1]);
		V[i]=(V[i-1]+K_1[0]*dt);
		n[0]=(n[0] +K_1[1]*dt);
		m[0]=(m[0] + K_1[2]*dt);
		//m[i]=model->constants->A_m(V[i])/(model->constants->A_m(V[i])+model->constants->B_m(V[i]));
		h[0]=(h[0]+K_1[3]*dt);
		if(model->eps!=0){
			//s1[i]=model->constants->SimplifiedCBM_s1(t[i],model->s1_init);
			s1[0]=(s1[0]+K_1[4]*dt);
			//s2[i]=(s2[i-1]+K_1[5]*dt);
		}
		else{
			s1[i]=model->s1_const;
			s2[i]=model->s2_const;
		}	
		//cout<<V[i]<<" "<<t[i]<<endl;
	};


};


void Solver_Stoch :: Euler_Guller(){
	//cout<<V[0]<<" "<<t[0]<<endl;
	for(int i=1;i<=num;i++){
		double * K_1=Calculator_Guler(V[i-1],n[0],m[0],h[0],s1[0],q_K[0],p_K[0],q_Na[0],p_Na[0],t[i-1]);
		V[i]=(V[i-1]+K_1[0]*dt);
		n[0]=(n[0] +K_1[1]*dt);
		//m[0]=(m[0] + K_1[2]*dt);
		//m[0]=model->M_stoc_gating(m[0],model->constants->A_m(V[i]),model->constants->B_m(V[i]),3*model->N_NA);
		m[0]=model->constants->A_m(V[i])/(model->constants->A_m(V[i])+model->constants->B_m(V[i]));
		//m[0]=K_1[2];
		h[0]=(h[0]+K_1[3]*dt);
		s1[0]=(s1[0]+K_1[4]*dt);
		q_K[0]=(q_K[0]+K_1[5]*dt);
		p_K[0]=(p_K[0]+K_1[6]*dt);
		q_Na[0]=(q_Na[0]+K_1[7]*dt);
		p_Na[0]=(p_Na[0]+K_1[8]*dt);
		//cout<<V[i]<<" "<<t[i]<<endl;
		/*if(i==num){
			cout<<V[i]<<endl;
			cout<<n[0]<<endl;
			cout<<m[0]<<endl;
			cout<<m[0]<<endl;
			cout<<s1[0]<<endl;
			cout<<q_Na[0]<<endl;
			cout<<q_K[0]<<endl;
			cout<<p_Na[0]<<endl;
			cout<<p_K[0]<<endl;
			cout<<t[i]<<endl;
		}*/
		//cout<<t[i]<<endl;
	};


};




