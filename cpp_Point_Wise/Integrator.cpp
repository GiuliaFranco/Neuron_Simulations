#include "Integrator.h"
#include <iostream>
#include <thread>  
#include <future>
using namespace std;

Solver :: Solver(double I,double Periodicity){
	I0=I;
	Period=Periodicity;
	model=new HH(I0,Period,T,4);
	//new_model=new Stoch_HH(I0,Period,T,1);
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
	if(model->eps!=0 && model->mode!=3){
		s2[0]=model->constants->D_s2(V[0])/(model->constants->D_s2(V[0])+model->constants->G_s2(V[0]));
		s1[0]=model->constants->D_s1(V[0])/(model->constants->D_s1(V[0])+model->constants->G_s1(V[0]));

		model->s1_init=s1[0];
	}
	if(model->eps!=0 && model->mode==3){
		s2[0]=model->constants->G_s2(V[0])/(model->constants->G_s2(V[0])+model->constants->D_s2(V[0]));
		s1[0]=model->constants->D_s1(V[0])/(model->constants->D_s1(V[0])+model->constants->G_s1(V[0]));
	}
	if(model->eps==0){
		s1[0]=model->s1_const;
		s2[0]=model->s2_const;
	}

};

/*
void Solver :: Update_Chi(double dt,double V,int i){
	this->Chi[i][0]=model->Update_Chi_0(dt,V,this->Chi[i-1][0]);
	this->Chi[i][1]=model->Update_Chi_1(dt,V,this->Chi[i-1][1]);
	this->Chi[i][2]=model->Update_Chi_2(dt,V,this->Chi[i-1][2]);
	this->Chi[i][3]=model->Update_Chi_3(dt,V,this->Chi[i-1][3]);
	this->Chi[i][4]=model->Update_Chi_4(dt,V,this->Chi[i-1][4]);
	this->Chi[i][5]=model->Update_Chi_5(dt,V,this->Chi[i-1][5]);
	this->Chi[i][6]=model->Update_Chi_6(dt,V,this->Chi[i-1][6]);
	double Na=0;
	for(int j=0;j<7;j++){Na+=this->Chi[i][j];};
	model->Na_noise=Na;
};

void Solver :: Update_Psi(double dt,double V,int i){
	this->Psi[i][0]=model->Update_Psi_0(dt,V,this->Psi[i-1][0]);
	this->Psi[i][1]=model->Update_Psi_1(dt,V,this->Psi[i-1][1]);
	this->Psi[i][2]=model->Update_Psi_2(dt,V,this->Psi[i-1][2]);
	this->Psi[i][3]=model->Update_Psi_3(dt,V,this->Psi[i-1][3]);
	double K=0;
	for(int j=0;j<4;j++){K+=this->Psi[i][j];};
	model->K_noise=K;
};*/

double * Solver :: Calculator_aux(double y_1,double y_2,double y_3,double y_4,double y_5,double y_6,double t){
	double y[4] = {y_1,y_2,y_3,y_4};
	double s[2] = {y_5,y_6};
	return model->model_call(y,s,t);

};

void Solver :: Euler(){
	//cout<<V[0]<<" "<<n[0]<<" "<<m[0]<<" "<<h[0]<<" "<<t[0]<<endl;
	for(int i=1;i<=num;i++){
		//this->Update_Chi(dt,this->V[i-1],i);
		//this->Update_Psi(dt,this->V[i-1],i);
		double * K_1=Calculator_aux(V[i-1],n[i-1],m[i-1],h[i-1],s1[i-1],s2[i-1],t[i-1]);
		V[i]=(V[i-1]+K_1[0]*dt);
		n[i]=(n[i-1] +K_1[1]*dt);
		//m[i]=(m[i-1] + K_1[2]*dt);
		m[i]=model->constants->A_m(V[i])/(model->constants->A_m(V[i])+model->constants->B_m(V[i]));
		h[i]=(h[i-1]+K_1[3]*dt);
		if(model->eps!=0){
			//s1[i]=model->constants->SimplifiedCBM_s1(t[i],model->s1_init);
			s1[i]=(s1[i-1]+K_1[4]*dt);
			s2[i]=(s2[i-1]+K_1[5]*dt);
		}
		else{
			s1[i]=model->s1_const;
			s2[i]=model->s2_const;
		}	
		//cout<<V[i]<<" "<<n[i]<<" "<<m[i]<<" "<<h[i]<<" "<<t[i]<<endl;
	};


};


void Solver ::	Runge_Kutta(){
	cout<<V[0]<<" "<<n[0]<<" "<<m[0]<<" "<<h[0]<<" "<<t[0]<<endl;
	for(int i=1;i<num;i++){
		double * K_1=this->Calculator_aux(this->V[i-1],this->n[i-1],this->m[i-1],this->h[i-1],this->s1[i-1],this->s2[i-1],this->t[i-1]);
		double * K_2=this->Calculator_aux(this->V[i-1]+(0.5*K_1[0]*dt),this->n[i-1]+(0.5*K_1[1]*dt),this->m[i-1]+(0.5*K_1[2]*dt),
			this->h[i-1]+(0.5*K_1[3]*dt),this->s1[i-1]+(0.5*K_1[4]*dt),this->s2[i-1]+(0.5*K_1[5]*dt),this->t[i-1]);
		double * K_3=this->Calculator_aux(this->V[i-1]+(0.5*K_2[0]*dt),this->n[i-1]+(0.5*K_2[1]*dt),this->m[i-1]+(0.5*K_2[2]*dt),
			this->h[i-1]+(0.5*K_2[3]*dt),this->s1[i-1]+(0.5*K_2[4]*dt),this->s2[i-1]+(0.5*K_1[5]*dt),t[i-1]);
		double * K_4=this->Calculator_aux(this->V[i-1]+K_3[0]*dt,this->n[i-1]+K_3[1]*dt,this->m[i-1]+K_3[2]*dt,
			this->h[i-1]+K_3[3]*dt,this->s1[i-1]+K_3[4]*dt,this->s2[i-1]+K_3[5]*dt,this->t[i-1]);
	
		
		// 2 thread
		V[i]=(V[i-1]+(K_1[0]+2*K_2[0]+2*K_3[0]+K_4[0])*dt/6);
		n[i]=(n[i-1]+(K_1[1]+2*K_2[1]+2*K_3[1]+K_4[1])*dt/6);
		m[i]=(m[i-1]+(K_1[2]+2*K_2[2]+2*K_3[2]+K_4[2])*dt/6);
		//m[i]=model->constants->A_m(V[i])/(model->constants->A_m(V[i])+model->constants->B_m(V[i]));
		h[i]=(h[i-1]+(K_1[3]+2*K_2[3]+2*K_3[3]+K_4[3])*dt/6);
		//3 thread
		/*if(model->eps!=0){
			s1[i]=(s1[i-1]+(K_1[4]+2*K_2[4]+2*K_3[4]+K_4[4])*dt/6);
			s2[i]=(s2[i-1]+(K_1[5]+2*K_2[5]+2*K_3[5]+K_4[5])*dt/6);
		}
		else{*/
			s1[i]=(model->s1_const);
			s2[i]=(model->s2_const);
	//	}
		cout<<V[i]<<" "<<n[i]<<" "<<m[i]<<" "<<h[i]<<" "<<t[i]<<endl;
	};	
};



void Solver :: Temporal_Integral(std :: string funz,double s){
	float sum=0;
	//float constant=1;
	if(funz=="gamma"){
		//constant=(-model->constants->G_s1(V[0])+model->constants->G_s1(V[num]))/num;
		for(int i=0;i<=num;i++){sum+=(model->constants->G_s1(V[i]));};
	}
	if(funz=="delta"){
		//sum=0.5*(model->constants->D_s1(V[0])+model->constants->D_s1(V[num]));
		for(int i=0;i<=num;i++){sum+=model->constants->D_s1(V[i]);};
	}
	if(funz=="volt"){
		//sum=0.5*(V[0]+V[num]);
		for(int i=0;i<num;i++){sum+=V[i];};
	}
	cout<<sum/num<<" "<<model->s1_const<<endl;
};