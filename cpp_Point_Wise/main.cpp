#include "Integrator.h"
#include "Stoch_Integrator.h"
#include "Statistics.h"
#include <iostream>
#include <thread> 
#include <cstdlib> 
#include <random>
#include <fstream>

using namespace std;
void GetnumSpikes(double T,double I){
	for(int i=0;i<10;i++){
	Solver * solv=new Solver(I,T);
	solv->Euler();
	Stat * calc=new Stat(I,solv->T,solv->num,solv->t,solv->Period,solv->V,solv->s1);
	cout<<calc->Peak()<<endl;
	}
};

double ReadLastvalues(Solver_Stoch * solv){
	ifstream ifile("res.txt", std::ios::in);
	float * vals=new float[10];
	float num; 
	int i=0;
    while (ifile >> num) {
        vals[i]=num;
        i++;
    }
    solv->V[0]=vals[0];
    solv->n[0]=vals[1];
    solv->m[0]=vals[2];
    solv->h[0]=vals[3];
    solv->s1[0]=vals[4];
    solv->q_Na[0]=vals[5];
    solv->q_K[0]=vals[6];
    solv->p_Na[0]=vals[7];
    solv->p_K[0]=vals[8];
    solv->t[0]=vals[9];
    solv->model->s1_init=vals[4];
    return vals[9];
};

int main(int argc,char* argv[]){
	double T=atof(argv[1]);   //ms periodicity
	double I=atof(argv[2]);
	/*thread th1(GetnumSpikes, T,I); 
	thread th2(GetnumSpikes, T,I); 
	thread th3(GetnumSpikes, T,I); 
	thread th4(GetnumSpikes, T,I); 
	th1.join();
	th2.join();
	th3.join();
	th4.join();*/

	Solver_Stoch * solv=new Solver_Stoch(I,T);
	//double start = ReadLastvalues(solv);
	//solv->InitTime(start);

	// CALL THE SIMULATOR FOR INTEGRATION
	//solv->Euler();
	solv->Euler_Guller();
	//solv->Predictor_Corrector();
	//solv->Runge_Kutta();
	Stat * calc=new Stat(I,solv->T,solv->num,solv->t,solv->Period,solv->V,solv->s1);
	//calc->Mean_Firing_rate();
	/*for(int s=1;s<=42;s+=5){
	//	solv->Period=1000/s;
		Solver * solv=new Solver(I,1000/s);
		solv->Euler();
		Stat * calc=new Stat(I,solv->T,solv->num,solv->t,solv->Period,solv->V,solv->s1);
		cout<<calc->Mean_Firing_rate()<<" "<<s<<endl;
		//calc->Latency();
		//cout<<" "<<s<<endl;
		//free(calc);
	}*/
	//calc->I_vect();
	//calc->Mean_s();

	// Calculate the latency
	calc->Latency(0);
	//calc->Peak();
	
	/*for(int s=0;s<=15;s++){
		solv->model->s1_const=1;
		solv->Euler();
		solv->Temporal_Integral("gamma",1);
	}*/
	//free(solv);

};
