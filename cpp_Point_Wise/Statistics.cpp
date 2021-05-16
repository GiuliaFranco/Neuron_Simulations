#include "Statistics.h"
#include "In_current.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <fstream>
using namespace std;

Stat ::Stat(double I,int tot,int nume,float * time,double Period,float * Volt,float * sod){
	I0=I;
	T=tot;
	num=nume;
	for(int i=0;i<num;i++){t.push_back(time[i]);}
	Periodicity=Period;
	V=Volt;
	s=sod;
	current=new In_current(I0,T,Periodicity);
};



void Stat :: Mean_s(){
	int tot_sec=T/1000.;
	int tot_val=T*100/tot_sec; 
	cout<<s[0]<<" "<<t[0]<<endl;   
	for(int i=0;i<tot_sec;i++){
		float sum=0;
		for(int k=0;k<tot_val;k++){
			sum+=s[k+i*tot_val];
		}
		cout<<sum/tot_val<<" "<<i+1<<endl;;
	};

};

float Stat :: I_vect(){
	std::vector<float> curr;
	//int  * Npeak=this->Peak();
	for(int i=0;i<num;i++){
		curr.push_back(current->Get_current(t[i])/I0);
		cout<<V[i]<<" "<<curr[i]<<" "<<t[i]/1000<<endl;
	};
	return curr.size();
};

int Stat :: max_index(int n,int k){
	float * sub_V=new float[n];
	for(int j=0;j<n;j++){
		sub_V[j]=V[k+j];
	}
    if(n <= 0) return -1;
    int i, max_i = 0;
    float max = sub_V[0];
    for(i = 1; i < n; ++i){
        if(sub_V[i] > max){
            max = sub_V[i];
            max_i = i;
            }
        }
    if(max>=0) return max_i+k;
    else return 0;
    };

void Stat :: Latency(double start){
	//ofstream myfile;
	//myfile.open ("Latency.txt");
	double delta=current->delta_t;
	int k=T/Periodicity;
	int sub=(Periodicity-delta)*100; //1000.01
	int sub1=(Periodicity-delta)*100;
	std::vector<float> Latency_val;
	for(int i=1;i<k;i++){
		float current_time=0;
		if(sub1-1<=t.size()){ current_time=t[sub1-1];}
		int V_indx=max_index(sub,sub1);
		if(V_indx!=0){
			float V_time=t[V_indx];
			Latency_val.push_back(fabs(V_time-current_time));
			// Latency, current time, spike time, spike amplitude
			cout<<fabs(V_time-current_time)<<" "<<t[sub1-1]/1000<<" "<<V_time/1000<<" "<<V[V_indx]<<endl;
		}
		sub1+=Periodicity*100;
	};
	//myfile.close();
	//float average = std:: accumulate( Latency_val.begin(), Latency_val.end(), 0.0) /Latency_val.size();
	//cout<<average;
};

float Stat :: Mean_Firing_rate(){
	//cout<<this->Peak()<<endl;
	//cout<<T/Periodicity<<endl;
	cout<<this->Peak()/(T/Periodicity)<<endl;
	return (1000/Periodicity)*this->Peak()/(T/Periodicity);
};

int  Stat :: Peak(){
	int count=0;
	std::vector<int> Npeak; 
	int Na=0;
	for(int i=0;i<num;i++){
		Npeak.push_back(0);
		if(count==0 && V[i]>=-10.){
			count=1;		
			Na+=1;
			Npeak.push_back(1);
		} 
		if(count==1 && V[i]<-10.){
			count=0;	
			Npeak.push_back(0);
		}
		//cout<<Npeak[i]<<" "<<t[i]<<endl;
	}
	return Na;
};

