#include "In_current.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

In_current :: In_current(double I,int tot,double Periodicity){
	I_0=I;   
	T=Periodicity;
	T_tot=tot;  
	k=T_tot/T;
	//i=0;
};

double In_current :: Get_current(double t){
	for(int i=1;i<=k;i++){
		//cout<<t<<" "<<((i*T)+1000.01)<<endl;
		if(fabs(t-i*T)<=delta_t){
			return I_0;
		}
	}
	return 0;

	
	/*if(t>delta_t+i*T && i<=k){
		i+=1;
	}
	if(fabs(t-i*T)<=delta_t && i<=k){
			return I_0;
		}
	return 0;*/
	
};
