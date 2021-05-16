#ifndef __Statistics_H__
#define __Statistics_H__
#include <vector>
#include "In_current.h"

class Stat{
public:
	std::vector<float> t; 
	float * s;
	double I0;
	int T;
	int num;
	double Periodicity;
	float * V;
	In_current * current;
	Stat(double I,int tot,int nume,float * time,double Period,float * Volt,float * sod);
	void Mean_s();
	float Mean_Firing_rate();
	float I_vect();
	int max_index(int n,int k);
	void Latency(double start);
	int  Peak();

};
#endif