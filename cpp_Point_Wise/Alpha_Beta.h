#ifndef __Alpha_Beta_H__
#define __Alpha_Beta_H__


class Alpha_Beta{

private:

	float  vtrap( float x,float y);

public:

	Alpha_Beta();
	float G_s1(float V);	//Hz

	float G_s2(float V);   //HZ

	float D_s1(float V);	//Hz

	float D_s2(float V);   //HZ

	//WB
	float A_n( float V);  //kHz
	float A_m(float V);  //kHz

	//HH
	float A_n_HH( float V);  //kHz
	float A_m_HH(float V);  //kHz

	float A_h(float V);   //kHz
	float B_n(float V); 
	float B_m(float V);
	float B_h(float V);

	float SimplifiedCBM_s1(double t,double init);
	float gamma_integral(double Periodicity,float * V);
	float delta_integral(double Periodicity,float * V);
};


#endif
