COMMENT
-----------------------------------------------------------------------------
ChR2H134R.mod

Model of Channelrhodopsin-2 (mutant H134R)
==========================================

from: 
Williams JC, Xu J, Lu Z, Klimas A, Chen X, et al. (2013) Computational Optogenetics: Empirically-Derived Voltage- and Light-Sensitive Channelrhodopsin- 2 Model. 
PLoS Comput Biol 9(9): e1003220. doi:10.1371/journal.pcbi.1003220

Original (MATLAB) code by: John C. Williams

Implemented by: Michele Giugliano, SISSA Trieste, 8/8/2018, mgiugliano@gmail.com

This mechanism includes the voltage, light, and temperature dependences of CHR2-H134R. 
Set the desired temperature by the hoc assignment statement ==> e.g. celsius = 37
-----------------------------------------------------------------------------
ENDCOMMENT

TITLE Channelrhodopsin-2 (mutant H134R) current density 


UNITS {
: Convenient aliases for the units...
	(mS) = (millisiemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}


NEURON {
: Public interface of the present mechanism...
	SUFFIX ChR2
	NONSPECIFIC_CURRENT i

	RANGE i, gmax, Irradiance
    RANGE light_intensity, light_delay, pulse_width, n,wavelength

    GLOBAL A, B, C, gamma
	GLOBAL hc, wloss, sigma_retinal
	GLOBAL Q10_Gd1, Q10_Gd2, Q10_Gr, Q10_e12dark, Q10_e21dark, Q10_epsilon1, Q10_epsilon2
}



PARAMETER {
: Below are constants, variables that are not changed within
: this mechanism, and other variables changed by the user through the hoc code...

    n=3.                                      : number of light pulse
	gmax  			= 0.4  		(mS/cm2)	  : maximal conductance
	EChR2 			= 0.     	(mV)          : reversal potential

	light_delay     = 100.		(ms)		  : initial delay, before switching on the light pulse
	pulse_width     = 100.		(ms)		  : width of the light pulse
	light_intensity = 5.					  : mW/mm^2, intensity of the light pulse
    Irradiance      = 0.
    gamma           = 0.05					  : ratio of conductances in states O2/O1, unit-less
	A 				= 10.6408   (mV)          : Be careful with implementing eqs. 1 and 12!
	B 				= -14.6408  (mV)		  :
	C 				= -42.7671  (mV)          :

	wavelength 		= 594		    		  : wavelength of max absorption for retinal, nm
	hc       		= 1.986446E-25  		  : Planck's constant * speed of light, kg m^3/s^2
	wloss    		= 1.3      : scaling factor for losses of photons due to scattering or absorption
	sigma_retinal 	= 12.E-20       		  : retinal cross-sectional area, m^2

	tauChR2  		= 3.2		(ms)          : time constant for ChR2 activation

	temp 			= 22	  (degC)		  : original temperature
	Q10_Gd1      	= 1.97					  : Q10 value for the temperature sensitivity
	Q10_Gd2      	= 1.77					  : Q10 value for the temperature sensitivity
	Q10_Gr       	= 2.56					  : Q10 value for the temperature sensitivity
	Q10_e12dark  	= 1.1					  : Q10 value for the temperature sensitivity
	Q10_e21dark  	= 1.95					  : Q10 value for the temperature sensitivity
	Q10_epsilon1 	= 1.46					  : Q10 value for the temperature sensitivity
	Q10_epsilon2 	= 2.77					  : Q10 value for the temperature sensitivity
    
    
}


ASSIGNED {
: variables calculated by the present mechanism or by NEURON
	v      		(mV)			: Membrane potential
	celsius		(degC)          : Temperature

	i    		(mA/cm2)		: Membrane current
    
	Gd1    		(1./ms)			: rate constant for O1->C1 transition
	Gd2    		(1./ms)			: rate constant for O2->C2 transition
	Gr     		(1./ms)			: rate constant for C2->C1 transition
	e12    		(1./ms)			: rate constant for O1->O2 transition
	e21    		(1./ms)			: rate constant for O2->O1 transition

	epsilon1 					: quantum efficiency for photon absorption from C1
	epsilon2 					: quantum efficiency for photon absorption from C2
	F   		(1./ms)			: photon flux: number of photons per molecule per second
	S0							: time- & irradiance-dep. activation func. (post-isomerization)
}


STATE {
: Let's declare the state variables
: Note: in any kinetic scheme with N state, you have N-1 *independent* state variables and
: this is the reason why C2 is not defined below, but rather expressed as (1 - O1 - O2 - C1). 
:
: O1, O2, C1 and C2 are the fractions of channels in open and closed states, while p is
: an additional state variable, capturing the kinetics of ChR2 activation 

	O1 O2 C1 p 
}


BREAKPOINT {
	SOLVE states METHOD cnexp

	Irradiance = 0.               : typically with values in the range 0 - 10 mW/mm^2
   
    : The control of the Irradiance waveform within Neuron is very rudimentary and
    : it is made explicit below. In other words, the time course of Irradiance is 
    : upfront defined - in this example - as a single-pulse photoactivation protocol.
    : Modify it to fit your needs!
COMMENT
    if (t < light_delay)                      { Irradiance = 0. }
    else if (t <= (light_delay + pulse_width)) { Irradiance = light_intensity }
    else if ( t >= (light_delay*2) && t <=(light_delay*2 + pulse_width) ) { Irradiance = light_intensity }
    else if(t >(light_delay*2 + pulse_width)) { Irradiance = 0. }

ENDCOMMENT

    FROM i = 0 TO n BY 1{
        if ( t >= (light_delay*i) && t <=(light_delay*i + pulse_width) ){ Irradiance = light_intensity }
    }


	i = (0.001) * (gmax * (A + B * exp(v /C))/v * (O1 + gamma * O2) * (v - EChR2))

	: note: the above expression includes the voltage-dep. rectification of ChR2
}


INITIAL {
: Let's set the state variables to their initial values... 
	rates(v)
	C1 = 1.
	O1 = 0.
	O2 = 0.
	p  = 0.
   
}


DERIVATIVE states {
: The state variables are computed here...
    rates(v)
    : Note: these are the full equations:
	:O1' = -(Gd1+e12)           * O1 + e21 * O2 + epsilon1*F*p * C1
	:O2' = -(Gd2+e21)           * O2 + e12 * O1 + epsilon2*F*p * C2
	:C1' = -epsilon1*F*p        * C1 + Gd1 * O1 + Gr  * C2   
	:C2' = -(epsilon2*F*p + Gr) * C2 + Gd2 * O2
	:p'  =  (S0 - p) / tauChR2 

	: However, only 3 of them are independent. Let's then define C2 as (1. - C1 - O1 - O2)
	O1' = -(Gd1+e12)           * O1 + e21 * O2 + epsilon1*F*p * C1
	O2' = -(Gd2+e21)           * O2 + e12 * O1 + epsilon2*F*p * (1. - C1 - O1 - O2)
	C1' = -epsilon1*F*p        * C1 + Gd1 * O1 + Gr  * (1. - C1 - O1 - O2)  

	p'  =  (S0 - p) / tauChR2
}


UNITSOFF
PROCEDURE rates(v (mV)) {
    LOCAL e12dark, e21dark, logphi0, Ephoton, flux

	: Values at 22 degrees celsius...
	e12dark  = 0.011                                 : ms^-1
	e21dark  = 0.008                                 : ms^-1
	epsilon1 = 0.8535
	epsilon2 = 0.14
	Gd1 = 0.37    	 : dark-adapted deactivation rate, ms^-1
	Gd2 = 0.01                                  	 : ms^-1
	Gr  = 0.00000667    	 : recovery rate ms^-1

	if (Irradiance>0) {
		logphi0  = log(1. + Irradiance / 0.024)      : unit-less
	}
	else {
		logphi0 = 0.	 							 : for consistency
	}
	e12      = e12dark + 0.005 * logphi0             : ms^-1
	e21      = e21dark + 0.004 * logphi0             : ms^-1
    

	S0       = 0.5 * (1. + tanh(120.*(100. * Irradiance - 0.1))) : unit-less
	Ephoton  = 1.E9 * hc / wavelength          : J, scaling wavelength from nm to m	
	flux     = 1000. * Irradiance / Ephoton    : 1/(s*m^2), scaling irradiance from mW/mm^2 to W/m^2
	F        = flux  * sigma_retinal / (wloss * 1000.) : ms^-1, scaling F from 1/s to 1/ms
}
UNITSON
