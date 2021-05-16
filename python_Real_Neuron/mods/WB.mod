TITLE WB.mod   squid sodium, potassium, and leak channels
 
COMMENT
Wang Buzsaki 
ENDCOMMENT
 
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
	SUFFIX WB
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el, NNa, NK, se, m, h, n
}
 
PARAMETER {
    se = -1
	gnabar = .08 (S/cm2)	<0,1e9>
	gkbar = .036 (S/cm2)	<0,1e9>
	gl = .0003 (S/cm2)	<0,1e9>
	el = -54.3 (mV)

}
 
ASSIGNED {
	v (mV)
	celsius (degC)
	ena (mV)
	ek (mV)
	dt (ms)
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	am	(/ms)
	ah	(/ms)
	an	(/ms)
	bm	(/ms)
	bh	(/ms)
	bn	(/ms)

}
 
STATE {	
	n
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    m=am/(am+bm)
	ina = gnabar*(m^3*h)*(v - ena)
	ik = gkbar*(n^4)*(v - ek)
	il = gl*(v - el)
}
 
INITIAL {
	func(v)
    if (se>0) {set_seed(se)} 
    n=an/(an+bn)
	m=am/(am+bm)
	h=ah/(ah+bh)
}

DERIVATIVE states {  
	LOCAL flag, n0, h0
	func(v)

	flag=0
	n0=n
    : If the variable leaves [0,1] then random numbers are drawn again
	while (flag==0) {
		n'=2*(an*(1-n)-bn*n)
		if (n<0 || n>1) {
			n=n0
		} else {flag=1}
	}
	
	flag=0
	h0=h
	while (flag==0) {
		h'=2*(ah*(1-h)-bh*h)
		if (h<0 || h>1) {
			h=h0
		} else {flag=1}
	}

}


PROCEDURE func(v(mV)) {  :Computes rate and other constants at current v.
	LOCAL q10
	UNITSOFF

	am = 0.1*(v+40)/(1-exp(-(v+40)/10))
	bm = 4*exp(-(v+65)/18)
	ah = 0.07*exp(-(v+65)/20) 
	bh = 1/(1+exp(-(v+35)/10))
	an = 0.01*(v+55)/(1-exp(-(v+55)/10))
	bn = 0.125*exp(-(v+65)/80)

	UNITSON 
}
