TITLE T-calcium channel



UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}



NEURON {
	SUFFIX cat
	USEION ca READ cai,cao WRITE ica
        RANGE gcat, gcatbar, hinf, minf, taum, tauh, iCa
}

PARAMETER {
	v (mV)
	
      tBase = 23.5  (degC)
	celsius = 22  (degC)
	gcatbar = 0   (mho/cm2)  : initialized conductance
	ki = 0.001    (mM)
        tfa = 1                  : activation time constant scaling factor
       : tfi = 0.68 
        tfi = 0.68               : inactivation time constant scaling factor
        eca = 140                : Ca++ reversal potential

}


STATE {
	m h 
}

ASSIGNED {
          ica (mA/cm2)
	  gcat (mho/cm2)
	  hinf
	  tauh
	  minf
	  taum
	  cai (mM)      : initial internal Ca++ concentration
	  cao (mM)      : initial external Ca++ concentration

}

INITIAL {
	rates(v)
	m = minf
	h = hinf
     gcat = gcatbar*m*m*h*h2(cai)

}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcat = gcatbar*m*m*h*h2(cai)
	ica = gcat*ghk(v,cai,cao)

}

DERIVATIVE states {	: exact when v held constant
	rates(v)
	m' = (minf - m)/taum
	h' = (hinf - h)/tauh
}


UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (DegC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}


FUNCTION alph(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alph = 1.6e-4*exp(-(v+57)/19)
}

FUNCTION beth(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
	:beth = 1/(exp((-v+15)/10)+1.0)
      beth = 1/(exp((-v+15)/10)+1.0)
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)
}

FUNCTION betm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	betm = 0.046*exp(-v/22.73)
}


PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a
        a = alpm(v)
        taum = 1/(tfa*(a + betm(v))) : estimation of activation tau
        minf =  a/(a+betm(v))        : estimation of activation steady state
        a = alph(v)
        tauh = 1/(tfi*(a + beth(v))) : estimation of inactivation tau
        hinf = a/(a+beth(v))         : estimation of inactivation steady state
        
}
