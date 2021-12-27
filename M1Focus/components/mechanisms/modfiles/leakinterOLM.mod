: passive leak current

NEURON {
	SUFFIX leakinterOLM
	NONSPECIFIC_CURRENT il
	RANGE il, el, glbar_inter
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	glbar_inter = 1.3e-4 (siemens/cm2) < 0, 1e9 >
	el = -68 (mV)
}

ASSIGNED {
	v (mV)
	il (mA/cm2)
}

BREAKPOINT { 
	il = glbar_inter*(v - el)
}
