#!/usr/bin/python
import numpy as np
xexmax=6.5
lgfmin=-6.
name='ElvisCantelli.013122'
fil=open(name, 'r')
tbdata=np.genfromtxt(fil, delimiter=',', autostrip=True, names="species, lambda, chiex, loggf, rad, stark, waals, lande, ref", usecols=("species", "lambda", "chiex", "loggf", "waals", "lande"), dtype=("|S6, f4, f2, f2, f2, f2, f2, f2, |S100"))
fil.close()
spc=[s[0].strip('\'').replace(" ","").upper() for s in tbdata]
lbd=[l[1] for l in tbdata]
xex=[x[2] for x in tbdata]
lgf=[g[3] for g in tbdata]
wls=[w[4] for w in tbdata]
lde=[d[5] for d in tbdata]
ofil=open('atom-uves-580.dat', 'w')
for i in range(len(lbd)):
	if (xex[i] < xexmax) and (lgf[i] > lgfmin) and (lde[i] != 99.000) and (i != lbd.index(lbd[-1])) and (spc[i][-1] < '3'):
		if wls[i] != 0.000:
			ofil.write("%3s%9.3f"%(spc[i],lbd[i])+'    0.\n'+"%8.3f%8.3f%2s%3.3E"%(xex[i],lgf[i],"  ",10**(2.5*wls[i]-12.32))+" 0.00E+00 0.00E+00  0.5 1. 0\n")
		else:
			ofil.write("%3s%9.3f"%(spc[i],lbd[i])+'    0.\n'+"%8.3f%8.3f"%(xex[i],lgf[i])+"  0.300E-31 0.00E+00 0.00E+00  0.5 1. 0\n")
	elif (xex[i] < xexmax) and (lgf[i] > lgfmin) and (lde[i] != 99.000) and (i == lbd.index(lbd[-1])) and (spc[i][-1] < '3'):
		if wls[i] != 0.000:
			ofil.write("%3s%9.3f"%(spc[i],lbd[i])+'    0.\n'+"%8.3f%8.3f%2s%3.3E"%(xex[i],lgf[i],"  ",10**(2.5*wls[i]-12.32))+" 0.00E+00 0.00E+00  0.5 1. 1\n")
		else:
			ofil.write("%3s%9.3f"%(spc[i],lbd[i])+'    0.\n'+"%8.3f%8.3f"%(xex[i],lgf[i])+"  0.300E-31 0.00E+00 0.00E+00  0.5 1. 1\n")
ofil.close()
