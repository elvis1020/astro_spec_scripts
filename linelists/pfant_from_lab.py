#!/usr/bin/python
import numpy as np
xexmax=6.
lgfmin=-5.
name='linelist_elvis'
#tbdata=np.genfromtxt(name, unpack=True)

spc = []
lbd = []
xex = []
lgf = []
c6  = []

fil=open(name, 'r')
for cols in (raw.strip().split() for raw in fil):
	spc.append(cols[1])
	lbd.append(float(cols[0]))
	xex.append(float(cols[2]))
	lgf.append(float(cols[4]))
	c6.append(cols[5])
ofil=open('atomgi.dat', 'w')
for i in range(len(lbd)):
	ofil.write("%3s%9.3f"%(spc[i],lbd[i])+'    0.\n'+"%8.3f%8.3f%2s%9s"%(xex[i],lgf[i],"  ",c6[i])+" 0.00E+00 0.00E+00  0.5 1. 0\n")
ofil.close()
