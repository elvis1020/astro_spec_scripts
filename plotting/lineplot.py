#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
from matplotlib import rc
import scipy.optimize as opt
import scipy.interpolate as itp

def bissec(A,x):
	n = len(A)-1
	if (x < A[0]):
		return 0
	elif (x > A[n]):
		return n+1
	n1 = 0
	n2 = n
	while ((n2-n1)>1):
		nm = (n1+n2)/2
		if ((x-A[nm]) > 0):
			n1 = nm
		else:
			n2 = nm
	return n1

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def readfits(file):
	import pyfits
	spectra = pyfits.open(file)
	try:
		bin = spectra[0].header['CDELT1']
	except:
		bin = spectra[0].header['CD1_1']
		print bin
	lambda0 = spectra[0].header['CRVAL1']
	flux = spectra[0].data
	x=[lambda0+i*bin for i in range(len(flux))]
	y=[i for i in flux]
	return x,y

def adjspec(x, flcor, lbcor):
	modf = itp.interp1d(wlmod[line], fluxmod[line])
	return modf(x-lbcor)*flcor


if len(sys.argv) > 1:
	inptfil = sys.argv[1]
else:
	inptfil = 'line_input.dat'

rc('font', family = 'serif', serif = 'cmr10')
rc('xtick', labelsize=14) 
rc('ytick', labelsize=14) 
command = 'x'
while command != 'q':
	lines=inptfil
	wlobs=[]
	fluxobs=[]
	wlmod=[]
	fluxmod=[]
	lnlst=[]
	elm=[]
	dspl=[]
	norm=[]
	ran=[]
	clr=[]
	stl=[]
	count=0
	multimod=[]
	flagfit = False
	with open(lines) as l:
		ln = l.readlines()[1:]
	for cols in (raw.strip().split() for raw in ln):
		if cols[0][0] != '#':
			wlobs.append([])
			fluxobs.append([])
			wlmod.append([])
			fluxmod.append([])
			if cols[0][-5:] != '.fits':
				obs = np.genfromtxt(cols[0], unpack=True)
				wlobs[count] = obs[0]
				fluxobs[count] = obs[1]
			else:
				wlobs[count],fluxobs[count]=readfits(cols[0])
			if (cols[1][0] != '@') and (cols[1] != 'none'):
				multimod.append(1)
				if cols[1][-5:] != '.fits':
					mod = np.genfromtxt(cols[1], unpack=True)
					wlmod[count] = mod[0]
					fluxmod[count] = mod[1]
				else:
					wlmod[count],fluxmod[count]=readfits(cols[1])
				if command == 'fit':
					flagfit = True
			elif cols[1] == 'none':
				multimod.append('x')
				wlmod[count],fluxmod[count]=[],[]
			else:
				multimod.append(0)
				lst = open(cols[1][1:])
				for coll in (raw.strip().split() for raw in lst):
					multimod[-1] += 1
					wlmod[count].append([])
					fluxmod[count].append([])
					if coll[0][-5:] != '.fits':
						mltspc = np.genfromtxt(coll[0], unpack=True)
						wlmod[count][multimod[-1]-1] = mltspc[0] 
						fluxmod[count][multimod[-1]-1] = mltspc[1]
					else:
						wlmod[count][multimod[-1]-1],fluxmod[count][multimod[-1]-1]=readfits(coll[0])
			elm.append(cols[2])
			lnlst.append(float(cols[3]))
			dspl.append(float(cols[4]))
			norm.append(float(cols[5]))
			ran.append(float(cols[6]))
			clr.append(cols[7])
			stl.append(cols[8])
			count += 1

#________plotting start________#

	fignum=0
	fig=plt.figure(figsize=(17,(3.6*(round(count/2.))+1.12)))
	for line in range(len(lnlst)):
		o1=bissec(wlobs[line],lnlst[line]-ran[line]-0.5)
		o2=bissec(wlobs[line],lnlst[line]+ran[line]+0.5)
		if multimod[line] == 1:
			m1=bissec(wlmod[line],lnlst[line]-ran[line]-0.5)
			m2=bissec(wlmod[line],lnlst[line]+ran[line]+0.5)
		else:
			m1=[]
			m2=[]
			if multimod[line] != 'x':
				for m in range(multimod[line]):
					m1.append(bissec(wlmod[line][m],lnlst[line]-ran[line]-0.5))
					m2.append(bissec(wlmod[line][m],lnlst[line]+ran[line]+0.5))
			else:
				m1.append(0)
				m2.append(0)
		oblinex = np.array(wlobs[line][o1:o2]) + dspl[line]
		obliney = np.array(fluxobs[line][o1:o2]) / norm[line]
		#oblinex=[lb+dspl[line] for lb in oblinex_raw]
		#obliney=[fx/norm[line] for fx in obliney_raw]
		if multimod[line] == 1:
			modlinex=(wlmod[line][m1:m2])
			modliney=(fluxmod[line][m1:m2])
		else: 
			modlinex=[]
			modliney=[]
			if multimod[line] != 'x':
				for m in range(multimod[line]):
					modlinex.append(wlmod[line][m][m1[m]:m2[m]])
					modliney.append(fluxmod[line][m][m1[m]:m2[m]])
			else:
				modlinex.append(0)
				modliney.append(0)
		A=int(fignum/2.)*(0.88*(round(count/2.))**(-0.95))
		B=(0.2*(round(count/2.))**(-0.61))
		ysize=0.741*(round(count/2.))**(-1.02)
		ypos=B+(1-B-ysize-A-0.017)
		fig.add_axes([0.09+(fignum%2)*0.47, ypos, 0.40, ysize])
		plt.ion()
		
		if flagfit:
			initial_guess = (1., 0.)
			popt, cov = opt.curve_fit(adjspec, oblinex, obliney, p0=initial_guess)
			plt.text(0.98, 0.06,popt[0],horizontalalignment='right',verticalalignment='bottom', fontsize=10, transform = plt.gca (). transAxes)
			plt.text(0.98, 0.16,popt[1],horizontalalignment='right',verticalalignment='bottom', fontsize=10, transform = plt.gca (). transAxes)
			
		plt.plot(oblinex,obliney,stl[line])
		plt.plot([oblinex[0], oblinex[-1]], [1, 1], 'b--', linewidth=0.5)
		if multimod[line] == 1:
			plt.plot(modlinex,modliney,'r-')
		else:
			if multimod[line] != 'x':
				colors = plt.get_cmap(clr[line])(np.linspace(0, 1.0, multimod[line]+2))
				for m in range(multimod[line]):             
					plt.plot(modlinex[m],modliney[m], color=colors[m+1])
		if elm[line] != 'none':
			if  elm[line][0] != '#':
				linename=elm[line] + ' ' + str(lnlst[line])
			else:
				linename=elm[line][1:]
			plt.text(0.04, 0.06,linename,horizontalalignment='left',verticalalignment='bottom', fontsize=16, transform = plt.gca (). transAxes)
		if multimod[line] == 1:
			modmin=min(modliney)
			modmax=max(modliney)
		else:
			modmin=min(obliney)
			modmax=max(obliney)
			if multimod[line] != 'x':
				for m in range(multimod[line]):
					if min(modliney[m]) < modmin:
						modmin = min(modliney[m])
					if max(modliney[m]) > modmax:
						modmax = max(modliney[m])
		plt.axis([lnlst[line]-ran[line], lnlst[line]+ran[line], (min([min(obliney),modmin]))-(max([max(obliney),modmax]))*0.05, (max([max(obliney),modmax]))+(max([max(obliney),modmax]))*0.05])
		plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
		mjtickx=(float(round((ran[line]*2./5.)*10)))/10
		if mjtickx == 0:
			mjtickx=(float(round((ran[line]*2./4.)*100)))/100
		mntickx=mjtickx/5
		plt.gca().xaxis.set_major_locator(MultipleLocator(mjtickx))
		plt.gca().xaxis.set_minor_locator(MultipleLocator(mntickx))
		yrange=(max([max(obliney),modmax]))*1.05 - (min([min(obliney),modmin]))*0.95
		plt.gca().yaxis.set_major_locator(MultipleLocator(yrange/4))
		plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%1.2f'))
		fignum += 1
	plt.figtext(0.5,0.02,"Wavelength (angstroms)",fontdict={'fontsize':18}, horizontalalignment='center')
	plt.figtext(0.01,0.5,"Arbitrary Flux",fontdict={'fontsize':18},rotation=90, verticalalignment='center')
	plt.show()
	command=raw_input('enter to refresh; q to quit\n')
	plt.close()
