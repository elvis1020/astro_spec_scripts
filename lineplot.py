import numpy
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
from matplotlib import rc

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

rc('font', family = 'serif', serif = 'cmr10')
rc('xtick', labelsize=14) 
rc('ytick', labelsize=14) 
command = 'x'
while command != 'q':
	lines='line_input.dat'
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
	with open(lines) as l:
		ln = l.readlines()[1:]
	for cols in (raw.strip().split() for raw in ln):
		wlobs.append([])
		fluxobs.append([])
		wlmod.append([])
		fluxmod.append([])
		if cols[0][-5:] != '.fits':
			with open(cols[0]) as o:
				obs = o.readlines()
			for coll in (raw.strip().split() for raw in obs):
				if is_number(coll[0]):
					wlobs[count].append(float(coll[0]))
					fluxobs[count].append(float(coll[1]))
		else:
			wlobs[count],fluxobs[count]=readfits(cols[0])
		if (cols[1][0] != '@') and (cols[1] != 'none'):
			multimod.append(1)
			if cols[1][-5:] != '.fits':
				with open(cols[1]) as m:
					mod = m.readlines()
				for coll in (raw.strip().split() for raw in mod):
					if is_number(coll[0]):
						wlmod[count].append(float(coll[0]))
						fluxmod[count].append(float(coll[1]))
			else:
				wlmod[count],fluxmod[count]=readfits(cols[1])
		elif cols[1] == 'none':
			multimod.append('x')
			wlmod[count],fluxmod[count]=[],[]
		else:
			multimod.append(0)
			with open(cols[1][1:]) as lf:
				lst = lf.readlines()
			for coll in (raw.strip().split() for raw in lst):
				multimod[-1] += 1
				wlmod[count].append([])
				fluxmod[count].append([])
				if coll[0][-5:] != '.fits':
					with open(coll[0]) as spf:
						mltspc = spf.readlines()
					for colu in (raw.strip().split() for raw in mltspc):
						if is_number(colu[0][0]):
							wlmod[count][multimod[-1]-1].append(float(colu[0]))
							fluxmod[count][multimod[-1]-1].append(float(colu[1]))
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
	fig=plt.figure(figsize=(12,(2.6*(round(count/2.))+1.12)))
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
		oblinex_raw=(wlobs[line][o1:o2])
		obliney_raw=(fluxobs[line][o1:o2])
		oblinex=[lb+dspl[line] for lb in oblinex_raw]
		obliney=[fx/norm[line] for fx in obliney_raw]
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
		plt.plot(oblinex,obliney,stl[line])
		if multimod[line] == 1:
			plt.plot(modlinex,modliney,'r-')
		else:
			if multimod[line] != 'x':
				colors = plt.get_cmap(clr[line])(numpy.linspace(0, 1.0, multimod[line]))
				for m in range(multimod[line]):
					plt.plot(modlinex[m],modliney[m], color=colors[m])
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
