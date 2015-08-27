##########################################################
# Multiplot: a simple, interactive tool for spectroscopy #
# Written by E. Cantelli, 2015                           #
# elvis.cantelli@usp.br                                  #
##########################################################

import matplotlib.pyplot as plt
import matplotlib
from numpy import exp, pi, sqrt
import numpy as np
from lmfit.models import VoigtModel, GaussianModel, LinearModel, ConstantModel
from scipy.special import wofz
from scipy.interpolate import splrep,splev
from scipy import integrate
import math
import os

####################################### Function definition

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
	lambda0 = spectra[0].header['CRVAL1']
	flux = spectra[0].data
	x=[lambda0+i*bin for i in range(len(flux))]
	y=[i for i in flux]
	return x,y

def onclick(event):
    xval,yval = event.xdata, event.ydata
    dift=9999
    idx=0
    for x in range(len(MX)):
    	    xloc=bissec(MX[x], xval)
    	    if math.sqrt((xval-MX[x][xloc])**2+(yval-MY[x][xloc])**2)<dift:
    	    	    dift=math.sqrt((xval-MX[x][xloc])**2+(yval-MY[x][xloc])**2)
    	    	    idx=x
    print xval,yval
    print A[idx]

####################################### Measure EW function definition

s2pi = sqrt(2*pi)
s2   = sqrt(2.0)

def linear(x, slope, intercept):
    return slope * x + intercept

def const(x, c):
    return c

def voigt(x, amplitude, center, sigma, gamma):
    if gamma is None:
        gamma = sigma
    z = (x-center + 1j*gamma)/ (sigma*s2)
    return amplitude*wofz(z).real / (sigma*s2pi)

def gaussian(x, amplitude, center, sigma):
    return (amplitude/(s2pi*sigma)) * exp(-(x-center)**2 /(2*sigma**2))

def call_cont(x, y):
	cont = LinearModel(prefix='cont_')
	pars = cont.guess(y, x=x)
	return cont, pars

def call_constant(x, y, ylock):
	const = ConstantModel(prefix='constant_')
	pars = const.guess(y, x=x)
	pars['constant_c'].set(ylock, min=ylock-0.01, max=ylock+0.0001)
	return const, pars

def call_gauss(x, y, cen, count, pars):
	label='g'+str(count)+'_'
	gauss = GaussianModel(prefix=label)
	pars.update(gauss.make_params())
	pars[label+'center'].set(cen, min=cen-0.01, max=cen+0.01)
	pars[label+'amplitude'].set(-0.5, min=-10., max=0.0001)
	pars[label+'sigma'].set(0.1, min=0.005, max=0.25)
	return gauss

def call_voigt(x, y, cen, count, pars):
	label='v'+str(count)+'_'
	voigt = VoigtModel(prefix=label)
	pars.update(voigt.make_params())
	pars[label+'center'].set(cen, min=cen-0.01, max=cen+0.01)
	pars[label+'amplitude'].set(-0.5, min=-10., max=0.0001)
	pars[label+'sigma'].set(0.1, min=0.005, max=0.25)
	pars[label+'gamma'].set(value=0.7, vary=True, expr='')
	return voigt

def build_voigt(number, fpars, x):
	label='v'+str(number)+'_'
	c=fpars['constant_c']
	lbd=fpars[label+'center']
	amp=fpars[label+'amplitude']
	sgm=fpars[label+'sigma']
	gma=fpars[label+'gamma']
	y=[const(i, c)+voigt(i, amp, lbd, sgm, gma) for i in x]
	return y

def build_gauss(number, fpars, x):
	label='g'+str(number)+'_'
	c=fpars['constant_c']
	lbd=fpars[label+'center']
	amp=fpars[label+'amplitude']
	sgm=fpars[label+'sigma']
	y=[const(i, c)+gaussian(i, amp, lbd, sgm) for i in x]
	return y

def return_ew_voigt(number, fpars, a, b):
	label='v'+str(number)+'_'
	c=fpars['constant_c']
	lbd=fpars[label+'center']
	amp=fpars[label+'amplitude']
	sgm=fpars[label+'sigma']
	gma=fpars[label+'gamma']
	def integrand(x):
		return (1-(((c)+(amp*wofz((x-lbd + 1j*gma)/ (sgm*s2)).real / (sgm*s2pi)))/(c)))
	return [integrate.quad(integrand, a, b)[0], lbd]

def return_ew_gauss(number, fpars, a, b):
	label='g'+str(number)+'_'
	c=fpars['constant_c']
	lbd=fpars[label+'center']
	amp=fpars[label+'amplitude']
	sgm=fpars[label+'sigma']
	def integrand(x):
		return 1-(((c)+((amp/(s2pi*sgm)) * exp(-(1.0*x-lbd)**2 /(2*sgm**2))))/(c))
	return [integrate.quad(integrand, a, b)[0], lbd]

def getclick(event):
	xval = event.xdata
	global yset
	global limits
	global limflag
	global lineflag
	if limflag:
		limits.append(xval)
	if len(limits) == 1:
		yset=event.ydata
		print 'click on the right limit for fitting\npress enter when finished'
	if (len(limits) == 2):
		yset=(yset+event.ydata)/2
		limflag = False
		plt.gcf().canvas.mpl_disconnect(clcm)

def getkey(event):
	keyp = event.key
	xval = event.xdata
	global lambdas
	global linefitype
	global returnk
	if keyp == 'enter':
		returnk=True
	else:
		print '%8.2f%3s%1s' %(xval, ' : ', keyp)
		lambdas.append(xval)
		linefitype.append(keyp)

####################################### Interactive normalize function definition

def nclick(event):
	toolbar = plt.get_current_fig_manager().toolbar
	if event.button==1 and toolbar.mode=='':
		plt.plot(event.xdata,event.ydata,'rs',ms=8,picker=5,label='cont_pnt')
	plt.draw()

def npick(event):
	if event.mouseevent.button==3:
		if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
			event.artist.remove()

def ntype(event):
	global nmflag
	global gcont
	if event.key=='enter':
		cont_pnt_coord = []
		for artist in plt.gca().get_children():
			if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
				cont_pnt_coord.append(artist.get_data())
			elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
				artist.remove()
		cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
		sort_array = np.argsort(cont_pnt_coord[:,0])
		xnp,ynp = cont_pnt_coord[sort_array].T
		spline = splrep(xnp,ynp,k=3)
		continuum = splev(xo,spline)
		gcont=continuum
		plt.plot(xo,continuum,'r-',lw=2,label='continuum')
	elif event.key=='n':
		continuum = None
		for artist in plt.gca().get_children():
			if hasattr(artist,'get_label') and artist.get_label()=='continuum':
				continuum = artist.get_data()[1]
				break
		if continuum is not None:
			plt.cla()
			plt.plot(xo,normy/continuum,'k-',label='normalized')
		gcont=continuum
	elif event.key=='r':
		plt.cla()
		plt.plot(xo,normy, obstyle, linewidth=0.5, label=obsname)
	
	elif event.key=='w':
		for artist in plt.gca().get_children():
			if hasattr(artist,'get_label') and artist.get_label()=='normalized':
				data = np.array(artist.get_data())
				np.savetxt(obsname+'_norm.dat',data.T)
				print('saved to '+obsname+'_norm.dat')
	elif event.key=='x':
		print 'exiting normalizing mode. press enter on terminal'
		nmflag = False
		plt.gcf().canvas.mpl_disconnect(nclick)
		plt.gcf().canvas.mpl_disconnect(npick)
		plt.gcf().canvas.mpl_disconnect(nclick)
	elif event.key == 'a':
		print 'accepting fit and exiting normalizing mode. press enter on terminal'
		global yo
		yo=normy/gcont
		nmflag = False
		plt.gcf().canvas.mpl_disconnect(nclick)
		plt.gcf().canvas.mpl_disconnect(npick)
		plt.gcf().canvas.mpl_disconnect(nclick)
	plt.draw()

####################################### Spectra file reading

flglst=False
file_flag=False
while file_flag == False:
	obs=raw_input("dynamic spectra: ")
	file_flag=os.path.exists(obs)
multispec=' '
while not(is_number(multispec) or os.path.exists(multispec)):
	multispec=raw_input("how many static spectra [list filename]? ")
if (not is_number(multispec)):
	obstyle='k.'
	flglst=True
	modlst=multispec
	multispec=0
multispec=int(multispec)
A=[]
if flglst == False:
	obstyle='k.'
	if (multispec == 0):
		obstyle='k-'
	for e in range(multispec):
		promp="static spectra "+str(e+1)+": "
		file_flag=False
		while file_flag == False:
			stc=raw_input(promp)
			file_flag=os.path.exists(stc)
		A.append(stc)

#######################################	Line list or manual mode definition

linhas=[]
list_lamb=' '
while not(is_number(list_lamb) or os.path.exists(list_lamb)):
	list_lamb=raw_input("line list [or lambda]: ")
if is_number(list_lamb):
	linhas.append(float(list_lamb))
	lin='manual'
else:
	lin='lista'
	ln=open(list_lamb,"r")
	for cols in (raw.strip().split() for raw in ln):
		linhas.append(float(cols[0]))
	ln.close()
	
####################################### Reading dyn spectrum

if obs[-5:] != '.fits':
	o=open(obs,"r")
	xo=[]
	yo=[]
	for cols in (raw.strip().split() for raw in o):
		if is_number(cols[0]):
			xo.append(float(cols[0]))
			yo.append(float(cols[1]))
	o.close()
else:
	xo,yo=readfits(obs)
if obs.find('.') != -1:
	obsname=obs[:obs.find('.')]
else:
	obsname=obs

####################################### Reading static spectra

if flglst == True:
	mod_filenames=open(modlst, 'r')
	for cols in (raw.strip().split() for raw in mod_filenames):
		if len(cols) != 0:
			A.append(cols[0])
	multispec=len(A)
MX=[]
MY=[]
for e in range(multispec):
	if A[e][-5:] != '.fits':
		mod=open(A[e],"r")
		xm=[]
		ym=[]
		for cols in (raw.strip().split() for raw in mod):
			if is_number(cols[0]):
				xm.append(float(cols[0]))
				ym.append(float(cols[1]))
		mod.close()
	else:
		xm,ym=readfits(A[e])
	MX.append(xm)
	MY.append(ym)

####################################### Extra function files reading

xtrafiles=[f for f in os.listdir('.') if (os.path.isfile(f) and (f[-4:-1] == '.mp'))]
flaglines=False
if (len(xtrafiles) != 0):
	for f in xtrafiles:
		if f[-1] == 'l':
			flaglines=True
			linames=[]
			linlamb=[]
			linenames=open(f, 'r')
			for cols in (raw.strip().split() for raw in linenames):
				if (float(cols[1]) > xo[0]) and (float(cols[1]) < xo[-1]):
					linames.append(cols[0]+' '+cols[1])
					linlamb.append(float(cols[1]))

####################################### Plotting the graphs and features: main loop

plt.ion()
k=0
unit='Wavelength [$\AA$]'
win=1
range_error=0
legends=True
legmod=True
annotation=[]
while k != (len(linhas)):
	n=linhas[k]
	opc="go"
	while (opc != "n" and opc != "k"):
		flag = True
		flg_free = False
		while (flag == True):
			if legmod:
				if legends:
					plt.figure(figsize=(10,6)).subplots_adjust(right=0.72, left=0.062)
				else:
					plt.figure(figsize=(8,6)).subplots_adjust(right=0.96, left=0.08)
			legmod=False
			plt.xlabel(unit)
			plt.ylabel("Arbitrary Flux")
			if multispec != 0:
				plt.plot(xo,yo, obstyle, linewidth=0.5, label=obs)
				colors = plt.get_cmap('jet')(np.linspace(0, 1.0, multispec))
			else:
				plt.plot(xo,yo, obstyle, linewidth=0.5, label=obs)
			for e in range(multispec):
				plt.plot(MX[e],MY[e], linewidth=0.7, label=A[e], color=colors[e])
			i1=bissec(xo,n-win)
			i2=bissec(xo,n+win)
			try:
				lmi=min(yo[i1:i2])
				lma=max(yo[i1:i2])
			except ValueError:
				print 'lambda value out of range'
			for l in range(multispec):
				i1=bissec(MX[l],n-win)
				i2=bissec(MX[l],n+win)
				try:
					dl=min(MY[l][i1:i2])
					if dl < lmi:
						lmi=dl
					ul=max(MY[l][i1:i2])
					if ul > lma:
						lma=ul
				except ValueError:
					range_error += 1
			ysize=math.fabs(lma)-math.fabs(lmi)
			linleng=0
			if flaglines:
				linleng=1
				for i in range(len(linlamb)):
					xloc=bissec(xo,linlamb[i])
					loclim=yo[xloc]
					for l in range(multispec):
						i3=bissec(MX[l],linlamb[i])
						try :
							locd=min(MY[l][i3-1:i3+1])
							if locd < loclim:
								loclim=locd
						except :
							pass
					anot=plt.annotate(linames[i], xy=(linlamb[i], loclim), xytext=(linlamb[i], loclim-math.fabs(0.15*ysize)), arrowprops=dict(width=0.02, headwidth=0.05, shrink=0.2, frac=0.01), rotation='vertical', size=10)
				linwin=linames[(bissec(linlamb, n-win)-1):(bissec(linlamb, n+win)+1)]
			if len(annotation) !=0:
				for a in range(len(annotation)):
					plt.text(0.04, 0.03+0.06*(a+1), annotation[-a-1],horizontalalignment='left',verticalalignment='center', fontsize=16, transform = plt.gca().transAxes)
			if flg_free == False:
				plt.axis([n-win, n+win, lmi-math.fabs((0.05+(linleng*0.1))*ysize), lma+math.fabs(0.05*ysize)])
			else:
				plt.axis(curraxlim)
			if legends:
				plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., prop={'size':10})
			plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
			plt.draw()

####################################### Interactive dialogue and corresponding functions

			print "\nline: ",n
			if range_error != 0:
				print str(range_error)+' spectra out of range'
				range_error = 0
			print "     change dynamic spectra: c ; interactive normalizing: in"
			if (multispec == 1 and flglst != True):
				print "      change static spectra: h ; normalize: a"
			if (flglst == True):
				print "change static spectra group: h ; normalize: a"
			if lin != 'manual':
				print "       displace dyn spectra: l ; next line: n\n         divide dyn spectra: d ; save image: s"
			else:
				print "       displace dyn spectra: l ; new line: k\n         divide dyn spectra: d ; save image: s"
			if flg_free == False:
				print " free-look: f ; change unit: u ; modify window: w"
			else:
				print " free-look: f ; change unit: u"
			if legends:
				print "  annotate ; clear ; dyn_style ; trim"
			else:
				print "  annotate;  clear ; dyn_style ; untrim"
			print "  EW measuring: ew ; quit: q"
			opc=raw_input()
			if flg_free == True:
				curraxlim = plt.axis()
			if (opc == "h" and multispec == 1 and flglst != True):
				file_flag=False
				while file_flag == False:
					A[0]=raw_input("new static spectra filename:")
					file_flag=os.path.exists(A[0])
				mod=open(A[0],"r")
				if A[0][(len(A[0])-5):] != '.fits':
					xm=[]
					ym=[]
					for cols in (raw.strip().split() for raw in mod):
						if is_number(cols[0]):
							xm.append(float(cols[0]))
							ym.append(float(cols[1]))
					mod.close()
				else:
					xm,ym=readfits(A[0])
				MX[0]=xm
				MY[0]=ym
			elif opc == 'in':
				plt.cla()
				print '\nleft-click on the spectrum to draw points for the continuum fit\nright-click to remove the point.\ngive commands on plotting window'
				print '\noptions:\nenter: fit continuum ; n: normalize to continuum ; r: clear\na: accept fit and exit ; w: write to a file ; x: exit normalizing mode'
				normy=[m for m in yo]
				plt.plot(xo,normy, obstyle, linewidth=0.5, label=obsname)
				plt.draw()
				nmflag=True
				plt.gcf().canvas.mpl_connect('key_press_event',ntype)
				plt.gcf().canvas.mpl_connect('button_press_event',nclick)
				plt.gcf().canvas.mpl_connect('pick_event',npick)
				while nmflag:
					whatever=raw_input()
			elif opc == 'ew':
				limits=[]
				yset=0
				lambdas=[]
				linefitype=[]
				limflag=True
				lineflag=False
				print 'click on the left limit for fitting'
				clcm=plt.gcf().canvas.mpl_connect('button_press_event', getclick)
				whatever=raw_input()
				print '%7s%8.2f%3s%8.2f' %('range: ',limits[0],' - ',limits[1])
				print 'now mark the bottom of the lines for the fit.\n\'g\' for gaussian profile, \'v\' for voigt profile\npress enter when finished'
				kprm=plt.gcf().canvas.mpl_connect('key_press_event', getkey)
				whatever=raw_input()
				plt.gcf().canvas.mpl_disconnect(kprm)
				if linefitype.count('g')+linefitype.count('v') != len(linefitype):
					print 'you typed an unsupported profile. exiting ew module'
					break
				lineflag=False
				x_trim=xo[bissec(xo, limits[0]):bissec(xo, limits[1])]
				y_trim=yo[bissec(xo, limits[0]):bissec(xo, limits[1])]
				modf,pars = call_constant(x_trim, y_trim, yset)
				for i in range(len(lambdas)):
					if linefitype[i] == 'v':
						modf = modf+call_voigt(x_trim, y_trim, lambdas[i], i+1, pars)
					if linefitype[i] == 'g':
						modf = modf+call_gauss(x_trim, y_trim, lambdas[i], i+1, pars)
				out = modf.fit(y_trim, pars, x=x_trim)
				alim,blim=x_trim[0],x_trim[-1]
				for i in range(len(lambdas)):
					if linefitype[i] == 'v':
						print '%5s%2d%5s%6.4f%11s%8.2f' %('line ', i+1, ': EW = ', return_ew_voigt(i+1, out.best_values, lambdas[i]-0.5, lambdas[i]+0.5)[0], ', center = ',return_ew_voigt(i+1, out.best_values, alim, blim)[1])
					if linefitype[i] == 'g':
						print '%5s%2d%5s%6.4f%11s%8.2f' %('line ', i+1, ': EW = ', return_ew_gauss(i+1, out.best_values, lambdas[i]-0.5, lambdas[i]+0.5)[0], ', center = ',return_ew_gauss(i+1, out.best_values, alim, blim)[1])
				plt.plot(x_trim, [yset for i in range(len(x_trim))], 'k-')
				plt.plot(x_trim, out.best_fit, 'b-')
				if len(lambdas) > 1:
					colors = plt.get_cmap('copper')(np.linspace(0.1, 0.9, len(lambdas)))
					for e in range(len(lambdas)):
						if linefitype[e] == 'v':
							plt.plot(x_trim, build_voigt(e+1, out.best_values, x_trim), color=colors[e])
						if linefitype[e] == 'g':
							plt.plot(x_trim, build_gauss(e+1, out.best_values, x_trim), color=colors[e])
				plt.draw()
				whatever=raw_input('press enter to continue')
			elif opc == 'dyn_style':
				temp=obstyle
				try:
					obstyle=(raw_input('new ob style: '))
					plt.plot(xo,yo, obstyle)
				except:
					obstyle=temp
					print 'invalid style'
			elif opc == 'deton':
				cid=plt.gcf().canvas.mpl_connect('button_press_event', onclick)
			elif opc == 'detoff':
				plt.gcf().canvas.mpl_disconnect(cid)
			elif opc == 'trim':
				legends=False
				legmod=True
				plt.close()
			elif opc == 'untrim':
				legends=True
				legmod=True
				plt.close()
			elif opc == 'annotate':
				annotation.append(raw_input('text: '))
			elif opc == 'clear':
				annotation=[]
			elif (opc == "h" and flglst == True):
				A=[]
				file_flag=False
				while file_flag == False:
					modlst=raw_input('new list: ')
					file_flag=os.path.exists(modlst)
				mod_filenames=open(modlst, 'r')
				for cols in (raw.strip().split() for raw in mod_filenames):
					if len(cols) != 0:
						A.append(cols[0])
				multispec=len(A)
				MX=[]
				MY=[]
				for e in range(multispec):
					if A[e][(len(A[e])-5):] != '.fits':
						mod=open(A[e],"r")
						xm=[]
						ym=[]
						for cols in (raw.strip().split() for raw in mod):
							if is_number(cols[0]):
								xm.append(float(cols[0]))
								ym.append(float(cols[1]))
						mod.close()
					else:
						xm,ym=readfits(A[e])
					MX.append(xm)
					MY.append(ym)
			elif (opc == "a"):
				normlist=[]
				for e in range(multispec):
					factor=max(MY[e])/2.
					MY[e] = [d/factor for d in MY[e]]
				factor=max(yo)/2.
				yo = [d/factor for d in yo]
			elif (opc == "u"):
				unf=False
				while unf == False:
					ut=raw_input("A for angstroms; nm for nanometers; um for micrometers: ")
					if ut == 'A':
						unit='Wavelength [$\AA$]'
						unf = True
					elif ut == 'nm':
						unit='Wavelength [nm]'
						unf = True
					elif ut == 'um':
						unit='Wavelength [$\mu$m]'
						unf = True
					else:
						print 'enter a valid option'
			elif (opc == 'c'):
				file_flag=False
				while file_flag == False:
					obs=raw_input("new dynamic spectra filename:")
					file_flag=os.path.exists(obs)
				if obs[(len(obs)-5):] != '.fits':
				 	o=open(obs,"r")
				 	xo=[]
				 	yo=[]
				 	for cols in (raw.strip().split() for raw in o):
				 		if is_number(cols[0]):
				 			xo.append(float(cols[0]))
				 			yo.append(float(cols[1]))
				 	o.close()
				else:
				 	xo,yo=readfits(obs)
			elif (opc == "d"):
				div=input("factor: ")
				yo = [d/div for d in yo]
			elif (opc == "l"):
				des=input("displacement: ")
				xo = [d+des for d in xo]
			elif (opc == "s"):
				nom=raw_input("output filename: ")
				plt.savefig(nom, dpi=300)
			elif (opc == "k"):
				lbini=' '
				while not is_number(lbini):
					lbini=raw_input('lambda: ') 
				linhas.append(float(lbini))
				flag = False
				flg_free = False
			elif ((opc == "n") and (lin != 'manual')):
				flag = False
			elif (opc == 'w'):
				win=input('new window size: ')/2.
			elif (opc == "q"):
				exit()
			elif (opc == "f"):
				flg_free=True
				ylower=min(yo)-math.fabs(0.1*min(yo))
				yupper=max(yo)+math.fabs(0.1*max(yo))
				for e in range(multispec):
					if min(MY[e])-math.fabs(0.1*min(MY[e])) < ylower:
						ylower = min(MY[e])-math.fabs(0.1*min(MY[e]))
					if max(MY[e])+math.fabs(0.1*max(MY[e])) > yupper:
						yupper = max(MY[e])+math.fabs(0.1*max(MY[e]))
				curraxlim=[xo[0],xo[len(xo)-1],ylower,yupper]
			else:
				print "enter a valid option"
			if (opc != 'trim') and (opc != 'untrim'):
				plt.clf()
	k += 1
