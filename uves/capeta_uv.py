import subprocess
import pyfits
import os

def wave_window(file):
	spectra = pyfits.open(file)
	try:
		bin = spectra[0].header['CDELT1']
	except:
		bin = spectra[0].header['CD1_1']
	start = spectra[0].header['CRVAL1']
	points = len(spectra[0].data)
	end=bin*points+start
	return start,end

fil=open('names.dat', 'r')
specnames=[]
for cols in (raw.strip().split() for raw in fil):
	specnames.append(cols[0])
fil.close()
lbd=[]
spc=[]
ew=[]
chiex=[]
loggf=[]
ref=[]
c6=[]
err=[]
qf=[]
for s in specnames:
	lbd.append([])
	spc.append([])
	ew.append([])
	chiex.append([])
	loggf.append([])
	ref.append([])
	c6.append([])
	err.append([])
	qf.append([])
	
	start,end=wave_window(s+'_l.fits')
	segs=int(round((end-start)/100.))
	seglen=(end-start)/segs
	for i in range(segs):
		pr=subprocess.Popen(['./daospec'], shell=True, close_fds=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		print 'sh='+str(round(i*seglen+start, 2))+'\nlo='+str(round((i+1)*seglen+start, 2))+'\n\n'+s+'_l.fits'+'\n\n'
		if i == 0:
			pr.communicate('fw=13\nor=3\nsh='+str(4810)+'\nlo='+str(round((i+1)*seglen+start, 2))+'\n\n'+s+'_l.fits'+'\n\n')
			#pr.communicate('or=2\nsh='+str(round(5837., 2))+'\nlo='+str(round((i+1)*seglen+start, 2))+'\n\n'+s+'_L.fits'+'\n\n')
		elif i != segs-1:
			pr.communicate('fw=13\nor=3\nsh='+str(round(i*seglen+start, 2))+'\nlo='+str(round((i+1)*seglen+start, 2))+'\n\n'+s+'_l.fits'+'\n\n')
		else:
			pr.communicate('fw=13\nor=3\nsh='+str(round(i*seglen+start, 2))+'\nlo='+str(5750.)+'\n\n'+s+'_l.fits'+'\n\n')
			#pr.communicate('or=2\nsh='+str(round(i*seglen+start, 2))+'\nlo='+str(round(6816.))+'\n\n'+s+'_L.fits'+'\n\n')
		results=open(s+'_l.daospec', 'r')
		for cols in (raw.strip().split() for raw in results):
			if (cols[0][0] != 'R') and (len(cols) > 5):
				lbd[-1].append(float(cols[5]))
				spc[-1].append(cols[6])
				ew[-1].append(float(cols[2]))
				chiex[-1].append(float(cols[7]))
				loggf[-1].append(float(cols[9]))
				ref[-1].append(cols[8])
				c6[-1].append(cols[10])
				err[-1].append(float(cols[3]))
				qf[-1].append(float(cols[4]))
		results.close()
		os.remove(s+'_l.daospec')
	
	start,end=wave_window(s+'_u.fits')
	segs=int((end-start)/80.)
	seglen=(end-start)/segs
	for i in range(segs):
		pr=subprocess.Popen(['./daospec'], shell=True, close_fds=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		print 'sh='+str(round(i*seglen+start, 2))+'\nlo='+str(round((i+1)*seglen+start, 2))+'\n\n'+s+'_u.fits'+'\n\n'
		if i == 0:
			pr.communicate('fw=16\nor=3\nsh='+str(round(5850., 2))+'\nlo='+str(round((i+1)*seglen+start, 2))+'\n\n'+s+'_u.fits'+'\n\n')
		elif i != segs-1:
			pr.communicate('fw=16\nor=3\nsh='+str(round(i*seglen+start, 2))+'\nlo='+str(round((i+1)*seglen+start, 2))+'\n\n'+s+'_u.fits'+'\n\n')
		else:
			pr.communicate('fw=16\nor=3\nsh='+str(round(i*seglen+start, 2))+'\nlo='+str(round(6780.))+'\n\n'+s+'_u.fits'+'\n\n')
		results=open(s+'_u.daospec', 'r')
		for cols in (raw.strip().split() for raw in results):
			if (cols[0][0] != 'R') and (len(cols) > 5):
				lbd[-1].append(float(cols[5]))
				spc[-1].append(cols[6])
				ew[-1].append(float(cols[2]))
				chiex[-1].append(float(cols[7]))
				loggf[-1].append(float(cols[9]))
				ref[-1].append(cols[8])
				c6[-1].append(cols[10])
				err[-1].append(float(cols[3]))
				qf[-1].append(float(cols[4]))
		results.close()
		os.remove(s+'_u.daospec')
	
	st = s.split('/')[-1]
	outfil=open('raies_'+st+'.dat', 'w')
	outfil.write('####\n')
	outfil.write('FORMAT(F7.2,1X,A2,I1,F5.2,I3,F6.2,3X,e9.3,3X,F6.1)\n')
	outfil.write('lambda   el kiex ref log gf      C6        EW    err     QF\n')
	for j in range(len(lbd[-1])):
		if spc[-1][j] == 'FE2': 
			outfil.write('%7.2f%4s%6.2f%3s%6.2f%12s%9.1f%6.1f%8.3f\n' %(lbd[-1][j], spc[-1][j], chiex[-1][j], ref[-1][j], loggf[-1][j], c6[-1][j], ew[-1][j], err[-1][j], qf[-1][j]))
	for j in range(len(lbd[-1])):
		if spc[-1][j] == 'FE1': 
			outfil.write('%7.2f%4s%6.2f%3s%6.2f%12s%9.1f%6.1f%8.3f\n' %(lbd[-1][j], spc[-1][j], chiex[-1][j], ref[-1][j], loggf[-1][j], c6[-1][j], ew[-1][j], err[-1][j], qf[-1][j]))
	outfil.close()
