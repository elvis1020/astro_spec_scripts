import subprocess
import pyfits
import os

setups = ['HR11', 'HR13', 'HR14']
rangei = [5594, 6115, 6301]
rangef = [5834, 6399, 6690]
subdir = '6522_giraffe'

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
    
    for n, g in enumerate(setups):
    	pr=subprocess.Popen(['./daospec'], shell=True, close_fds=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    	pr.communicate('fw=10.0\nor=3\nsh='+str(rangei[n])+'\nlo='+str(rangef[n])+'\n\n'+'6522_giraffe/'+g+'/'+s+'_'+g+'.fits'+'\n\n')
    	with open(subdir+'/'+g+'/'+s+'_'+g+'.daospec', 'r') as results:
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
		#os.remove('6522_giraffe/'+g+'/'+s+'_HR11.daospec')
    
    #pr=subprocess.Popen(['./daospec'], shell=True, close_fds=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #pr.communicate('fw=10.0\nor=3\nsh='+str(5595.5)+'\nlo='+str(5836.)+'\n\n'+'6522_giraffe/HR11/'+s+'_HR11.fits'+'\n\n')
    #results=open('6522_giraffe/HR11/'+s+'_HR11.daospec', 'r')
    #for cols in (raw.strip().split() for raw in results):
    #    if (cols[0][0] != 'R') and (len(cols) > 5):
    #        lbd[-1].append(float(cols[5]))
    #        spc[-1].append(cols[6])
    #        ew[-1].append(float(cols[2]))
    #        chiex[-1].append(float(cols[7]))
    #        loggf[-1].append(float(cols[9]))
    #        ref[-1].append(cols[8])
    #        c6[-1].append(cols[10])
    #        err[-1].append(float(cols[3]))
    #        qf[-1].append(float(cols[4]))
    #results.close()
    ##os.remove('6522_giraffe/HR11/'+s+'_HR11.daospec')
    #           
    #pr=subprocess.Popen(['./daospec'], shell=True, close_fds=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #pr.communicate('fw=10.0\nor=3\nsh='+str(5836.06)+'\nlo='+str(6140.83)+'\n\n'+'6522_giraffe/HR12/'+s+'_HR12.fits'+'\n\n')
    #results=open('6522_giraffe/HR12/'+s+'_HR12.daospec', 'r')
    #for cols in (raw.strip().split() for raw in results):
    #    if (cols[0][0] != 'R') and (len(cols) > 5):
    #        lbd[-1].append(float(cols[5]))
    #        spc[-1].append(cols[6])
    #        ew[-1].append(float(cols[2]))
    #        chiex[-1].append(float(cols[7]))
    #        loggf[-1].append(float(cols[9]))
    #        ref[-1].append(cols[8])
    #        c6[-1].append(cols[10])
    #        err[-1].append(float(cols[3]))
    #        qf[-1].append(float(cols[4]))
    #results.close()
    ##os.remove('6522_giraffe/HR12/'+s+'_HR12.daospec')
    
    outfil=open('raies_'+s+'.dat', 'w')
    outfil.write(s+'\n')
    outfil.write('FORMAT(F7.2,1X,A2,I1,F5.2,I3,F6.2,3X,e9.3,3X,F6.1)\n')
    outfil.write('lambda   el kiex ref log gf      C6        EW    err     QF\n')
    for j in range(len(lbd[-1])):
        if spc[-1][j] == 'FE2': 
            outfil.write('%7.2f%4s%6.2f%3s%6.2f%12s%9.1f%6.1f%8.3f\n' %(lbd[-1][j], spc[-1][j], chiex[-1][j], ref[-1][j], loggf[-1][j], c6[-1][j], ew[-1][j], err[-1][j], qf[-1][j]))
    for j in range(len(lbd[-1])):
        if spc[-1][j] == 'FE1': 
            outfil.write('%7.2f%4s%6.2f%3s%6.2f%12s%9.1f%6.1f%8.3f\n' %(lbd[-1][j], spc[-1][j], chiex[-1][j], ref[-1][j], loggf[-1][j], c6[-1][j], ew[-1][j], err[-1][j], qf[-1][j]))
    outfil.close()
