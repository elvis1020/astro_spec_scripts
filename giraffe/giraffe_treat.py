###############################################################
#GIRAFFE Medusa mode spectra treatment, written by E. Cantelli#
###############################################################
import pyfits
import os
import sys
from pyraf import iraf

lista='pack_names'
gi_setups = ['H379.0', 'H395.8', 'H412.4', 'H429.7', 'H447.1', 'H465.6', 'H484.5', 'H504.8', 'H525.8',
            'H548.8', 'H572.8', 'H599.3', 'H627.3', 'H651.5', 'H665.0', 'H679.7', 'H710.5', 'H737.0']
gi_templates = ['sun_HR1_379.0.fits', 'sun_HR2_395.8.fits', 'sun_HR3_412.4.fits', 'sun_HR4_429.7.fits',
               'sun_HR5_447.1.fits', 'sun_HR6_465.6.fits', 'sun_HR7_484.5.fits', 'sun_HR8_504.8.fits',
               'sun_HR9_525.8.fits', 'sun_HR10_548.8.fits', 'sun_HR11_572.8.fits', 'sun_HR12_599.3.fits',
               'sun_HR13_627.3.fits', 'sun_HR14_651.5A.fits', 'sun_HR15_665.0.fits', 'sun_HR15_679.7.fits',
               'sun_HR16_710.5.fits', 'sun_HR17_737.0.fits']

def clean_cr(lista):
	inp=lista
	iraf.images()
	iraf.imfit()
	iraf.lineclean(inp, inp, interac='no', low_rej='4.5', high_rej='2.5')

def extract(packnam):
	obs=packnam[-6]
	pack=pyfits.open(packnam)
	table=pack[1].data
	head=pack[0].header
	setup=head['HIERARCH ESO INS EXP MODE']
	fib=[]
	obj=[]
	valid_obj=[]
	for i in range(len(table)):
		fib.append(table[i][0])
		obj.append(table[i][7])
	iraf.noao()
	iraf.onedspec()
	os.makedirs(setup+'/'+obs)
	for k in range(len(fib)):
		if (obj[k] != 'CALSIM') and (obj[k][0:4] != 'Grid'):
			valid_obj.append(setup+'/'+obs+'/'+obj[k])
			iraf.scopy(packnam+'['+str(fib[k])+',*]', setup+'/'+obs+'/'+obj[k])
	return valid_obj, setup

def run_fxcor(st_names, templ, setup):
	arq=open(setup+'/'+'speclist', 'w')
	for i in range(len(st_names)):
		arq.write(st_names[i]+'.fits\n')
	arq.close()
	lista='@'+setup+'/'+'speclist'
	clean_cr(lista)
	outname=setup+'/'+'rvout_'+setup
	iraf.noao()
	iraf.rv()
	iraf.fxcor(lista, templ, output=outname, verbose='txtonly', interac='no', observa='paranal')
	rv=[]
	rvlist=open(outname+'.txt', "r")
	for cols in (raw.strip().split() for raw in rvlist):
		vobs_index=len(cols)-3
		if cols[0][0] != '#':
			rv.append(float(cols[vobs_index]))
	rvlist.close()
	return rv

def include_rv(setup, rv):
	namefil=setup+'/speclist'
	namelist=[]
	names=open(namefil,"r")
	for cols in (raw.strip().split() for raw in names):
		namelist.append(cols[0])
	names.close()
	for i in range(len(namelist)):
		spectra = pyfits.open(namelist[i], mode='update')
		prihdr = spectra[0].header
		prihdr.set('RV', rv[i])
		spectra.close()
		
def nm_to_A(setup, filelist):
	for i in filelist:
		namefil=setup+'/'+i
		spec=pyfits.open(namefil, mode='update')
		prihdr=spec[0].header
		prihdr['CRVAL1']=prihdr['CRVAL1']*10
		prihdr['CDELT1']=prihdr['CDELT1']*10
		try:
			prihdr['CD1_1']=prihdr['CD1_1']*10
		except:
			pass
		prihdr['CUNIT2']='angstroms'
		spec.close()

def run_dopcor(lista):
	inp='@'+lista
	iraf.noao()
	iraf.onedspec()
	iraf.dopcor(inp, inp, 'RV', isveloc='yes', apertur='*')

def run_scombine(listin, fn, setup):
	namefil=[]
	namelist=open(listin,"r")
	for cols in (raw.strip().split() for raw in namelist):
		namefil.append(cols[0])
	namelist.close()
	obs=[]
	count=-1
	for i in range(len(namefil)):
		if namefil[i-1][namefil[i-1].index('/')+1:][0] != namefil[i][namefil[i].index('/')+1:][0]:
			obs.append([])
			count += 1 
		obs[count].append(namefil[i])
	spec=[]
	for i in range(len(obs[0])):
		spec.append([])
		for j in range(len(obs)):
			Flag=True
			try:
				val=obs[j][i]
				Flag=True
			except IndexError:
				Flag=False
			if Flag:
				spec[i].append(obs[j][i])
		if len(spec[i]) != len(obs):
			spec.pop(i)
	iraf.noao()
	iraf.onedspec()
	for i in range(len(spec)):
		temp=open(setup+'/temp_list', 'w')
		for j in range(len(spec[i])):
			temp.write(spec[i][j]+'\n')
		temp.close()
		longn=spec[i][0]
		specname=longn[(longn.index('/')+1):][(longn[(longn.index('/')+1):]).index('/')+1:]
		iraf.scombine('@'+setup+'/temp_list', setup+'/'+specname, logfile=setup+'/combine_log')
		fn.append(specname)

def normalize(setup, filelist):
	iraf.noao()
	iraf.onedspec()
	for i in (filelist):
		name=setup+'/'+i
		iraf.continuum(name, name, order=1, ask='no', logfile=setup+'/norm_log')

def gauss2(setup, filelist):
	iraf.images()
	iraf.imfilter()
	for i in (filelist):
		name=setup+'/'+i
		iraf.gauss(name, name, '2')

packnames=[]
packlist=open(lista, 'r')
for cols in (raw.strip().split() for raw in packlist):
	packnames.append(cols[0])
file_names=[]
setups=[]
stp_ind=0
print '----------Extracting packs----------'
for i in range(len(packnames)):
	fn, st = extract(packnames[i])
	if stp_ind != st:
		file_names.append(fn)
		setups.append(st)
	else:
		for k in fn:
			file_names[-1].append(k)
	stp_ind=st
rv=[]
for i in range(len(setups)):
	print '----------Cleaning CR\'s and measuring RV\'s from setup '+setups[i]+'----------'
	st_rv=run_fxcor(file_names[i], gi_templates[gi_setups.index(setups[i])], setups[i])
	rv.append(st_rv)
filens=[]
for i in range(len(setups)):
	include_rv(setups[i], rv[i])
	print '----------Doppler correcting setup '+setups[i]+'----------'
	run_dopcor(setups[i]+'/speclist')
	filens.append([])
	print '----------Combining setup '+setups[i]+'----------'
	run_scombine(setups[i]+'/speclist', filens[i], setups[i])
	print '----------Converting setup '+setups[i]+' to angstroms----------'
	nm_to_A(setups[i], filens[i])
	print '----------Normalizing setup '+setups[i]+'----------'
	normalize(setups[i], filens[i])
	print '----------Smoothing setup '+setups[i]+'----------'
	gauss2(setups[i], filens[i])
	print 'Setup '+setups[i]+' completed successfully\n'
print 'Done'
