###########################################################                         
#UVES fibre mode spectra treatment, written by E. Cantelli#
###########################################################
import pyfits
import sys
from pyraf import iraf

#Is this OB using SimCal? [y/n] 
simcal = 'n'

#Enter list filenames for all spectra, separated by chips and templates
#Choose which chip will be used to calculate the RV's on line 148 : 'l' or 'u'
lista_lower = 'names_l.dat'
lista_upper = 'names_u.dat'
tpl_l = 'sun_580_lower.fits'
tpl_u = 'sun_580_upper.fits'

#Enter infotab file
infotab = 'infotab.fits'

#OR

#Enter used fiber positions starting by zero and respective object names
#fibnum = ['0', '1', '2', '3', '4', '5', '6', '7']
#ob_id = ['8326', '2', '4152', '8025', '8734', '6159', '8233', '8731']

#Do NOT FORGET to comment line 123 if you use this 'manual' mode.

def fibnamer(infotab):
	pack = pyfits.open(infotab)
	table = pack[1].data
	fib = []
	obj = []
	for i in range(len(table)):
		if (table[i][2] != 0) and (table[i][3][0:4] != 'Grid') and (table[i][3] != ''):
			fib.append(str(table[i][2] - simcal))
                        obj.append(str(table[i][3]))
        return fib, obj

def clean_cr(lista):
	inp = '@' + lista
	iraf.images()
	iraf.imfit()
	iraf.lineclean(inp, inp, interac = 'no', low_rej = '5.6', high_rej = '2.5')

def run_fxcor(lstl, lstu, tpl, tpu, chip):
	if chip == 'l':
		templ = tpl
		lst = lstl
	else:
		templ = tpu
		lst = lstu
	lista = '@' + lst
	chip_par = 'rvout_' + chip
	iraf.noao()
	iraf.rv()
	iraf.fxcor(lista, templ, output=chip_par, verbose = 'txtonly', interac = 'no', observa = 'paranal')
	rv = []
	rvlist = open(chip_par + '.txt', "r")
	for cols in (raw.strip().split() for raw in rvlist):
		vobs_index = len(cols) - 3
		if cols[0][0] != '#':
			rv.append(float(cols[vobs_index]))
	rvlist.close()
	return rv

def include_rv(namefil, rv):
	namelist = []
	names = open(namefil, "r")
	for cols in (raw.strip().split() for raw in names):
		namelist.append(cols[0])
	names.close()
	for i in range(len(namelist)):
		spectra = pyfits.open(namelist[i], mode = 'update')
		prihdr = spectra[0].header
		prihdr.set('RV', rv[i])
		spectra.close()

def run_dopcor(lista):
	inp='@' + lista
	iraf.noao()
	iraf.onedspec()
	iraf.dopcor(inp, inp, 'RV', isveloc = 'yes', apertur = '*')

def run_scombine(listin, fn):
	namefil = []
	namelist = open(listin, "r")
	for cols in (raw.strip().split() for raw in namelist):
		namefil.append(cols[0])
	namelist.close()
	obs = []
	count = -1
	for i in range(len(namefil)):
		if namefil[i - 1][0:2] != namefil[i][0:2]:
			obs.append([])
			count += 1
		obs[count].append(namefil[i])
	spec = []
	for i in range(len(obs[0])):
		spec.append([])
		for j in range(len(obs)):
			spec[i].append(obs[j][i])
	iraf.noao()
	iraf.onedspec()
	for i in range(len(spec)):
		temp = open('temp_list', 'w')
		for j in range(len(spec[i])):
			temp.write(spec[i][j] + '\n')
		temp.close()
		if (spec[i][0][-6] == 'L') or (spec[i][0][-6] == 'U'):
			name_index = '0_' + spec[i][0][-6]
		else:
			name_index = spec[i][0][-6] + '_' + spec[i][0][-8]
		if len(fibnum) != 0:
			newname = ob_id[fibnum.index(name_index[0])] + '_' + name_index[-1]
			name_index = newname
		iraf.scombine('@temp_list', name_index, logfile = 'combine_log')
		fn.append(name_index)

def normalize(fn):
	iraf.noao()
	iraf.onedspec()
	for i in range(len(fn)):
		if fn[i][-1] == 'L':
			nord = '10'
		else:
			nord = '7'
		iraf.continuum(fn[i], fn[i], order=nord, ask='no', logfile='norm_log')

def gauss2(fn):
	iraf.images()
	iraf.imfilter()
	for i in (fn):	
		iraf.gauss(i, i, '2')

if simcal == 'n':
	simcal = 1
else:
	simcal = 0
fibnum, ob_id = fibnamer(infotab)
print ('Fibre association made:')
print fibnum
print ob_id
print('----------Cleaning CR\'s----------')
clean_cr(lista_lower)
clean_cr(lista_upper)
print('----------Measuring RV\'s----------')
rv_file = run_fxcor(lista_lower, lista_upper, tpl_l, tpl_u, 'u')
include_rv(lista_lower, rv_file)
include_rv(lista_upper, rv_file)
print('----------Performing Doppler correction----------')
run_dopcor(lista_lower)
run_dopcor(lista_upper)
out_names = []
print('----------Combining spectra----------')
run_scombine(lista_lower, out_names)
run_scombine(lista_upper, out_names)
print('----------Normalizing spectra----------')
normalize(out_names)
print('----------Smoothing spectra----------')
gauss2(out_names)
print('Done')
