#!/usr/bin/python
#ESO ADP extracting task for GIRAFFE
import multiprocessing as mp
import numpy as np
import os
import pickle

from tasks import *

global prog
global xdir
global tr_mag
global tp
global hkeys

rdir = '.'
xdir = 'extracted'
hkeys = ['PROG_ID',
         'OBJECT',
         'RA',
         'DEC',
         'HIERARCH ESO INS FILT NAME',
         'HIERARCH ESO INS SLIT NAME',
         'HELICORR',
         'DATE-OBS',
         'HIERARCH ESO OBS NAME',
         'EXPTIME',
         'HIERARCH ESO TEL AIRM START',
         'HIERARCH ESO TEL AIRM END',
         'HIERARCH ESO TEL AMBI FWHM START',
         'HIERARCH ESO TEL AMBI FWHM END',
         'SNR',
         'SPEC_RES']

templates = [['HR11', 'sunh11.dat'],
             ['HR12', 'sunh12.dat'],
             ['HR13', 'sunh13.dat'],
             ['HR14', 'sunh14.dat'],
             ['HR15', 'sunh15.dat'],
             ['HR18', 'sunh18.dat']]

tp = {}
for i in templates:
    tp[i[0]] = np.genfromtxt(i[1], unpack=True)

filstruc, dsetid, prog_mags = structurate(rdir, xdir)

with open('fistruc.dp', 'w') as filstrucdp:
	pickle.dump(filstruc, filstrucdp)
with open('dsetid.dp', 'w') as dsetiddp:
	pickle.dump(dsetid, dsetiddp)
with open('prog_mags.dp', 'w') as prog_magsdp:
	pickle.dump(prog_mags, prog_magsdp)

print '\nValidated objects:'
for n in filstruc.keys():
    print n + ': '+str(len(filstruc[n]))+' files'

pool = mp.Pool(processes=12)

combine = {}
for p in filstruc.keys():
    combine[p] = {}


for prog in filstruc.keys():
    dsetid[prog] = []
    tr_mag = [list(l) for l in list(np.transpose(prog_mags[prog]))]
    print 'Writing spectra to directory '+xdir+'/'+prog
    
    xtractions = [pool.apply_async(filxtract, args=(f, prog, hkeys, dsetid,
                 tr_mag, xdir, tp)) for f in filstruc[prog]]
    obpars = [p.get() for p in xtractions]
    for i in obpars:
        if not i[0] in combine[prog]:
            combine[prog][i[0]] = {}
        if not i[1] in combine[prog][i[0]]:
            combine[prog][i[0]][i[1]] = []
        combine[prog][i[0]][i[1]].append(i[2])

print 'Extracting done'

with open('combine.dp', 'w') as combinedp:
	pickle.dump(combine, combinedp)

combs = []
for prog in combine:
    print 'Combining program '+prog 
    for setup in combine[prog]:
    	try:
        	direc = 'extracted/combined/'+prog+'/'+setup
        	try:
        	    os.makedirs(direc)
        	except:
        	    pass
        	for obj in combine[prog][setup]:
        	    cbname = direc+'/'+obj
        	    combs.append(pool.apply_async(combiner,
        	                 args=(combine[prog][setup][obj], cbname)))
        except:
            print("error on setup "+setup)
            pass
info = [p.get() for p in combs]

vrtab = open('extracted/combined/vrtab.dat', 'w')
vrtab.write('%10s  %7s\n'%('Object', 'VHelio'))
for i in info:
	vrtab.write('%10s  %5.2f\n'%(i[0].split('/')[-1], i[1]))
vrtab.close()

print 'Combining done'
