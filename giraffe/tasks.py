import numpy as np
from PyAstronomy import pyasl
from lmfit.models import GaussianModel, ConstantModel, SkewedGaussianModel,\
                         PolynomialModel
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter
import pyfits as pf
from astropy.io import fits
import re
import os
import matplotlib.pyplot as plt
import traceback
import sys
import warnings

__all__ = ["deg2HMS", "fit", "clean", "pyfxcor", "combiner", "validate",
           "incorporate", "structurate", "filxtract", "fits_save",
           "measure_SNR", "cnt_regions"]


def fits_save(spec, filen, stup, obnm, vh, n, snr):
	
	filtspec = gaussian_filter(np.copy(spec[1]), 0.8).astype('float32')
	hdu = fits.PrimaryHDU(filtspec)
	delt = np.mean([spec[0][i+1]-spec[0][i] for i in range(len(spec[0])-1)])
	hdu.header['CTYPE1'] = 'LINEAR'
	hdu.header['CRPIX1'] = 1.
	hdu.header['CRVAL1'] = spec[0][0]
	hdu.header['CDELT1'] = delt
	hdu.header['CD1_1'] = delt
	hdu.header['OBJ_NAME'] = obnm
	hdu.header['SETUP'] = stup
	hdu.header['VHELIO'] = vh
	hdu.header['N_SPEC'] = n
	hdu.header['SNR'] = snr
	name = filen+stup+'.fits'
	try:
		hdu.writeto(name)
	except:
		os.remove(name)
		hdu.writeto(name)

def measure_SNR(spec, region):
    
    idx1 = (np.abs(spec[0]-region[0])).argmin()
    idx2 = (np.abs(spec[0]-region[1])).argmin()
    reg = np.array(spec[1][idx1:idx2])
    if (len(reg) != 0) and (np.std(reg) > 1e-10):
        return np.mean(reg)/np.std(reg)
    else:
        return 0


cnt_regions = [(5672.5, 5674.8),
               (5763.9, 5765.8),
               (6069.5, 6075.8),
               (6184.3, 6184.9),
               (6328.0, 6329.5),
               (6371.9, 6372.5),
               (6422.4, 6423.6),
               (6426.6, 6427.8),
               (6577.8, 6578.4),
               (6615.1, 6616.4),
               (6626.4, 6627.1),
               (6691.15, 6691.95),
               (6718.1, 6719.15),
               (6778.1, 6779.3),
               (7496.5, 7497.0),
               (7519.9, 7520.3)]


def deg2HMS(ra='', dec='', round=False):
    
    RA, DEC, rs, ds = '', '', '', ''
    if dec:
        if str(dec)[0] == '-':
            ds, dec = '-', abs(dec)
        deg = int(dec)
        decM = abs(int((dec-deg)*60))
        if round:
            decS = int((abs((dec-deg)*60)-decM)*60)
        else:
            decS = (abs((dec-deg)*60)-decM)*60
        DEC = '{0}{1}:{2}:{3}'.format(ds, deg, decM, decS)
    if ra:
        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra/15)
        raM = int(((ra/15)-raH)*60)
        if round:
            raS = int(((((ra/15)-raH)*60)-raM)*60)
        else:
            raS = ((((ra/15)-raH)*60)-raM)*60
        RA = '{0}{1}:{2}:{3}'.format(rs, raH, raM, raS)
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC


def fit(spectra, obj, sigma=2.0, ord=4, iter=4):
    
    poly = PolynomialModel(3)
    pars = poly.make_params()
    for p in range(4):
        label = 'c'+str(p)
        pars[label].set(value=1., vary=True)
    wkcopy = np.copy(spectra[1])
    truesp = [i for i in wkcopy if i > 5]
    truex = [spectra[0][i] for i in range(len(spectra[1])) if spectra[1][i] > 5]
    outcont = poly.fit(truesp, pars, x=truex)
    firstcont = outcont.eval(x=spectra[0])
    
    xn = np.copy(spectra[0])
    yn = np.copy(spectra[1])/firstcont
    
    pl1=plt.subplot((iter+1)*100+11)
    pl1.plot(xn, spectra[1], 'k-', linewidth=0.3)
    pl1.plot(xn, firstcont, 'r-', linewidth=0.6)
    pl1.set_ylim([0, np.mean(firstcont)*1.5])
    
    for i in range(iter):
        i_=np.copy(i)
        niter=str(i_+1)
        sigma = sigma-i*0.21*sigma
        
        md = np.median(yn)
        n = len([i for i in yn if i > 0.1])
        offset = (len(xn)-n)/2
        absor = md - min(yn[offset:n-offset])
        freq, bin = np.histogram(yn, bins=50, range=(md-absor, md+absor))
        rebin = [(bin[b+1]+bin[b])/2 for b in range(len(bin)-1)]
        
        
        gauss = SkewedGaussianModel()
        pars = gauss.make_params()
        pars['center'].set(value=md, vary=True)
        pars['amplitude'].set(vary=True)
        pars['sigma'].set(vary=True)
        pars['gamma'].set(vary=True)
        out = gauss.fit(freq, pars, x=rebin)
        
        var = sigma*out.best_values['sigma']
        xrbn = np.linspace(rebin[0], rebin[-1], num=100)
        yrbn = list(out.eval(x=xrbn))
        mode = xrbn[yrbn.index(max(yrbn))]
        
        ync = np.copy(spectra[1])
        xnc = np.copy(spectra[0])
        
        mask = []
        for j in range(len(yn)):
            if (yn[j] > mode+var/2) or (yn[j] < mode-var/2):
                mask.append(False)
            else:
                mask.append(True)
        mask = np.array(mask)
        ync = ync[mask]
        xnc = xnc[mask]
        
        poly2 = PolynomialModel(ord)
        pars2 = poly2.make_params()
        for p in range(ord+1):
            label = 'c'+str(p)
            pars2[label].set(value=1., vary=True)
        outcont2 = poly2.fit(ync, pars2, x=xnc)
        
        contf = outcont2.eval(x=xn)
        yn = spectra[1]/contf
        err = spectra[2]/contf
        
        pln=plt.subplot(int((iter+1)*100+10+(i_+2)))
        pln.plot(xn, yn*(np.mean(contf)*0.8), 'k-', linewidth=0.3)
        pln.plot(xnc, ync, 'r-', linewidth=0.3)
        pln.plot(xn, contf, 'b-', linewidth=0.6)
        pln.set_ylim([0, np.mean(contf)*1.2])
        
    plt.savefig(obj[0]+'_fit.png', dpi=300)
    plt.clf()
        
    return np.array([xn, yn, err])


def clean(spectra, sigma=2.6):
    
    md = np.median(spectra[1])
    n = int(len(spectra[0])*0.8)
    offset = (len(spectra[0])-n)/2
    absor = md - min(spectra[1][offset:n-offset])
    freq, bin = np.histogram(spectra[1], bins=50, range=(md-absor, md+absor))
    rebin = [(bin[b+1]+bin[b])/2 for b in range(len(bin)-1)]
    
    gauss = SkewedGaussianModel()
    pars = gauss.make_params()
    pars['center'].set(value=md, vary=True)
    pars['amplitude'].set(vary=True)
    pars['sigma'].set(vary=True)
    pars['gamma'].set(vary=True)
    out = gauss.fit(freq, pars, x=rebin)
    
    var = sigma*out.best_values['sigma']
    xrbn = np.linspace(rebin[0], rebin[-1], num=100)
    yrbn = list(out.eval(x=xrbn))
    mode = xrbn[yrbn.index(max(yrbn))]
    
    xn = np.copy(spectra[0])
    yn = np.copy(spectra[1])
    
    ist=0
    pts=[]
    errflg = False
    for i in range(len(xn)):
        if (yn[i] > mode+var) and (errflg == False):
            cnt = 0
            errflg = True
            ist = np.copy(i)
            cnt += 1
        if (yn[i] > mode+var) and (errflg == True):
            cnt += 1
        if (yn[i] < mode+var) and (errflg == True):
            pts = np.linspace(yn[ist-1], yn[i], cnt+2)[1:-1]
            for p in range(ist, i):
                yn[p] = pts[p-ist]
            errflg = False
            
    return np.array([xn, yn])


def pyfxcor(inspec, template, obj, vmin=-400., vmax=400., res=3, rej=200):
    
    rv, cc = pyasl.crosscorrRV(inspec[0], inspec[1], template[0], template[1],
                               vmin, vmax, res, skipedge=rej)
    
    cen_gs = np.argmax(cc)
    perfx, perfy = rv[cen_gs-5:cen_gs+6], cc[cen_gs-5:cen_gs+6]
    
    try:
        gauss = ConstantModel() + GaussianModel()
        pars = gauss.make_params()
        pars['center'].set(value=rv[np.argmax(cc)], vary=True)
        pars['amplitude'].set(value=max(cc), vary=True)
        pars['sigma'].set(vary=True)
        pars['c'].set(value=0, vary=True)
        out = gauss.fit(perfy, pars, x=perfx)
        ct = out.best_values['center']
        cterr = out.params['center'].stderr
    except:
        plt.plot(inspec[0], inspec[1])
        plt.savefig(obj+'.png', dpi=300)
        pl.clf()
        return 'error', ''
    
    plt.subplot(311)
    plt.plot(rv, cc)
    curraxlim = plt.axis()
    out.plot_fit(numpoints=100)
    plt.axis(curraxlim)
    plt.subplot(312)
    plt.plot(obj[1][0], obj[1][1], 'r-', linewidth=0.5)
    plt.subplot(313)
    plt.plot(inspec[0], inspec[1], 'b-', linewidth=0.5)
    plt.savefig(obj[0]+'.png', dpi=300)
    plt.clf()
    
    return ct, cterr


def combiner(specs, out):
    
    try:
    	sp = []
    	vrs = {'rel': [], 'corr': [], 'run': []}
    	plt.subplot(311)
    	for s in specs:
    	    try:
    	    	spectra = np.genfromtxt(s, comments='#', unpack=True)
    	    	hd1 = open(s)
    	    	run = s.split('/')[-2]
    	    	hd = hd1.readlines()[1].split()
    	    	vrs['rel'].append(float(hd[-1]))
    	    	vrs['corr'].append(float(hd[7]))
    	    	vrs['run'].append(run)
    	    	hd1.close()
    	    	sp.append(spectra)
    	    	plt.plot(spectra[0], spectra[1], linewidth=0.3)
    	    except:
    	    	print("problem with file "+s)
    	    	pass
    	
    	N = np.ceil(np.mean([len(s[0]) for s in sp]))
    	lmin = max([s[0][0] for s in sp])
    	lmax = min([s[0][-1] for s in sp])
    	rblbd = np.linspace(lmin, lmax, N)
    	rbsp = []
    	errors = np.sum(np.asarray(sp), axis=0)[2]/(len(sp)*(3./2.))
    	
    	for s in sp:
    	    rbsp.append(interp1d(s[0], s[1])(rblbd))
    	
    	median = np.median(rbsp, 0)
    	mean = np.mean(rbsp, 0)
    	
    	var_mn = np.var(mean)
        var_md = np.var(median)
        if var_mn < var_md:
            finspec = np.array([rblbd, mean, errors])
        else:
    	   finspec = np.array([rblbd, median, errors])
        
    	snr_ms = []
    	warnings.filterwarnings("ignore")
        for r in cnt_regions:
            try:
                snrval = float(measure_SNR(finspec, r))
                snr_ms.append(snrval)
            except:
                pass
        if len(snr_ms) != 0:
            snr_f = max(snr_ms)
        else:
            snr_f = 0
    	
    	plt.subplot(312)
    	plt.plot(rblbd, mean, 'k-', linewidth=0.3)
    	plt.text(0.05, 0.05, 'mean: var='+str(var_mn), transform=plt.gca().transAxes)
    	
    	plt.subplot(313)
    	plt.plot(rblbd, median, 'k-', linewidth=0.3)
    	plt.text(0.05, 0.05, 'median: var='+str(var_md), transform=plt.gca().transAxes)
    	
    	binf = np.mean([rblbd[i+1]-rblbd[i] for i in range(len(rblbd)-1)])
    	
    	stupname = '_HR'+re.split("hr|HR",(out.split('/')[-2]))[-1]
    	
    	plt.savefig(out+stupname+'.png', dpi=500)
    	plt.clf()
    	
    	obn = out.split('/')[-1]
    	vh = np.mean(vrs['rel'])
    	
    	cbn_mn = open(out+stupname+'.dat', 'w')
        for i in range(len(rblbd)):
            cbn_mn.write('%10.4f%10.4f%10.4f\n' %(finspec[0][i], finspec[1][i], finspec[2][i]))
        cbn_mn.close()
        fits_save(np.array([finspec[0], finspec[1]]), out, stupname, obn, vh, len(specs), snr_f)
    	
    	vr_out = open(out+stupname+'.vr', 'w')
    	vr_out.write(obn)
    	vr_out.write('\nRun             Rel. velocity  Helioc. velocity\n')
    	for i in range(len(vrs['rel'])):
    	    vr_out.write('%-16s%13s%16s\n'%(vrs['run'][i], (vrs['rel'][i]-
    	    	                            vrs['corr'][i]), vrs['rel'][i]))
    	vr_out.write('%-16s%13s%16s'%('mean', '', vh))
    	vr_out.close()
    	
    except:
    	print("problem combining "+out)
    	traceback.print_exc(file=sys.stdout)
    	pass
    
    return out+stupname, vh

def validate(ext, dirname, names):
    
    ext = ext.lower()
    for name in names:
        if name.lower().endswith(ext):
            subdir = dirname.strip('./')
            if not subdir in filstruc:
                filstruc[subdir] = []
                print 'Validating objects directory '+subdir
            with pf.open(subdir+'/'+name) as tstf:
                pid = tstf[0].header['HIERARCH ESO OBS PROG ID']
                try:
                    ix = tstf.index_of('FIBER_SETUP')
                    ob_tab = [list(l) for l in tstf[ix].data]
                    if 'M' in [i[11] for i in ob_tab]:
                        if pid not in prog_mags:
                            prog_mags[pid] = [] 
                        for f in ob_tab:
                            if f[11] == 'M':
                                prog_mags[pid].append([f[7], f[14]])
                except:
                    pass


def incorporate(ext, dirname, names):
    
    ext = ext.lower()
    for name in names:
        if name.lower().endswith(ext):
            subdir = dirname.strip('./')
            with pf.open(subdir+'/'+name) as tstf:
                if ('PROG_ID' in tstf[0].header) and\
                   ('OBJECT' in tstf[0].header):
                    pid = tstf[0].header['PROG_ID']
                    valid_obj = np.transpose(prog_mags[pid])[0]
                    if tstf[0].header['OBJECT'] in valid_obj:
                        filstruc[subdir].append(name)


def structurate(path, xdir):
    
    global filstruc
    global dset_id
    global prog_mags
    
    filstruc = {}
    dsetid = {}
    prog_mags = {}
    exten = '.fits'
    
    os.path.walk(path, validate, exten)
    
    os.path.walk(path, incorporate, exten)
    
    try:
        [os.makedirs(xdir+'/'+a) for a in filstruc.keys()]
    except:
        pass
    
    return filstruc, dsetid, prog_mags


def filxtract(gfspecname, prog, hkeys, dsetid, tr_mag, xdir, tp):
    
    gfspec = pf.open(prog+'/'+gfspecname)
    header = gfspec[0].header
    data = gfspec[1].data
    if not header[hkeys[8]] in dsetid:
        dsetid[prog].append(header[hkeys[8]])
    
    dset = header[hkeys[8]]
    obj = header[hkeys[1]].strip()
    ob_idx = tr_mag[0].index(obj)
    mag = tr_mag[1][ob_idx]
    setup = header[hkeys[4]].strip()
    
    try:
        os.makedirs(xdir+'/'+prog+'/'+dset)
    except:
        pass
    coords = deg2HMS(ra=float(header[hkeys[2]]), dec=float(header[hkeys[3]]))
    
    inspec = np.copy(data[0])
    imge = [xdir+'/'+prog+'/'+dset+'/'+header[hkeys[1]], inspec]
    
    try:
    	#print("fitting "+imge[0]+" attempt 1")
        nr_spec = fit(inspec, imge, ord=6)
    except:
        try:
            #print("fitting "+imge[0]+" attempt 2")
            nr_spec = fit(inspec, imge, ord=4, iter=3)
        except:
            try:
            	#print("fitting "+imge[0]+" attempt 3")
                nr_spec = fit(inspec, imge, ord=2, iter=2)
            except:
                nr_spec = inspec
    
    try:
        cl_spec = clean(nr_spec)
    except:
        try:
            cl_spec = clean(nr_spec, sigma=4)
        except:
            cl_spec = nr_spec
    
    try:
    	#print("start pyfxcor"+dset)
        vr, vrerr = pyfxcor(cl_spec, tp[setup], imge)
    except:
    	#print("something went bad here "+dset)
        vr, vrerr = 0, '0'
    
    c = 299792.458
    try:
        shift = (-vr/c)*np.mean(cl_spec[0])
        wlen = cl_spec[0]+shift
    except:
        wlen = cl_spec[0]
    
    flux = cl_spec[1]
    err = nr_spec[2]
    
    newname = os.path.join(xdir, prog, dset, obj+'.dat')
    xspec = open(newname, 'w')
    xspec.write('#%-15s%-12s%-16.14s%-16.14s%-8.3s%-8s%-10s%-14.12s%-12.9s'
                '%-20.17s%-20.16s%-10.8s%-12.10s%-16.13s\n'%('ProgID', 'Object',
                'RA', 'Dec', 'Mag', 'Setup', 'Plate', 'VHelio corr.',
                'Exp. time', 'Airmass start-end', 'Seeing start-end',
                'SN ratio', 'Resolution', 'Rel. velocity'))
    xspec.write('#%-15s%-12s%-16.14s%-16.14s%-8.5s%-8s%-10s%-14.10s%-12.8s'
                '%9s - %-7s%11s - %-7s%-10.8s%-12s%-16.9s\n'%(
                 header[hkeys[0]].strip(), obj, coords[0], coords[1], mag,
                 setup, header[hkeys[5]].strip(), header[hkeys[6]],
                 header[hkeys[9]], header[hkeys[10]], header[hkeys[11]],
                 header[hkeys[12]], header[hkeys[13]], header[hkeys[14]],
                 header[hkeys[15]], vr))
    xspec.write('#\n#%9s%12s%10s\n' %('lambda', 'flux', 'error'))
    
    for i in range(len(wlen)):
        xspec.write('%10.4f%12.4f%10.4f\n' %(wlen[i]*10, flux[i], err[i]))
    print 'wrote '+newname
    
    try:
        run = str(int(dset[-3:].split('_')[-1]))
    except:
    	try:
        	run = str(int(dset[-3:].split('-')[-1]))
        except:
        	run = "1"
    xspec.close()
    gfspec.close()
    
    sub_ob = "".join(re.split("_|-",header[hkeys[8]])[:-2])+"_"
    
    return [sub_ob+setup, obj, newname]
