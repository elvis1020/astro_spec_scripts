import os
import numpy as np
from astropy.io import fits
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel, ConstantModel, SkewedGaussianModel, PolynomialModel
from scipy.interpolate import splrep, splev, interp1d
from scipy.ndimage import gaussian_filter

rej_ob = []


def readfits(name):
    hdu = fits.open(name)
    head = hdu[0].header
    flux = hdu[0].data
    wavest = head['CRVAL1']
    wbin = head['CDELT1']
    wave = [wavest + wbin * i for i in range(len(flux))]
    return [wave, flux]


def fits_save(spec, filen, obnm, n):
    btspec = spec[1].astype('float32')
    hdu = fits.PrimaryHDU(btspec)
    delt = np.mean([spec[0][i + 1] - spec[0][i] for i in range(len(spec[0]) - 1)])
    hdu.header['CTYPE1'] = 'Wavelength'
    hdu.header['BUNIT'] = 'Flux'
    hdu.header['CRPIX1'] = 1.
    hdu.header['CRVAL1'] = spec[0][0]
    hdu.header['CDELT1'] = delt
    hdu.header['CD1_1'] = delt
    hdu.header['OBJ_NAME'] = obnm
    hdu.header['N_SPEC'] = n
    # hdu.header['SNR'] = snr
    name = filen + '.fits'
    try:
        hdu.writeto(name)
    except:
        os.remove(name)
        hdu.writeto(name)


def bissec(A, x):
    n = len(A) - 1
    if (x < A[0]):
        return 0
    elif (x > A[n]):
        return n + 1
    n1 = 0
    n2 = n
    while ((n2 - n1) > 1):
        nm = (n1 + n2) / 2
        if ((x - A[nm]) > 0):
            n1 = nm
        else:
            n2 = nm
    return n1


def measure_SNR(spec, region):
    idx1 = (np.abs(spec[0] - region[0])).argmin()
    idx2 = (np.abs(spec[0] - region[1])).argmin()
    reg = np.array(spec[1][idx1:idx2])
    if (len(reg) != 0) and (np.std(reg) > 1e-10):
        return np.mean(reg) / np.std(reg)
    else:
        return 0


def deg2HMS(ra='', dec='', round=False):
    RA, DEC, rs, ds = '', '', '', ''
    if dec:
        if str(dec)[0] == '-':
            ds, dec = '-', abs(dec)
        deg = int(dec)
        decM = abs(int((dec - deg) * 60))
        if round:
            decS = int((abs((dec - deg) * 60) - decM) * 60)
        else:
            decS = (abs((dec - deg) * 60) - decM) * 60
        DEC = '{0}{1}:{2}:{3}'.format(ds, deg, decM, decS)
    if ra:
        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra / 15)
        raM = int(((ra / 15) - raH) * 60)
        if round:
            raS = int(((((ra / 15) - raH) * 60) - raM) * 60)
        else:
            raS = ((((ra / 15) - raH) * 60) - raM) * 60
        RA = '{0}{1}:{2}:{3}'.format(rs, raH, raM, raS)
    if ra and dec:
        return RA, DEC
    else:
        return RA or DEC


def fit(spectra, sigma=4.0, ptreg=8., ord=4, iter=2):
    rej = []
    for i, w in enumerate(spectra[1]):
        if w <= 0:
            rej.append(i)
    spectra[0] = np.delete(spectra[0], rej)
    spectra[1] = np.delete(spectra[1], rej)

    # prepare first kick
    poly = PolynomialModel(3)
    pars = poly.make_params()
    for p in range(4):
        label = 'c' + str(p)
        pars[label].set(value=1., vary=True)
    wkcopy = np.copy(spectra[1])
    truesp = [i for i in wkcopy if i >= 0]
    truex = [spectra[0][i] for i in range(len(spectra[1])) if spectra[1][i] >= 0]
    outcont = poly.fit(truesp, pars, x=truex)
    firstcont = outcont.eval(x=spectra[0])

    xn = np.copy(spectra[0])
    yn = np.copy(spectra[1]) / firstcont

    # start cont. cleaning iterations
    for i in range(iter):
        i_ = np.copy(i)
        niter = str(i_ + 1)
        sigma = sigma - i * 0.21 * sigma

        md = np.median(yn)
        n = len([i for i in yn if i > 0.1])
        offset = (len(xn) - n) / 2
        absor = md - min(yn[offset:n - offset])
        freq, bin = np.histogram(yn, bins=50, range=(md - absor, md + absor))
        rebin = [(bin[b + 1] + bin[b]) / 2 for b in range(len(bin) - 1)]

        gauss = SkewedGaussianModel()
        pars = gauss.make_params()
        pars['center'].set(vary=True)
        pars['amplitude'].set(vary=True)
        pars['sigma'].set(vary=True)
        pars['gamma'].set(vary=True)
        out = gauss.fit(freq, pars, x=rebin)

        var = sigma * out.best_values['sigma']
        xrbn = np.linspace(rebin[0], rebin[-1], num=100)
        yrbn = list(out.eval(x=xrbn))
        mode = xrbn[yrbn.index(max(yrbn))]

        # clean cont.
        ync = np.copy(spectra[1])
        xnc = np.copy(spectra[0])

        mask = []
        for j in range(len(yn)):
            if (yn[j] > mode + var / 2) or (yn[j] < mode - var / 2):
                mask.append(False)
            else:
                mask.append(True)
        mask = np.array(mask)
        ync = ync[mask]
        xnc = xnc[mask]

        # re-fitting
        poly2 = PolynomialModel(ord)
        pars2 = poly2.make_params()
        for p in range(ord + 1):
            label = 'c' + str(p)
            pars2[label].set(value=1., vary=True)
        try:
            outcont2 = poly2.fit(ync, pars2, x=xnc)
        except:
            plt.plot(xn, yn, 'k-')
            plt.plot([xn[0], xn[-1]], [mode, mode], 'b-')
            plt.plot([xn[0], xn[-1]], [mode + var / 2, mode + var / 2], 'r-')
            plt.plot([xn[0], xn[-1]], [mode - var / 2, mode - var / 2], 'r-')
            plt.show()

        contf = outcont2.eval(x=xn)
        yn = spectra[1] / contf

    clspec = [xnc, ync]

    # start slicing
    firstv = clspec[0][0]
    wavrange = clspec[0][-1] - firstv
    sliceno = wavrange / ptreg
    slisize = wavrange / sliceno
    points = [[], []]

    # continuum point definition
    for s in range(int(sliceno)):
        i = bissec(clspec[0], firstv + s * slisize)
        f = bissec(clspec[0], firstv + (s + 1) * slisize)
        slc = [clspec[0][i:f], clspec[1][i:f]]
        if len(slc[1]) > 2.:
            md = np.median(slc[1])
            absor = min(slc[1])
            high = max(slc[1])
            freq, bin = np.histogram(slc[1], bins=20, range=(absor, high))
            rebin = [(bin[b + 1] + bin[b]) / 2 for b in range(len(bin) - 1)]

            fmode = rebin[list(freq).index(max(freq))]
            fsigma = rebin[-1] - rebin[0]

            gauss = GaussianModel()
            pars = gauss.make_params()
            pars['center'].set(value=fmode, vary=True)
            pars['amplitude'].set(value=max(freq), vary=True)
            pars['sigma'].set(value=fsigma, vary=True)
            out = gauss.fit(freq, pars, x=rebin)

            xrbn = np.linspace(rebin[0], rebin[-1], num=100)
            yrbn = list(out.eval(x=xrbn))
            mode = xrbn[yrbn.index(max(yrbn))]
            xp = slc[0][len(slc[0]) / 2]
            points[0].append(xp)
            points[1].append(mode)

    spline = splrep(points[0], points[1], k=3)
    contx = splev(clspec[0], spline)
    continuum = splev(spectra[0], spline)

    return [spectra[0], spectra[1] / continuum]


def clean(spectra, sigma=4.5):
    rej = []
    for i, w in enumerate(spectra[1]):
        if w <= 0:
            rej.append(i)
    spectra[0] = np.delete(spectra[0], rej)
    spectra[1] = np.delete(spectra[1], rej)

    md = np.median(spectra[1])
    n = int(len(spectra[0]) * 0.8)
    offset = (len(spectra[0]) - n) / 2
    absor = md - min(spectra[1][offset:n - offset])
    freq, bin = np.histogram(spectra[1], bins=30, range=(md - absor, md + absor))
    rebin = [(bin[b + 1] + bin[b]) / 2 for b in range(len(bin) - 1)]

    gauss = SkewedGaussianModel()
    pars = gauss.make_params()
    pars['center'].set(vary=True)
    pars['amplitude'].set(vary=True)
    pars['sigma'].set(vary=True)
    pars['gamma'].set(vary=True)
    out = gauss.fit(freq, pars, x=rebin)

    var = sigma * out.best_values['sigma']
    xrbn = np.linspace(rebin[0], rebin[-1], num=100)
    yrbn = list(out.eval(x=xrbn))
    mode = xrbn[yrbn.index(max(yrbn))]

    xn = np.copy(spectra[0])
    yn = np.copy(spectra[1])

    var = mode

    for i, w in enumerate(yn):
        if w > mode + var:
            rej.append(i)
    xn = np.delete(xn, rej)
    yn = np.delete(yn, rej)

    return np.array([xn, yn])


def pyfxcor(inspec, template, vmin=-200., vmax=200., res=3, rej=1000):
    rv, cc = pyasl.crosscorrRV(inspec[0], inspec[1], template[0], template[1], vmin, vmax, res, skipedge=rej)

    cen_gs = np.argmax(cc)
    perfx, perfy = rv[cen_gs - 5:cen_gs + 6], cc[cen_gs - 5:cen_gs + 6]

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
        return 'error', ''

    return ct, cterr


def combiner(specs):
    sp = specs

    N = np.ceil(np.mean([len(s[0]) for s in sp]))
    lmin = max([s[0][0] for s in sp])
    lmax = min([s[0][-1] for s in sp])
    rblbd = np.linspace(lmin, lmax, N)
    rbsp = []

    for s in sp:
        rbsp.append(interp1d(s[0], s[1])(rblbd))

    median = np.median(rbsp, 0)
    mean = np.mean(rbsp, 0)

    var_mn = np.var(mean)
    var_md = np.var(median)
    if var_mn < var_md:
        finspec = [rblbd, mean]
    else:
        finspec = [rblbd, median]

    return finspec


def structurate(ext, dirname, names):
    ext = ext.lower()
    for filein in names:
        if filein.lower().endswith(ext):
            filein = os.path.join(dirname, filein).strip('./')
            f = fits.open(filein)
            head = f[0].header

            prog = head['HIERARCH ESO OBS PROG ID']
            ob = head['HIERARCH ESO OBS NAME']

            if ob not in rej_ob:

                pplfil = head['PIPEFILE']
                airm = np.mean([head['HIERARCH ESO TEL AIRM START'], head['HIERARCH ESO TEL AIRM END']])
                sng = np.mean([head['HIERARCH ESO TEL AMBI FWHM START'], head['HIERARCH ESO TEL AMBI FWHM END']])
                pla = head['HIERARCH ESO INS OBSPLATE']
                expt = head['EXPTIME']
                dat = head['DATE-OBS'].split('T')[0]
                tim = head['DATE-OBS'].split('T')[1]
                lat = head['HIERARCH ESO TEL GEOLAT']
                lon = head['HIERARCH ESO TEL GEOLON']
                elv = head['HIERARCH ESO TEL GEOELEV']
                ptra = head['RA']
                ptdec = head['DEC']
                jld = head['MJD-OBS']
                helicorr = pyasl.helcorr(lon, lat, elv, ptra, ptdec, jld + 2.4e6)[0]

                if prog not in filstruc.keys():
                    filstruc[prog] = {}
                    infos[prog] = {}

                if ob not in filstruc[prog].keys():
                    filstruc[prog][ob] = {}
                    infos[prog][ob] = {}

                if 'bin_table_info' in pplfil:
                    tab = f[1].data
                    assoc = {}
                    coords = {}
                    for l in tab:
                        if l[3] != '':
                            assoc[l[2] + 1] = l[3]
                            coords[l[3]] = deg2HMS(l[4], l[5])
                    filstruc[prog][ob]['fibassoc'] = assoc
                    infos[prog][ob]['coords'] = coords
                    infos[prog][ob]['airmass'] = airm
                    infos[prog][ob]['seeing'] = sng
                    infos[prog][ob]['plate'] = pla
                    infos[prog][ob]['exptime'] = expt
                    infos[prog][ob]['date'] = dat
                    infos[prog][ob]['time'] = tim
                    infos[prog][ob]['vhelio'] = helicorr
                    infos[prog][ob]['rv'] = {}

                elif ('mwfxb' in pplfil) and ('sigma' not in pplfil):
                    if head['EXTNAME'].strip() == 'CCD-44':
                        chip = 'lower'
                    else:
                        chip = 'upper'
                    fibre = int(pplfil.strip('.fits').split('_')[2])
                    if fibre not in filstruc[prog][ob].keys():
                        filstruc[prog][ob][fibre] = {}

                    flux = f[0].data
                    wavest = head['CRVAL1']
                    wbin = head['CDELT1']
                    wave = [wavest + wbin * i for i in range(len(flux))]
                    filstruc[prog][ob][fibre][chip] = [wave, flux]


filstruc = {}
infos = {}
comb = {}
exten = '.fits'
path = '.'
suntpl = np.genfromtxt('solar580.dat', unpack=True)
os.path.walk(path, structurate, exten)

for pid in filstruc.keys():
    print "working on dataset ", pid
    for ob in filstruc[pid].keys():
        print "OB name: ", ob
        for fibre in filstruc[pid][ob].keys():
            if fibre != 'fibassoc':
                fibnam = filstruc[pid][ob]['fibassoc'][fibre]
                print "fibre no. ", fibre, "; object: ", fibnam

                # lower chip
                print "doing lower part"
                spectra = filstruc[pid][ob][fibre]['lower']
                fitspec = fit(spectra, ptreg=8.)
                cleanspec = clean(fitspec)
                rv, err = pyfxcor(cleanspec, suntpl)
                c = 299792.458
                shift = (-rv / c) * np.mean(cleanspec[0])
                wlen = cleanspec[0] + shift
                lower = [wlen, cleanspec[1]]

                # upper chip
                print "doing upper part"
                spectra = filstruc[pid][ob][fibre]['upper']
                fitspec = fit(spectra, ptreg=20.)
                cleanspec = clean(fitspec)
                rv, err = pyfxcor(cleanspec, suntpl)
                c = 299792.458
                shift = (-rv / c) * np.mean(cleanspec[0])
                wlen = cleanspec[0] + shift
                upper = np.array([wlen, cleanspec[1]])
                if fibnam not in infos[pid][ob]['rv']:
                    infos[pid][ob]['rv'][fibnam] = rv

                if fibnam not in comb.keys():
                    comb[fibnam] = {'lower': [], 'upper': []}

                comb[fibnam]['lower'].append(lower)
                comb[fibnam]['upper'].append(upper)
        print "OB", ob, "done."

for star in comb.keys():

    infof = open(star + '_infotab.dat', 'w')
    pid = infos.keys()[0]
    infof.write('%-16s%-12s%-14s%-14s\n' % ('Program ID', 'Object', 'RA', 'DEC'))
    infof.write('%-16s%-12s%-14.12s%-14.12s\n\n' % (pid, star, infos[pid][ob]['coords'][star][0],
                                                    infos[pid][ob]['coords'][star][1]))

    infof.write(
        '%-19s%-9s%-9s%-9s%-8s%-7s%-11s%-25s\n' % ('Observation block', 'RV', 'Vhelio', 'airmass', 'seeing', 'plate',
                                                   'exp.time', 'exposure date/time'))
    vhelios = []
    for ob in infos[pid].keys():
        rv = infos[pid][ob]['rv'][star]
        vhl = rv + infos[pid][ob]['vhelio']
        infof.write('%-19s%-9.2f%-9.2f%-9.3f%-8.3f%-7s%-11.4f%-25s\n' % (ob, rv, vhl, infos[pid][ob]['airmass'],
                                                                         infos[pid][ob]['seeing'],
                                                                         str(infos[pid][ob]['plate']),
                                                                         infos[pid][ob]['exptime'],
                                                                         infos[pid][ob]['date'] + '/' + infos[pid][ob][
                                                                             'time']))
        vhelios.append(vhl)
    infof.write('Mean vhelio: %7.2f' % (np.mean(vhelios)))
    infof.close()

for star in comb.keys():
    finl = combiner(comb[star]['lower'])
    fillg = gaussian_filter(finl[1], 1.0)
    finu = combiner(comb[star]['upper'])
    filug = gaussian_filter(finu[1], 1.0)
    num = len(comb[star]['lower'])
    fits_save([finl[0], fillg], star + '_L', star, num)
    fits_save([finu[0], filug], star + '_U', star, num)
    with open(star + '.dat', 'w') as otp:
        otp.write('#Object: %-10s' % (star))
        for i in range(len(finl[0])):
            otp.write('\n%9.4f  %6.4f' % (finl[0][i], fillg[i]))
        for i in range(len(finu[0])):
            otp.write('\n%9.4f  %6.4f' % (finu[0][i], filug[i]))
