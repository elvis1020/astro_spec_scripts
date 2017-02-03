#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import pyfits


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


def readfits(file):
    
    spectra = pyfits.open(file)
    try:
        bin = spectra[0].header['CDELT1']
    except:
        bin = spectra[0].header['CD1_1']
        print bin
    lambda0 = spectra[0].header['CRVAL1']
    flux = spectra[0].data
    x = np.array([lambda0 + i * bin for i in range(len(flux))])
    y = np.array([i for i in flux])
    
    return [x, y]


def rec_pars(line_id, cols):
    
    line_data[line_id]['lambda'] = float(line_id.split()[1])
    line_data[line_id]['range'] = float(cols[6])
    line_data[line_id]['obpars']['wlcorr'] = float(cols[4])
    line_data[line_id]['obpars']['flxcorr'] = float(cols[5])
    line_data[line_id]['obpars']['style'] = cols[7]
    line_data[line_id]['stpars']['map'] = cols[8]
    line_data[line_id]['stpars']['width'] = float(cols[9])
    line_data[line_id]['stpars']['alpha'] = float(cols[10])


def rec_spec(line_id, cols):
    
    if cols[0] not in file_library.keys():
        file_library[cols[0]] = {'type': 'spectrum', 'status': 'closed'}
    line_data[line_id]['obspec'] = file_library[cols[0]]

    if (cols[1][0] != '@') and (cols[1] != 'none'):
        if cols[1] not in file_library.keys():
            file_library[cols[1]] = {'type': 'spectrum', 'status': 'closed'}
        line_data[line_id]['stspec'] = {cols[1]: file_library[cols[1]]}
    elif cols[1] == 'none':
        line_data[line_id]['stspec'] = 'none'
    else:
        line_data[line_id]['stspec'] = {}
        lst = open(cols[1][1:])
        for coll in (raw.strip().split() for raw in lst):
            if coll[0] not in file_library.keys():
                file_library[coll[0]] = {'type': 'spectrum', 'status': 'closed'}
            line_data[line_id]['stspec'][coll[0]] = file_library[coll[0]]


def build_data(inptfil):
    
    global line_data
    global par_rec
    line_data = {}
    par_rec = []
    with open(inptfil) as l:
        ln = l.readlines()[1:]
    for cols in (raw.strip().split() for raw in ln):
        line_id = cols[2]+' '+cols[3]
        if cols[0][0] != '#':
            par_rec.append(cols)
            line_data[line_id] = {'obpars': {}, 'stpars': {}}
            rec_spec(line_id, cols)
            rec_pars(line_id, cols)
    line_data.keys().sort()
    
    return line_data


def reload_data(inptfil):
    
    changed_lines = []
    with open(inptfil) as l:
        ln = l.readlines()[1:]
    for cols in (raw.strip().split() for raw in ln):
        line_id = cols[2]+' '+cols[3]
        if cols[0][0] != '#' and line_id in line_data.keys() and cols not in par_rec:
            rec_spec(line_id, cols)    
            rec_pars(line_id, cols)
            changed_lines.append(line_id)
            for n, i in enumerate(par_rec):
                if i[3] == cols[3]:
                    par_rec[n] = cols

    return changed_lines


def reload_pars(inptfil):
    
    changed_lines = []
    with open(inptfil) as l:
        ln = l.readlines()[1:]
    for cols in (raw.strip().split() for raw in ln):
        line_id = cols[2]+' '+cols[3]
        if cols[0][0] != '#' and line_id in line_data.keys() and cols not in par_rec:
            rec_pars(line_id, cols)
            changed_lines.append(line_id)
            for n, i in enumerate(par_rec):
                if i[3] == cols[3]:
                    par_rec[n] = cols

    return changed_lines


def update_library(file_library):
    
    for entry in file_library.keys():
        if file_library[entry]['status'] == 'closed':
            if file_library[entry]['type'] == 'spectrum':
                if entry[-5:] == '.fits':
                    file_library[entry]['type'] = 'fits_spectrum'
                else:
                    file_library[entry]['type'] = 'ascii_spectrum'
                    
            if file_library[entry]['type'] == 'fits_spectrum':
                try:
                    file_library[entry]['data'] = {}
                    spec = readfits(entry)
                    file_library[entry]['data']['x'] = spec[0]
                    file_library[entry]['data']['y'] = spec[1]
                    file_library[entry]['status'] = 'open'
                except:
                    print('file '+entry+' not found. skipping...\n')
                    pass
            
            if file_library[entry]['type'] == 'ascii_spectrum':
                try:
                    file_library[entry]['data'] = {}
                    spec = np.genfromtxt(entry, unpack=True)
                    file_library[entry]['data']['x'] = spec[0]
                    file_library[entry]['data']['y'] = spec[1]
                    file_library[entry]['status'] = 'open'
                except:
                    print('file '+entry+' not found. skipping...\n')
                    pass
        
        file_library.keys().sort()


def slicer(line, spec, ctr, rge):
    
    obsx = spec['x']
    x_idx = bissec(obsx, ctr)
    x_m = bissec(obsx, ctr - rge * 0.7)
    x_p = bissec(obsx, ctr + rge * 0.7)
    xslice = spec['x'][x_m:x_p]
    yslice = spec['y'][x_m:x_p]
    
    return [xslice, yslice]


def draw_subplot(line, axn):
    line_data[line]['subplot'] = {'axis': axn}
    ctr = line_data[line]['lambda']
    rge = line_data[line]['range']
    
    obslice = slicer(line, line_data[line]['obspec']['data'], ctr, rge)
    line_data[line]['subplot']['obsplot'] = axn.plot(obslice[0] + line_data[line]['obpars']['wlcorr'],
                                                     obslice[1] / line_data[line]['obpars']['flxcorr'],
                                                     line_data[line]['obpars']['style'])[0]
    
    nplots = len(line_data[line]['stspec'].keys())
    colors = plt.get_cmap(line_data[line]['stpars']['map'])(np.linspace(0, 1.0, nplots+2))
    line_data[line]['subplot']['stplot'] = {}
    for n, static in enumerate(line_data[line]['stspec'].keys()):
        stslice = slicer(line, line_data[line]['stspec'][static]['data'], ctr, rge)
        line_data[line]['subplot']['stplot'][static] = axn.plot(stslice[0], stslice[1], ls='-',
                                                                 lw=line_data[line]['stpars']['width'],
                                                                 alpha=line_data[line]['stpars']['alpha'],
                                                                 color=colors[n+1])[0]
    
    ymin = 999999.
    ymax = -999999.
    for plot in axn.lines:
        plty = plot.get_ydata()
        ymn = min(plty)
        ymx = max(plty)
        if ymn < ymin:
            ymin = ymn
        if ymx > ymax:
            ymax = ymx
    
    axn.set_xlim([ctr - rge / 2., ctr + rge / 2.])
    yrange = ymax - ymin
    axn.set_ylim([ymin - yrange*0.06, ymax + yrange*0.06])
    line_data[line]['subplot']['text'] = axn.text(0.04, 0.06, line ,horizontalalignment='left',
                                                  verticalalignment='bottom', fontsize=12,
                                                  transform = plt.gca().transAxes)


def draw_fig():
    fig = plt.figure()
    nrows = len(line_data.keys()) / 2 + 1
    for n, line in enumerate(sorted(line_data.keys())):
        axn = fig.add_subplot(nrows, 2, n+1)
        draw_subplot(line, axn)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)

    #plt.xlabel('Wavelength [$\AA$]')
    #plt.ylabel('Flux')
    plt.tight_layout()
    plt.ion()
    plt.show()


def redrawer(line):
    line_data[line]['subplot']['obsplot'].remove()
    for static in line_data[line]['subplot']['stplot'].keys():
        line_data[line]['subplot']['stplot'][static].remove()
    line_data[line]['subplot']['text'].remove()
    draw_subplot(line, line_data[line]['subplot']['axis'])


if len(sys.argv) > 1:
    inptfil = sys.argv[1]
else:
    inptfil = 'line_input.dat'
file_library = {}
line_data = build_data(inptfil)
update_library(file_library)
draw_fig()

option = ''
while option != 'q':
    option = raw_input('enter to refresh parameters\nr to reload line data\n'+
    	                'b to rebuild line data (necessary if lambda is change'+
    	                'd)\nq to quit\n')
    if option == '':
        changes = reload_pars(inptfil)
        for line in changes:
            redrawer(line)
    elif option == 'r':
        changes = reload_data(inptfil)
        update_library(file_library)
        for line in changes:
            redrawer(line)
    elif option == 'b':
        plt.close()
        line_data = build_data(inptfil)
        update_library(file_library)
        draw_fig()
    elif option == 'q':
        plt.close()
        exit()
    option = ''
