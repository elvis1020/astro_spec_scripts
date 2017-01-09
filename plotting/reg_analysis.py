#!/usr/bin/python

import pyfant as pf
import numpy as np
import matplotlib.pyplot as plt
import copy
import sys, getopt
import os


################################################################################
#defining tasks
################################################################################

def drw_ln(axis, info, y, ytx):
    lbd = info[1]
    plc = (lbd, y)
    txt = (lbd, y - ytx)
    arw = {'width': 0.1, 'headwidth': 0.01, 'headlength': 0.01, 'shrink': 0.1}
    line = axis.annotate('%s %s' % (info[0], lbd), plc, txt, arrowprops=arw, rotation='vertical',
                         ha='right', va='top', size=8)
    props = axis.annotate('%s %s' % (info[2], info[3]), plc, txt, arrowprops=arw, rotation='vertical',
                          ha='left', va='top', size=8)

################################################################################

def drawmarks(ax, lpars):
    deflims = ax.axis()
    y_l = []
    xlims = ax.get_xlim()
    ltrim = []  # list containing line data in the current axis range
    for l in lpars:
        if float(l[1]) > xlims[0] and float(l[1]) < xlims[1]:
            ltrim.append(l)

    for n, l in enumerate(ltrim):
        lbf = float(l[1])
        y = 999999
        for s in ax.lines:
            xdt = s.get_xdata()
            ydt = s.get_ydata()
            if (lbf > min(xdt)) and (lbf < max(xdt)):
                ix = bissec(xdt, lbf)
                yn = min(ydt[ix - 1:ix + 1])
                if yn < y:
                    y = yn  # lowest flux value for the nearest lambda point (nearest to the theoretical lambda)
        ytx = 0.1  # offset of text relative to y value cited above
        if n >= 1:
            span = xlims[1] - xlims[0]  # x axis range
            dst = float(ltrim[n][1]) - float(ltrim[n - 1][1])  # distance between 2 markeable lines
            if (dst / span < 0.02) and (y - ytx - 0.18 < y_l[-1][0] - y_l[-1][1]):
                # if the distance between the current plotting line and the last line is too small, thus making
                # them overlap, then we push down the info text. these values are relative to the window size
                # and a calibration could be useful.
                # in addition, if the line info gets too low, it makes room for the next line in the upper part
                # of the plotting area (y - ytx - 0.28 > y_l[-1][0] - y_l[-1][1])
                ytxn = y - (y_l[-1][0] - y_l[-1][1] - 0.18)  # y - (ylast - textlast - 0.28)
                if ytxn > 0.1:
                    ytx = ytxn
        y_l.append([y, ytx])
        drw_ln(ax, l, y, ytx=ytx)
    ax.axis(deflims)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)

################################################################################

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

################################################################################

def el_multirun(elm):
    
    #separating each line from atom data
    atomtrim = pf.FileAtoms()
    atomtrim.load(atomdata)
    atomtrim.cut(mainfile.llzero, mainfile.llfin)
    atomel = pf.FileAtoms()
    atomel.atoms = [a for a in atomtrim.atoms if a.elem == pf.adjust_atomic_symbol(elm)]
    
    speclines = []
    line_num_specs = []
    
    #configuring multi-running for each line
    el_combos = []
    for atom in atomel.atoms:
        for line in atom.lines:
            
            singline = pf.FileAtoms()
            a = copy.copy(atom)
            a.lines = [line]
            singline.atoms = [a]
    
            cm = pf.Combo([pf.FOR_PFANT, pf.FOR_NULBAD])
            cm.conf.flag_output_to_dir = True
            cm.conf.flag_log_console = False
            cm.conf.opt.no_molecules = True
            cm.conf.opt.absoru = True
            cm.conf.opt.opa = False
            cm.conf.file_main = mainfile
            cm.conf.file_abonds = abondsfile
            cm.conf.file_atoms = singline
            cm.conf.file_main = mainfile
            cm.conf.opt.fn_cv = cm.conf.file_main.flprefix + ".norm.nulbad.%5.3f" % mainfile.fwhm
            
            el_combos.append(cm)
    
    #run pfant-multi
    pf.run_parallel(el_combos, flag_console=False)
    
    #load results and save spectra into lists
    i = 0
    for atom in atomel.atoms:
        for line in atom.lines:
            pl = el_combos[i].pfant
            n = el_combos[i].nulbad
    
            n.load_result()
            n.convolved.filename = str(line.lambda_)
            
            speclines.append(n.convolved)
            
            lbdfluxpec = [n.convolved.x, n.convolved.y]
            line_num_specs.append(lbdfluxpec)
            
            pl.conf.sid.clean()
            
            i += 1
    
    return speclines, line_num_specs

################################################################################

def pfant_cut(no_at=False, no_mol=False, fname='flux'):
    
    #runs pfant cutting down atoms, molecules or neither
    pa = pf.Pfant()
    pa.conf.flag_log_console = False
    pa.conf.opt.no_molecules = no_mol
    pa.conf.opt.no_atoms = no_at
    pa.conf.opt.absoru = True
    pa.conf.opt.opa = False
    pa.conf.file_main = mainfile
    pa.conf.file_abonds = abondsfile
    pa.conf.file_atoms = atomsfile
    pa.run()
    
    na = pf.Nulbad()
    na.conf.flag_log_console = False
    na.conf.file_main = mainfile
    na.run()
    
    na.conf.opt.fn_cv = na.conf.file_main.flprefix + ".norm.nulbad.%5.3f" % mainfile.fwhm
    
    na.load_result()
    spec_cut = na.convolved
    spec_cut.filename = fname
    
    pa.conf.sid.clean()
    
    return spec_cut

################################################################################

def synth_plot_build(ax, name=None, st_dt=None):
#takes a name for the synthesis and the star data dictionary
    
    global mainfile
    global abondsfile
    global atomsfile
    
    #initialize main parameters
    mainfile = pf.FileMain()
    mainfile.load()
    mainfile.pas = 0.01
    
    #setting atmospherical parameters
    if pflag or multiflag:
        mainfile.fwhm = st_dt['fwhm']
        teff = st_dt['teff']
        logg = st_dt['logg']
        met = st_dt['met']
        vt = st_dt['vt']
    if fflag:
        mainfile.fwhm = fw
    
    if pflag or multiflag:
        mainfile.titrav = name
        mainfile.teff = teff
        mainfile.glog = logg
        mainfile.asalog = mainfile.afstar = met
        mainfile.vvt = [vt]
        mainfile.flprefix = name[:3]
    mainfile.nhe = 0.085
    mainfile.llzero = linelambda - linerange / 2 - 3
    mainfile.llfin = linelambda + linerange / 2 + 3
    
    #initializing elemental abundance file and values
    #note that if there's a file "abonds.dat" in the directory, it will be used
    #if not, default values (solar) will be loaded
    #additionally, if multirun is being used, the specified abonds file will be used
    abondsfile = pf.FileAbonds()
    if multiflag:
    	abondsfile.load(st_dt['ab_file'])
    else:
    	abondsfile.load()
    idx = abondsfile.ele.index(pf.adjust_atomic_symbol(element))
    if pf.adjust_atomic_symbol(element) != "FE":
        abondsfile.abol[abondsfile.ele.index("FE")] = -99.99
    abondsfile.abol[idx] = -99.99
    
    #initializing atomic data, based on filename passed or not
    atomsfile = pf.FileAtoms()
    atomsfile.load(atomdata)
    
    #setting and interpolating atmospherical model
    m = pf.Innewmarcs()
    m.conf.flag_log_console = False
    m.conf.file_main = mainfile
    m.conf.opt.absoru = True
    m.conf.opt.opa = False
    m.run()
    
    #calculating hydrogen lines
    h = pf.Hydro2()
    h.conf.flag_log_console = False
    h.conf.file_main = mainfile
    h.run()
    
    #running pfant only for atoms (except chosen element and Fe)
    spec_atoms = pfant_cut(no_at=False, no_mol=True, fname='only_atoms')
    
    
    #running pfant only for molecules
    spec_molecules = pfant_cut(no_at=True, no_mol=False, fname='only_mols')
    
    #running pfant for separate chosen element lines
    #loading back true abundances
    abondsfile = pf.FileAbonds()
    if multiflag:
    	abondsfile.load(st_dt['ab_file'])
    else:
    	abondsfile.load()
    
    if element != 'Fe':
        el_speclines, el_specval = el_multirun(element)
    
    #running pfant for separate iron lines
    fe_speclines, fe_specval = el_multirun("Fe")
    
    #running pfant for the full spectrum
    spec_full = pfant_cut(no_at=False, no_mol=False, fname='full')
    
    #START PLOTTING
    value_specat = [spec_atoms.x, spec_atoms.y]
    value_specmol = [spec_molecules.x, spec_molecules.y]
    value_total = [spec_full.x, spec_full.y]
    
    ax.plot(value_total[0], value_total[1], 'b-')
    ax.plot(value_specat[0], value_specat[1], 'g-')
    ax.plot(value_specmol[0], value_specmol[1], 'c-')
    
    if sflag or multiflag:
        if multiflag:
            global ospec
            ospec = np.genfromtxt(st_dt['spec'], unpack=True)
            ax.plot(ospec[0]+st_dt['sp_adj'][0], ospec[1]/st_dt['sp_adj'][1], 'k--')
        elif rflag:
            ax.plot(ospec[0]+rel[0], ospec[1]/rel[1], 'k--')
        else:
            ax.plot(ospec[0], ospec[1], 'k--')
    for l in fe_specval:
        ax.plot(l[0], l[1], color='#FF8C00', ls='-.', lw=1.5)
    if element != 'Fe':
        for l in el_specval:
            ax.plot(l[0], l[1], 'm-.', lw=1.5)
    if multiflag:
    	plt.text(0.04, 0.06, name, horizontalalignment='left', 
    	     	verticalalignment='bottom', fontsize=16,
    	     	transform=plt.gca().transAxes)
    plt.ylim((0., 1.05))


################################################################################
#pipeline start
################################################################################

global pflag
global sflag
global rflag
global fflag
global multiflag
pflag = False
sflag = False
rflag = False
fflag = False
multiflag = False

################argument setting definition
opts, args = getopt.getopt(sys.argv[1:], "s:r:f:", ["obspec=", "relocation=", "fwhm="])

argno = len(args)

if argno == 3:
    element = args[0]
    linelambda = float(args[1])
    linerange = float(args[2])
elif argno < 3 and argno != 1:
    print("must give at least 3 arguments:\nelement, lambda, lambda range")
    exit()
elif argno < 3 and argno == 1:
    print("reading parameters and file lists from " + args[0])
    parfile = args[0]
    multiflag = True
elif (argno > 3 and argno < 8) or argno > 8:
    print(
        "if you are giving more than 3 arguments (element, lambda, range), please supply atmospherical parameters:\n\
        synth name, Teff, logG, [Fe/H], m.turb.velocity")
    exit()
elif argno == 8:
    star_data = {}
    st_name = args[3]
    star_data[st_name] = {}
    star_data[st_name]['teff'] = float(args[4])
    star_data[st_name]['logg'] = float(args[5])
    star_data[st_name]['met'] = float(args[6])
    star_data[st_name]['vt'] = float(args[7])
    element = args[0]
    linelambda = float(args[1])
    linerange = float(args[2])
    pflag = True

################options settings
for opt, arg in opts:
    if opt == '-s':
        global ospec
        ospec = np.genfromtxt(arg, unpack=True)
        sflag = True
    if opt == '-r':
        rel = [float(arg.split(',')[0]), float(arg.split(',')[1])]
        rflag = True
    if opt == '-f':
        fw = float(arg)
        fflag = True

################multi synth mode settings
if multiflag:
    star_data = {}
    pardata = open(parfile).readlines()
    global_pars = pardata[1].strip().split()
    element = global_pars[0]
    linelambda = float(global_pars[1])
    linerange = float(global_pars[2])
    atomdata = global_pars[3]
    for cols in [raw.strip().split() for raw in pardata[3:]]:
        if len(cols) > 3:
            st_name = cols[0]
            star_data[st_name] = {}
            star_data[st_name]['teff'] = float(cols[1])
            star_data[st_name]['logg'] = float(cols[2])
            star_data[st_name]['met'] = float(cols[3])
            star_data[st_name]['vt'] = float(cols[4])
            star_data[st_name]['fwhm'] = float(cols[5])
            star_data[st_name]['spec'] = cols[6]
            star_data[st_name]['sp_adj'] = (float(cols[7]), float(cols[8]))
            star_data[st_name]['ab_file'] = cols[10]
            if cols[9] == 'T':
                star_data[st_name]['mark'] = True
            else:
                star_data[st_name]['mark'] = False
else:
    atomdata = "atoms.dat"

################generating line properies for marking
atfil = open(atomdata)
cols = [raw.strip().split() for raw in atfil]
trim = cols
lineprops = []

for n, i in enumerate(trim):
    if (len(i) <= 4):
        el = i[0][:-1]
        if el == pf.adjust_atomic_symbol(element):
            lineprops.append([i[0], i[1], trim[n + 1][0], trim[n + 1][1]])

################plotting
fig = plt.figure()

if multiflag:
    stars = star_data.keys()
    subcode = len(stars)*100
    for n, s in enumerate(stars):
        print('working on '+s)
        if n == 0:
            star_data[s]['subpt'] = fig.add_subplot(subcode+n+11)
        else:
            star_data[s]['subpt'] = fig.add_subplot(subcode+n+11, sharex=star_data[stars[0]]['subpt'])
        synth_plot_build(star_data[s]['subpt'], name=s, st_dt=star_data[s])
        plt.xlim((mainfile.llzero+2, mainfile.llfin-2))
        if star_data[s]['mark']:
            drawmarks(star_data[s]['subpt'], lineprops)
elif pflag:
    ax = fig.add_subplot(111)
    synth_plot_build(ax, name=st_name, st_dt=star_data[st_name])
    plt.xlim((mainfile.llzero+2, mainfile.llfin-2))
    drawmarks(ax, lineprops)
else:
    ax = fig.add_subplot(111)
    synth_plot_build(ax)
    plt.xlim((mainfile.llzero+2, mainfile.llfin-2))
    drawmarks(ax, lineprops)

mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.tight_layout()
plt.show()

os.system("rm -r session*")
