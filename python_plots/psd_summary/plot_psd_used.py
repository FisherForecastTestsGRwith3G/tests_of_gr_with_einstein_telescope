#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H1, L1, V1 sensitivity curves and ranges during O3b run. 
ET sensitivity curves and ranges for 10km and 15km arms.  
"""

import numpy as np
from numpy import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.size'] = 9
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['legend.fontsize'] = 9
import matplotlib.pyplot as plt
plt.rc('axes', axisbelow=True)

colors = {'L1': '#4AA553', 'H1': "#BF70CF", 'V1':"#A98E2A", "ET10":"#DD6E4D", "ET15":"#339AF6"}

################################################################################
# Load PSD data 
################################################################################
# O3b data was downloaded from the zenodo page, figure2.tar
# See here: https://zenodo.org/records/7997424 

asd_data = {}
fvb = 'O3-V1-1265246178_sensitivity_strain_asd.txt'
asd_data["V1"] = np.loadtxt(fvb)

# Import H1 sensitivity (C01_CLEAN_SUB60HZ)
fhb = 'O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_strain_asd.txt'
asd_data["H1"]= np.loadtxt(fhb) 

# Import L1 sensitivity (C01_CLEAN_SUB60HZ)
flb = 'O3-L1-C01_CLEAN_SUB60HZ-1262141640.0_sensitivity_strain_asd.txt'
asd_data["L1"] = np.loadtxt(flb)


psd_data = {}
# Import ET10km sensitivity 
f10 = 'ET10km.txt'
psd_data["ET10"] = np.loadtxt(f10)

# Import ET15km sensitivity
f15 = 'ET15km.txt'
psd_data["ET15"] = np.loadtxt(f15)

y_lims = [1e-25, 1e-16]
lw = 1.5
anot = 0.4


fig = plt.figure()
ax = fig.add_axes([0,0,1,0.8])

xrange_lvk = [20,5000]
idx_select = asd_data["L1"][:,0]>20
plt.loglog(
    asd_data["L1"][:,0][idx_select], 
    asd_data["L1"][:,1][idx_select], 
    label=r'$\mathrm{LIGO\,O3b}\,\mathtt{LHO}$', color=colors['L1'], linewidth=lw, alpha=1)
plt.loglog(
    asd_data["L1"][:,0][np.invert(idx_select)], 
    asd_data["L1"][:,1][np.invert(idx_select)], 
    label=None, color=colors['L1'], linewidth=lw, alpha=anot)

idx_select = asd_data["H1"][:,0]>20
plt.loglog(
    asd_data["H1"][:,0][idx_select], 
    asd_data["H1"][:,1][idx_select], 
    label=r'$\mathrm{LIGO\,O3b}\,\mathtt{LHO}$', color=colors['H1'], linewidth=lw, alpha=1)
plt.loglog(
    asd_data["H1"][:,0][np.invert(idx_select)], 
    asd_data["H1"][:,1][np.invert(idx_select)], 
    label=None,color=colors['H1'], linewidth=lw, alpha=anot)

idx_select = np.logical_and(asd_data["V1"][:,0] > 20,asd_data["V1"][:,0] < 5000 )
plt.loglog(
    asd_data["V1"][:,0][idx_select], 
    asd_data["V1"][:,1][idx_select], 
    label=r'$\mathrm{Virgo\,O3b}\,\mathtt{V}$', color=colors['V1'], linewidth=lw, alpha=1)
idx_select = np.logical_and(asd_data["V1"][:,0] > 10,asd_data["V1"][:,0] < 10 )
plt.loglog(
    asd_data["V1"][:,0][idx_select], 
    asd_data["V1"][:,1][idx_select], 
    label=None, color=colors['V1'], linewidth=lw, alpha=anot)
idx_select = asd_data["V1"][:,0] > 5000
plt.loglog(
    asd_data["V1"][:,0][idx_select], 
    asd_data["V1"][:,1][idx_select], 
    label=None, color=colors['V1'], linewidth=lw, alpha=anot)

idx_select = psd_data["ET10"][:,0] > 2
plt.loglog(
    psd_data["ET10"][:,0][idx_select], 
    np.sqrt(psd_data["ET10"][:,1][idx_select]), 
    label=r'$\mathrm{ET} 10\mathrm{km}$', color=colors['ET10'], linewidth=lw, alpha=1)
plt.loglog(
    psd_data["ET10"][:,0][np.invert(idx_select)], 
    np.sqrt(psd_data["ET10"][:,1][np.invert(idx_select)]), 
    label=None, color=colors['ET10'], linewidth=lw, alpha=anot)

idx_select = psd_data["ET15"][:,0] > 2
plt.loglog(
    psd_data["ET15"][:,0][idx_select], 
    np.sqrt(psd_data["ET15"][:,1][idx_select]), 
    label=r'$\mathrm{ET} 15\mathrm{km}$', color=colors['ET15'], linewidth=lw, alpha=1)
plt.loglog(
    psd_data["ET15"][:,0][np.invert(idx_select)], 
    np.sqrt(psd_data["ET15"][:,1][np.invert(idx_select)]), 
    label=None, color=colors['ET15'], linewidth=lw, alpha=anot)

plt.legend(loc=(0.7, 0.64), frameon=True, fontsize=12)
plt.xlabel(r'$\mathrm{Frequency}\,\mathrm{[Hz]}$',fontsize=12)
plt.ylabel(r'$\mathrm{Strain}\,[1/\sqrt{\mathrm{Hz}}]$', fontsize=12)
plt.xlim([0.8, 12000])
plt.ylim(y_lims[0],y_lims[1])
ax = plt.gca()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.grid(which='both', axis='both',linewidth=0.2)
plt.savefig(f'psds_used.pdf', bbox_inches='tight')
plt.close(fig)



