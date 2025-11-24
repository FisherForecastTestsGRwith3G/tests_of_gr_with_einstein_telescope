#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H1, L1, V1 sensitivity curves and ranges during O3a and O3b run
"""

import numpy as np
from numpy import *

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['font.size'] = 9
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
matplotlib.rcParams['legend.fontsize'] = 9
import matplotlib.pyplot as plt
plt.rc('axes', axisbelow=True)

colors = {'L1': '#4ba6ff', 'H1': '#ee0000', 'V1': '#9b59b6'}

################################################################################
# Load PSD data
################################################################################
# O3b data was downloaded from the zenodo page, figure2.tar
# See here: https://zenodo.org/records/7997424
#
# O3a data was downloaded from here: https://dcc.ligo.org/LIGO-P2000251/public

psd_data = {}
# Import V1 sensitivity
psd_data["V"] = {}
fvb = 'O3-V1-1265246178_sensitivity_strain_asd.txt'
psd_data["V"]["O3b"] = np.loadtxt(fvb)
fva = 'O3-V1_sensitivity_strain_asd.txt'
psd_data["V"]["O3a"] = np.loadtxt(fva)

psd_data["H"] = {}
# Import H1 sensitivity (C01_CLEAN_SUB60HZ)
fhb = 'O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_strain_asd.txt'
psd_data["H"]["O3b"]= np.loadtxt(fhb) 
fha = "O3-H1-C01_CLEAN_SUB60HZ-1251752040.0_sensitivity_strain_asd.txt"
psd_data["H"]["O3a"] = np.loadtxt(fha)

psd_data["L"] = {}
# Import L1 sensitivity (C01_CLEAN_SUB60HZ)
flb = 'O3-L1-C01_CLEAN_SUB60HZ-1262141640.0_sensitivity_strain_asd.txt'
psd_data["L"]["O3b"] = np.loadtxt(flb)
fla = "O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_strain_asd.txt"
psd_data["L"]["O3a"] = np.loadtxt(fla)

# Sensitivity curves for O3b and O3a
for run in ["O3b", "O3a"]:
    fig = plt.figure()
    ax = fig.add_axes([0,0,1.2,0.8])
    plt.loglog(
        psd_data["L"][run][:,0], 
        psd_data["L"][run][:,1], 
        label=r'$\mathrm{LIGO}\,\mathrm{Livingston}$', color=colors['L1'], linewidth=1, alpha=0.7)
    plt.loglog(
        psd_data["H"][run][:,0], 
        psd_data["H"][run][:,1], 
        label=r'$\mathrm{LIGO}\,\mathrm{Hanford}$', color=colors['H1'], linewidth=1, alpha=0.7)
    plt.loglog(
        psd_data["V"][run][:,0], 
        psd_data["V"][run][:,1], 
        label=r'$\mathrm{Virgo}$', color=colors['V1'], linewidth=1, alpha=0.7)
    plt.legend(loc=(0.12, 0.73), frameon=False, fontsize=12)
    plt.xlabel(r'$\mathrm{Frequency}\,\mathrm{[Hz]}$')
    plt.ylabel(r'$\mathrm{Strain}\,[1/\sqrt{\mathrm{Hz}}]$')
    plt.xlim([10, 4000])
    plt.ylim(1e-24, 1e-18)
    plt.grid(False)
    ax = plt.gca()
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.set_title(f"Sensitivity curves {run}")
    plt.grid(which='both', axis='both',linewidth=0.2)
    plt.savefig(f'LVK_psds_{run}.pdf', bbox_inches='tight')
    plt.close(fig)



