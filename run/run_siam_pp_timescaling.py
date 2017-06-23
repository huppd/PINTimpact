""" runner for SIAM PP timescaling """
import os
import xml.etree.ElementTree as ET
from math import pi
import numpy as np

import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3D.xml')
ROOT = TREE.getroot()

ma.setParameter(ROOT, 'withoutput', 0)
ma.setParameter(ROOT, 'refinement level', 1)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*10

RUNS = range(10, 12)

RE = 100
ST = 0.1
NFS = 2**np.arange(1, 5) - 1
NFS = [15]
NX = 2
NPX = 4


CASE_PATH[0] = '/scale_time_re_' + str(RE) + '_a2_' + str(ST) + '_nx_' + \
    str(NX) + '_npx_' + str(NPX)

pp.mkdir(CASE_PATH, 0)

for nf in NFS:
    CASE_PATH[3] = '/nf_'+str(nf)
    pp.mkdir(CASE_PATH, 3)
    for npf in 2**np.arange(4, 5):
        if npf <= (nf+1) and nf/npf <= 15:
            CASE_PATH[6] = '/npf_'+str(npf)
            pp.mkdir(CASE_PATH, 6)
            #
            pp.chdir(CASE_PATH, 6)
            #
            ma.setParameter(ROOT, 'Re', RE)
            ma.setParameter(ROOT, 'alpha2', 2.*pi*ST*RE)
            ma.setParameter(ROOT, 'nx', 64*NX+1)
            ma.setParameter(ROOT, 'ny', 16*NX+1)
            ma.setParameter(ROOT, 'nz', 32*NX+1)
            ma.setParameter(ROOT, 'nf', nf)
            ma.setParameter(ROOT, 'npx', NPX)
            ma.setParameter(ROOT, 'npy', max(NPX/2, 1))
            ma.setParameter(ROOT, 'npz', NPX)
            ma.setParameter(ROOT, 'npf', npf)
            TREE.write('parameter3D.xml')
            nptot = NPX*NPX*max(NPX/2, 1)*npf
            for run in RUNS:
                print()
                print(CASE_PATH)
                exe_str = pp.exe_pre(nptot,
                                     ' -N -R "select[model==Opteron6174"]  ',
                                     run) + pp.EXE_PATH+'/'+EXE
                print(exe_str)
                os.system(exe_str)
