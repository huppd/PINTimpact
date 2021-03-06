""" runner for Rayleigh streaming """
import os
import xml.etree.ElementTree as ET
from math import pi
import numpy as np
import manipulator as ma
import platform_paths as pp


# load parameter file
ma.set_ids('../XML/parameterStreaming2D.xml')
TREE = ET.parse('../XML/parameterStreaming2D.xml')
ROOT = TREE.getroot()


ma.set_parameter(ROOT, 'withoutput', 1)
ma.set_parameter(ROOT, 'refinement step', 2)
ma.set_parameter(ROOT, 'max refinement', 1)
ma.set_parameter(ROOT, 'refinement tol', 1.e-4)


ma.set_parameter(ROOT, 'lx', 2.)
ma.set_parameter(ROOT, 'ly', 2.)
ma.set_parameter(ROOT, 'nf', 14)

NP = 2

ma.set_parameter(ROOT, 'npx', NP)
ma.set_parameter(ROOT, 'npy', NP)


NXS = [129]

RES = 10**np.linspace(0, 2, 9)[-1:]
STS = 10**np.linspace(-2, 0, 9)[1:2]
STS = STS[::-1]


# make executable ready
EXE = 'peri_navier2D'
os.chdir(pp.EXE_PATH)
os.system('make ' + EXE + ' -j4')

# work direcotries


for tol in [2]:
    CASE_PATH = ['']*4
    CASE_PATH[0] = pp.DATA_PATH + '/streaming_' + str(tol)
    pp.mkdir(CASE_PATH, 0)
    # ma.set_parameter(ROOT, 'Convergence Tolerance', 10**(-tol))
    for nx in NXS:
        CASE_PATH[1] = '/nx_'+str(nx)
        pp.mkdir(CASE_PATH, 1)
        for st in STS:
            CASE_PATH[2] = '/alpha2_'+str(st)
            pp.mkdir(CASE_PATH, 2)
            for re in RES:
                CASE_PATH[3] = '/re_'+str(re)
                pp.mkdir(CASE_PATH, 3)
                pp.chdir(CASE_PATH, 3)
                #
                ma.set_parameter(ROOT, 'nx', nx)
                ma.set_parameter(ROOT, 'ny', nx)
                ma.set_parameter(ROOT, 'Re', re)
                ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*re)
                print(CASE_PATH)
                TREE.write('parameter3D.xml')
                nptot = NP**2
                memtot = int(1024.*max(60./nptot, 2))
                # memtot = 6144
                print()
                print(CASE_PATH)
                EXE_STRING = pp.exe_pre(
                    nptot, ' -N -W 42:00 ' +
                    '-R "rusage[mem=' + str(memtot) +
                    ']" ') + pp.EXE_PATH + '/'+EXE
                # EXE_STRING = pp.exe_pre(
                # nptot, ' -N -W 6:00 ') + pp.EXE_PATH + '/'+EXE
                print(EXE_STRING)
                os.system(EXE_STRING)
