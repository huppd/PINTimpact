""" runner for swept Hiemenz flow """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterSHLabs.xml')
ROOT = TREE.getroot()

ma.setParameter(ROOT, 'withoutput', 1)
ma.setParameter(ROOT, 'max refinement', 3)
ma.setParameter(ROOT, 'refinement step', 1)
ma.setParameter(ROOT, 'refinement tol', 1.e-4)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*3

RUNS = range(1)

RES = [300]
ALPHA2S = [1./60., 1./30., 1./10.]
ALPHA2S = [1./60., 1./10.]


NPX = 2
NPY = 3
NPZ = 4
NPF = 1


CASE_PATH[0] = pp.DATA_PATH + '/ultimate2'
pp.mkdir(CASE_PATH, 0)

for re in RES:
    CASE_PATH[1] = '/re_'+str(re)
    pp.mkdir(CASE_PATH, 1)
    for alpha2 in ALPHA2S:
        CASE_PATH[2] = '/a2_'+str(int(1./alpha2))
        pp.mkdir(CASE_PATH, 2)
        # for nf in nfs:
        # CASE_PATH[3] = '/nf_'+str(nf)
        # mkdir(CASE_PATH, 3)
        # for nx in nxs:
        # CASE_PATH[4] = '/nx_'+str(nx)
        # mkdir(CASE_PATH, 4)
        # for npx in npxs:
        # CASE_PATH[5] = '/npx_'+str(npx)
        # mkdir(CASE_PATH, 5)
        # for npf in range(1, 2):
        # # for npf in range(1, nf+2):
        # # for npf in range(nf+1, nf+2):
        # CASE_PATH[6] = '/npf_'+str(npf)
        # mkdir(CASE_PATH, 6)
        #
        pp.chdir(CASE_PATH, 2)
        #
        ma.setParameter(ROOT, 'Re', re)
        ma.setParameter(ROOT, 'alpha2', 2.*pi*alpha2*re)
        # ma.setParameter(ROOT, 'nx', 48*NX+1)
        # ma.setParameter(ROOT, 'ny', 96*NX+1)
        # ma.setParameter(ROOT, 'nz', 64*NX+1)
        # ma.setParameter(ROOT, 'nf', NF)
        ma.setParameter(ROOT, 'npx', NPX)
        ma.setParameter(ROOT, 'npy', NPY)
        ma.setParameter(ROOT, 'npz', NPZ)
        ma.setParameter(ROOT, 'npf', NPF)
        TREE.write('parameter3D.xml')
        nptot = NPX*NPY*NPZ*NPF
        memtot = max(1024*8/nptot, 1024)
        print()
        print(CASE_PATH)
        EXE_STRING = pp.exe_pre(nptot, ' -N -W 20:00 -R "rusage[mem=' +
                                str(memtot) + ']" ') + \
            pp.EXE_PATH + '/'+EXE
        print(EXE_STRING)
        os.system(EXE_STRING)
