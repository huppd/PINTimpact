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
ma.setParameter(ROOT, 'refinement step', 1)
ma.setParameter(ROOT, 'refinement tol', 1.e-4)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*3

RUNS = range(1)

RES = [300]
STS = [1./60., 1./30., 1./10.]
STS = [1./60., 1./10.]
STS = [1./60.]

NFS = [0, 1, 2]

NPX = 1
NPY = 2
NPZ = 2
NPF = 1


CASE_PATH[0] = pp.DATA_PATH + '/ultimate3'
pp.mkdir(CASE_PATH, 0)

for re in RES:
    for nf in NFS:
        CASE_PATH[1] = '/nf_'+str(nf)
        pp.mkdir(CASE_PATH, 1)
        for st in STS:
            CASE_PATH[2] = '/a2_'+str(int(1./st))
            pp.mkdir(CASE_PATH, 2)
            pp.chdir(CASE_PATH, 2)
            #
            ma.setParameter(ROOT, 'Re', re)
            ma.setParameter(ROOT, 'alpha2', 2.*pi*st*re)
            ma.setParameter(ROOT, 'lx', 18.)
            ma.setParameter(ROOT, 'ly', 100.)
            ma.setParameter(ROOT, 'lz', 40.)
            ma.setParameter(ROOT, 'origin z', 20.)
            ma.setParameter(ROOT, 'nx', 31)
            ma.setParameter(ROOT, 'ny', 73)
            ma.setParameter(ROOT, 'nz', 65)
            ma.setParameter(ROOT, 'nf', nf)
            ma.setParameter(ROOT, 'max refinement', 3-nf)
            ma.setParameter(ROOT, 'npx', NPX)
            ma.setParameter(ROOT, 'npy', NPY)
            ma.setParameter(ROOT, 'npz', NPZ)
            ma.setParameter(ROOT, 'npf', NPF)
            TREE.write('parameter3D.xml')
            nptot = NPX*NPY*NPZ*NPF
            memtot = max(1024*8/nptot, 1024)
            print()
            print(CASE_PATH)
            EXE_STRING = pp.exe_pre(nptot, ' -N -W 24:00 -R "rusage[mem=' +
                                    str(memtot) + ']" ') + \
                pp.EXE_PATH + '/'+EXE
            print(EXE_STRING)
            os.system(EXE_STRING)
