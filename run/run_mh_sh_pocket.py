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

NPX = 3
NPY = 2
NPZ = 2
NPF = 1

LXO = 22.5
LYO = 600.
LZO = 150.
NXO = 97
NYO = 1537
NZO = 513

# NX = NXO
# NX = (65-1)/2 +1
# NY = (129-1)/2 + 1
# NZ = (65-1)/2 + 1
NX = 31
NY = 73
NZ = 65

LX = round(2.5*LXO/(NXO-1)*(NX-1))
LY = round(4.*LYO/(NYO-1)*(NY-1))
LZ = round(2.*LZO/(NZO-1)*(NZ-1))

print('LX', LX)
print('LY', LY)
print('LZ', LZ)

CASE_PATH[0] = pp.DATA_PATH + '/ultimate'
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
            ma.setParameter(ROOT, 'lx', LX)
            ma.setParameter(ROOT, 'ly', LY)
            ma.setParameter(ROOT, 'lz', LZ)
            ma.setParameter(ROOT, 'origin z', LZ/2.)
            ma.setParameter(ROOT, 'nx', NX)
            ma.setParameter(ROOT, 'ny', NY)
            ma.setParameter(ROOT, 'nz', NZ)
            ma.setParameter(ROOT, 'nf', nf)
            ma.setParameter(ROOT, 'max refinement', 4-nf)
            ma.setParameter(ROOT, 'npx', NPX)
            ma.setParameter(ROOT, 'npy', NPY)
            ma.setParameter(ROOT, 'npz', NPZ)
            ma.setParameter(ROOT, 'npf', NPF)
            TREE.write('parameter3D.xml')
            nptot = NPX*NPY*NPZ*NPF
            memtot = int(1024.*max(8/nptot, 2))
            print()
            print(CASE_PATH)
            EXE_STRING = pp.exe_pre(nptot, ' -N -W 12:00 -R "rusage[mem=' +
                                    str(memtot) + ']" ') + \
                pp.EXE_PATH + '/'+EXE
            print(EXE_STRING)
            os.system(EXE_STRING)
