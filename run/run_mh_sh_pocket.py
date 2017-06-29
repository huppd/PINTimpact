""" runner for swept Hiemenz flow """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterSHLabs.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 1)
ma.set_parameter(ROOT, 'refinement step', 1)
ma.set_parameter(ROOT, 'refinement tol', 1.e-4)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*5

RUNS = range(1)

RES = [300]
STS = [1./60., 1./30., 1./10.]
STS = [1./60., 1./10.]
STS = [1./60.]

NFS = [0, 1, 2]
NFS = [0]

NPX = 1
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

JACOBIAN = [True, False]
NITERS = [2, 4, 8]

CASE_PATH[0] = pp.DATA_PATH + '/ultimate'
pp.mkdir(CASE_PATH, 0)

for re in RES:
    for nf in NFS:
        for st in STS:
            for jacobi in JACOBIAN:
                CASE_PATH[1] = '/jacobi_'+str(jacobi)
                pp.mkdir(CASE_PATH, 1)
                for n_iters in NITERS:
                    CASE_PATH[2] = '/sweeps_'+str(n_iters)
                    pp.mkdir(CASE_PATH, 2)
                    pp.chdir(CASE_PATH, 2)
                    #
                    ma.set_parameter(ROOT, 'Re', re)
                    ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*re)
                    ma.set_parameter(ROOT, 'lx', LX)
                    ma.set_parameter(ROOT, 'ly', LY)
                    ma.set_parameter(ROOT, 'lz', LZ)
                    ma.set_parameter(ROOT, 'origin z', LZ/2.)
                    ma.set_parameter(ROOT, 'nx', NX)
                    ma.set_parameter(ROOT, 'ny', NY)
                    ma.set_parameter(ROOT, 'nz', NZ)
                    ma.set_parameter(ROOT, 'nf', nf)
                    ma.set_parameter(ROOT, 'max refinement', 4-nf)
                    ma.set_parameter(ROOT, 'npx', NPX)
                    ma.set_parameter(ROOT, 'npy', NPY)
                    ma.set_parameter(ROOT, 'npz', NPZ)
                    ma.set_parameter(ROOT, 'npf', NPF)
                    ma.set_parameter(ROOT, 'Jacobi', jacobi)
                    ma.set_insublist(ROOT, 'DivGrad', 'numIters', n_iters)
                    TREE.write('parameter3D.xml')
                    nptot = NPX*NPY*NPZ*NPF
                    memtot = int(1024.*max(8/nptot, 2))
                    print()
                    print(CASE_PATH)
                    EXE_STRING = pp.exe_pre(nptot, ' -N -W 12:00 ' +
                                            '-R "rusage[mem=' + str(memtot) +
                                            ']" ') + pp.EXE_PATH + '/'+EXE
                    print(EXE_STRING)
                    os.system(EXE_STRING)
