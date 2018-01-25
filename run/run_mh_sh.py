""" runner for swept Hiemenz flow """
from math import pi
import os
import xml.etree.ElementTree as ET
import manipulator as ma
import platform_paths as pp


# load parameter file
ma.set_ids('../XML/parameterSHL.xml')
TREE = ET.parse('../XML/parameterSHL.xml')
ROOT = TREE.getroot()

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


ma.set_parameter(ROOT, 'withoutput', 1)
ma.set_parameter(ROOT, 'refinement step', 1)
ma.set_parameter(ROOT, 'refinement tol', 1.e-4)


STS = [1./60., 1./30., 1./10.]

st = 1./30.

NF = 2


CASE_PATH = ['']*3


NPX = 1
NPY = 2
NPZ = 4
NPF = 2
#
LXO = 22.5
LYO = 600.
LZO = 150.
#
NXO = 97
NYO = 1537
NZO = 513
#
NX = 65
NY = 65
NZ = 65
# NY = 129
# NZ = 129
#
LX = LXO
LX = round(LXO*1.5/(NXO-1)*(NX-1), 1)
LY = round(LYO*1.5/(NYO-1)*(NY-1), 1)
LZ = round(LZO*1.5/(NZO-1)*(NZ-1), 1)
#
print('NX', NX)
print('NY', NY)
print('NZ', NZ)
#
print('LX', LX)
print('LY', LY)
print('LZ', LZ)
#
print('DX', LX/LXO*(NXO-1)/(NX-1))
print('DY', LY/LYO*(NYO-1)/(NY-1))
print('DZ', LZ/LZO*(NZO-1)/(NZ-1))
#
ma.set_parameter(ROOT, 'Re', 300.)
ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*300.)
ma.set_parameter(ROOT, 'lx', LX)
ma.set_parameter(ROOT, 'ly', LY)
ma.set_parameter(ROOT, 'lz', LZ)
ma.set_parameter(ROOT, 'origin z', LZ/2.)
ma.set_parameter(ROOT, 'nx', NX)
ma.set_parameter(ROOT, 'ny', NY)
ma.set_parameter(ROOT, 'nz', NZ)
ma.set_parameter(ROOT, 'nf', NF)
# ma.set_parameter(ROOT, 'max refinement', 5-NF)
ma.set_parameter(ROOT, 'max refinement', 1)
ma.set_parameter(ROOT, 'npx', NPX)
ma.set_parameter(ROOT, 'npy', NPY)
ma.set_parameter(ROOT, 'npz', NPZ)
ma.set_parameter(ROOT, 'npf', NPF)
ma.set_parameter(ROOT, 'initial guess', 'base')
# ma.set_insublist(ROOT, 'Line Search', 'Method', 'Full Step')
# ma.set_parameter(ROOT, 'Convergence Tolerance', 0.01)
# ma.set_insublist(ROOT, 'Coarse Grid Solver', 'numIters', 4)
# ma.set_insublist(ROOT, 'Coarse Grid Solver', 'Jacobi',
        # True)


CASE_PATH[0] = pp.DATA_PATH + '/shbl_' + str(LY)
pp.mkdir(CASE_PATH, 0)
pp.chdir(CASE_PATH, 0)

TREE.write('parameter3D.xml')
nptot = int(NPX*NPY*NPZ*NPF)
memtot = int(1024.*2)
print()
print(CASE_PATH)
EXE_STRING = pp.exe_pre(nptot, ' -N -W 24:00 ' +
                        '-R "rusage[mem=' + str(memtot) +
                        ']" ') + pp.EXE_PATH + '/'+EXE
print(EXE_STRING)
os.system(EXE_STRING)
