""" runner script to investigate mode preconditioner """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterSHL.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 1)

# make executable ready
EXE = 'modeConvDiff'
os.chdir(pp.EXE_PATH)
os.system('make ' + EXE + ' -j4')


st = 1./30.
re = 300.
# re = 100.

DX = 1
DY = 1
DZ = 1
#
NPX = 1
NPY = 4
NPZ = 1
NPF = 1
#
LXO = 22.5
LYO = 600.
LZO = 150.
#
NXO = 97
NYO = 1537
NZO = 513
#
NX = 97
NY = (1025-1)/4/DY+1
NZ = (385-1)/4/DZ+1
#
LX = LXO
LY = LYO*2./3./4./DY
LZ = round(LZO/(NZO-1)*(NZ-1), 1)
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
ma.set_parameter(ROOT, 'nf', 1)
ma.set_parameter(ROOT, 'npx', NPX)
ma.set_parameter(ROOT, 'npy', NPY)
ma.set_parameter(ROOT, 'npz', NPZ)
ma.set_parameter(ROOT, 'npf', 1)
ma.set_parameter(ROOT, 'Maximum Iterations', 40)
ma.set_parameter(ROOT, 'Convergence Tolerance', 1.e-6)
ma.set_parameter(ROOT, 'Output Frequency', 1)


PRECS = [1, 2, 3, 4, 5]
PRECS = [2, 3, 4, 5]
PRECS = [3, 4, 5]
# PRECS = [4, 5]
# PRECS = [4]

CYCLES = [2, 4, 8, 16]
CYCLES = [4, 8, 16]
# CYCLES = [4]
CYCLES = [8]

SWEEPS = [1, 2, 4, 8, 16]
SWEEPS = [2, 4, 8]
# SWEEPS = [1, 2]
SWEEPS = [8]

MAXGRIDS = [2, 3, 4]

CASE_PATH = ['']*6


# for side in ['left', 'right']:
# for side in ['left', 'none', 'right']:
# for side in ['left']:
# for side in ['none']:
for side in ['right']:
    # CASE_PATH[0] = pp.DATA_PATH + '/SHL_mode_prec_' + side
    CASE_PATH[0] = pp.DATA_PATH + '/SHL_mode_prec_none'
    pp.mkdir(CASE_PATH, 0)
    for prec in PRECS:
        CASE_PATH[1] = '/prec_'+str(prec)
        pp.mkdir(CASE_PATH, 1)
        pp.chdir(CASE_PATH, 1)
        for cycle in CYCLES:
            CASE_PATH[2] = '/cycle_'+str(cycle)
            pp.mkdir(CASE_PATH, 2)
            pp.chdir(CASE_PATH, 2)
            for sweep in SWEEPS:
                CASE_PATH[3] = '/sweep_'+str(sweep)
                pp.mkdir(CASE_PATH, 3)
                pp.chdir(CASE_PATH, 3)
                for max_grids in MAXGRIDS:
                    CASE_PATH[5] = '/maxGrids_'+str(max_grids)
                    pp.mkdir(CASE_PATH, 5)
                    pp.chdir(CASE_PATH, 5)
                    #
                    ma.set_parameter(ROOT, 'numGrids', max_grids)
                    # ma.set_parameter(ROOT, 'numGrids', 3)
                    ma.set_parameter(ROOT, 'type', prec)
                    # ma.set_parameter(ROOT, 'type', 5)
                    # ma.set_parameter(ROOT, 'cycle type', prec-4)
                    ma.set_parameter(ROOT, 'preconditioner', side)
                    ma.set_parameter(ROOT, 'numCycles', cycle)
                    ma.set_parameter(ROOT, 'numIters', sweep)
                    TREE.write('parameter3D.xml')
                    nptot = NPX*NPY*NPZ
                    print()
                    print(CASE_PATH)
                    exe_str = \
                        pp.exe_pre(nptot,
                                   ' -N -W 1:00 -R "rusage[mem=' +
                                   str(1024*2) + ']" ') + \
                        pp.EXE_PATH+'/'+EXE+' --realCase=1 '
                    print(exe_str)
                    os.system(exe_str)
