""" runner script to investigate mode preconditioner """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterSHLabs.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 0)

# make executable ready
EXE = 'modeConvDiff'
os.chdir(pp.EXE_PATH)
os.system('make ' + EXE + ' -j4')


st = 1./30.
#
NPX = 1
NPY = 2
NPZ = 2
NPF = 1
#
LXO = 22.5
LYO = 600.
LZO = 150.
NXO = 97
NYO = 1537
NZO = 513
#
NX = 65
# NY = 193
NZ = 129
#
LX = round(1.2*LXO/(NXO-1)*(NX-1), 1)
LZ = round(1.2*LZO/(NZO-1)*(NZ-1), 1)
#
print('LX', LX)
# print('LY', LY)
print('LZ', LZ)
#
ma.set_parameter(ROOT, 'Re', 300.)
ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*300.)
ma.set_parameter(ROOT, 'lx', LX)
ma.set_parameter(ROOT, 'lz', LZ)
ma.set_parameter(ROOT, 'origin z', LZ/2.)
ma.set_parameter(ROOT, 'nx', NX)
# ma.set_parameter(ROOT, 'ny', NY)
ma.set_parameter(ROOT, 'nz', NZ)
ma.set_parameter(ROOT, 'nf', 1)
ma.set_parameter(ROOT, 'npx', NPX)
ma.set_parameter(ROOT, 'npy', NPY)
ma.set_parameter(ROOT, 'npz', NPZ)
ma.set_parameter(ROOT, 'npf', NPF)

NYS = [5]

PRECS = [1, 2, 3, 4]
PRECS = [2, 3, 4]
PRECS = [3, 4]
PRECS = [4]

CYCLES = [1, 2, 4, 8, 16]
CYCLES = [2, 3, 4, 6, 8]
CYCLES = [1, 2, 4]

SWEEPS = [1, 2, 4, 8, 16]
SWEEPS = [1, 2, 4]

MAXGRIDS = [1, 3, 5]
# MAXGRIDS = [1, 2, 3]

CASE_PATH = ['']*6


# for side in ['left', 'right']:
for side in ['left']:
    CASE_PATH[0] = pp.DATA_PATH + '/SHL_mode_prec3_' + side
    pp.mkdir(CASE_PATH, 0)
    for y in NYS:
        CASE_PATH[1] = '/ny_'+str(y)
        pp.mkdir(CASE_PATH, 1)
        pp.chdir(CASE_PATH, 1)
        for prec in PRECS:
            # CASE_PATH[1] = '/prec_'+str(prec)
            # pp.mkdir(CASE_PATH, 1)
            # pp.chdir(CASE_PATH, 1)
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
                        NY = y*64 + 1
                        LY = round(1.2*LYO/(NYO-1)*(NY-1), 1)
                        ma.set_parameter(ROOT, 'ly', LY)
                        ma.set_parameter(ROOT, 'maxGrids', max_grids)
                        ma.set_parameter(ROOT, 'type', prec)
                        ma.set_parameter(ROOT, 'preconditioner', side)
                        ma.set_parameter(ROOT, 'numCycles', cycle)
                        ma.set_parameter(ROOT, 'numIters', sweep)
                        TREE.write('parameter3D.xml')
                        nptot = NPX*NPY*NPZ*NPF
                        print()
                        print(CASE_PATH)
                        exe_str = \
                            pp.exe_pre(nptot,
                                       ' -N -W 1:00 -R "rusage[mem=' +
                                       str(max(1024*4, 1024)) + ']" ') + \
                            pp.EXE_PATH+'/'+EXE+' --realCase=1 '
                        print(exe_str)
                        os.system(exe_str)
