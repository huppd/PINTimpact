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


ST = 1./30.
RE = 300.

DX = 1
DY = 4
DZ = 1
#
NPX = 1
NPY = 2
NPZ = 4
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
NX = 65
NY = (1025-1)/2/DY+1
# NZ = (385-1)/2/DZ+1
NZ = 129
NZ = 129
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
#
ma.set_parameter(ROOT, 'Re', RE)
ma.set_parameter(ROOT, 'alpha2', 2.*pi*ST*RE)
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
ma.set_parameter(ROOT, 'Maximum Iterations', 20)
# ma.set_parameter(ROOT, 'Convergence Tolerance', 1.e-6)
ma.set_parameter(ROOT, 'Output Frequency', 1)


PRECS = [1]
# PRECS = [6, 7]
# PRECS = [2, 6, 7]
# PRECS = [5, 3, 6, 7]


# OMEGAS = [0.5, 1., 1.5]
# OMEGAS = [0.1, 0.2, 0.3, 0.4, 0.5]
OMEGAS = [0.1, 0.3, 0.5]
# OMEGAS = [0.5]


CASE_PATH = ['']*3


for side in ['right']:
    CASE_PATH[0] = pp.DATA_PATH + '/SHL_mode_prec'
    pp.mkdir(CASE_PATH, 0)
    for prec in PRECS:
        CASE_PATH[1] = '/prec_'+str(prec)
        pp.mkdir(CASE_PATH, 1)
        pp.chdir(CASE_PATH, 1)
        for omega in OMEGAS:
            CASE_PATH[2] = '/omega_'+str(omega)
            pp.mkdir(CASE_PATH, 2)
            pp.chdir(CASE_PATH, 2)
            #
            # ma.set_parameter(ROOT, 'preconditioner', 'none')
            ma.set_parameter(ROOT, 'type', prec)
            ma.set_parameter(ROOT, 'omega', omega)
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
