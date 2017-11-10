""" runner script to investigate mode preconditioner for Rayleigh streaming """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterStreaming2D.xml')
ROOT = TREE.getroot()

# ma.set_parameter(ROOT, 'withoutput', 0)

# make executable ready
EXE = 'modeConvDiff'
os.chdir(pp.EXE_PATH)
os.system('make ' + EXE + ' -j4')


st = 0.001
#
NPX = 2
NPY = 2
NPZ = 1
NPF = 1
#
#
#
ma.set_parameter(ROOT, 'Re', 100.)
ma.set_parameter(ROOT, 'alpha2', 2.*pi*0.01*100.)

ma.set_parameter(ROOT, 'lx', 2.)
ma.set_parameter(ROOT, 'ly', 2.)

ma.set_parameter(ROOT, 'npx', 2)
ma.set_parameter(ROOT, 'npy', 2)

ma.set_parameter(ROOT, 'nf', 1)


NXS = [33, 65, 129]

PRECS = [1, 2, 3, 4]
# PRECS = [2, 3, 4]
# PRECS = [3, 4]
# PRECS = [4]

# CYCLES = [1, 2, 4, 8, 16]
# CYCLES = [2, 3, 4, 6, 8]
CYCLES = [1, 2, 4]

# SWEEPS = [1, 2, 4, 8, 16]
SWEEPS = [1, 2, 4]

MAXGRIDS = [1, 3, 5]
# MAXGRIDS = [1, 2, 3]

CASE_PATH = ['']*6


# for side in ['left', 'right']:
for side in ['right']:
    CASE_PATH[0] = pp.DATA_PATH + '/stream_mode_prec_' + side
    pp.mkdir(CASE_PATH, 0)
    for nx in NXS:
        CASE_PATH[1] = '/nx_'+str(nx)
        pp.mkdir(CASE_PATH, 1)
        pp.chdir(CASE_PATH, 1)
        for prec in PRECS:
            CASE_PATH[2] = '/prec_'+str(prec)
            pp.mkdir(CASE_PATH, 2)
            pp.chdir(CASE_PATH, 2)
            for cycle in CYCLES:
                CASE_PATH[3] = '/cycle_'+str(cycle)
                pp.mkdir(CASE_PATH, 3)
                pp.chdir(CASE_PATH, 3)
                for sweep in SWEEPS:
                    CASE_PATH[4] = '/sweep_'+str(sweep)
                    pp.mkdir(CASE_PATH, 4)
                    pp.chdir(CASE_PATH, 4)
                    for max_grids in MAXGRIDS:
                        CASE_PATH[5] = '/maxGrids_'+str(max_grids)
                        pp.mkdir(CASE_PATH, 5)
                        pp.chdir(CASE_PATH, 5)
                        #
                        ma.set_parameter(ROOT, 'nx', nx)
                        ma.set_parameter(ROOT, 'ny', nx)
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
