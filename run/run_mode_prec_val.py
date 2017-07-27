""" runner script to investigate mode preconditioner """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma
import numpy as np


# load parameter file
TREE = ET.parse('../XML/parameterModeConvDiff.xml')
ROOT = TREE.getroot()

# ma.set_parameter( ROOT, 'withoutput', 0 )
# ma.set_parameter( ROOT, 'refinement level', 1 )

# make executable ready
EXE = 'modeConvDiff'
os.chdir(pp.EXE_PATH)
os.system('make ' + EXE + ' -j4')


RUNS = range(1)

RES = 3.*10.**np.linspace(0, 2, 3)
ALPHA2S = 3.*10.**np.linspace(0, 3, 4)
# PRECS = [2, 3, 4]
PRECS = [4]
NXS = [3]


CASE_PATH = ['']*6


for side in ['left', 'right']:
    CASE_PATH[0] = pp.DATA_PATH + '/mode_prec_' + side
    pp.mkdir(CASE_PATH, 0)
    for i_re, re in enumerate(RES):
        CASE_PATH[1] = '/re_'+str(i_re)
        pp.mkdir(CASE_PATH, 1)
        for i_a2, alpha2 in enumerate(ALPHA2S):
            CASE_PATH[2] = '/a2_'+str(i_a2)
            pp.mkdir(CASE_PATH, 2)
            for prec in PRECS:
                CASE_PATH[3] = '/prec_'+str(prec)
                pp.mkdir(CASE_PATH, 3)
                for nx in NXS:
                    CASE_PATH[4] = '/nx_'+str(nx)
                    pp.mkdir(CASE_PATH, 4)
                    # for NPX in NPXs:
                    # CASE_PATH[5] = '/npx_'+str(NPX)
                    # mkdir( CASE_PATH, 5 )
                    # for NPF in range( 1, 2 ):
                    # # for NPF in range( 1, NF+2 ):
                    # # for NPF in range( NF+1, NF+2 ):
                    # CASE_PATH[6] = '/npf_'+str(NPF)
                    # mkdir( CASE_PATH, 6 )
                    #
                    pp.chdir(CASE_PATH, 4)
                    #
                    ma.set_parameter(ROOT, 'preconditioner', side)
                    ma.set_parameter(ROOT, 'type', prec)
                    ma.set_parameter(ROOT, 'Re', re)
                    ma.set_parameter(ROOT, 'alpha2', 2.*pi*alpha2)
                    ma.set_parameter(ROOT, 'nx', 16*nx+1)
                    ma.set_parameter(ROOT, 'ny', 32*nx+1)
                    ma.set_parameter(ROOT, 'nz', 6*nx+1)
                    ma.set_parameter(ROOT, 'nf', 1)
                    ma.set_parameter(ROOT, 'NPX', 1)
                    ma.set_parameter(ROOT, 'npy', 1)
                    ma.set_parameter(ROOT, 'npz', 1)
                    ma.set_parameter(ROOT, 'npf', 1)
                    TREE.write('parameter3D.xml')
                    nptot = 1
                    for run in RUNS:
                        print()
                        print(CASE_PATH)
                        exe_str = \
                            pp.exe_pre(nptot,
                                       ' -N -W 4:00 -R "rusage[mem=' +
                                       str(max(1024*4, 1024)) + ']" ', run) + \
                            pp.EXE_PATH+'/'+EXE
                        print(exe_str)
                        os.system(exe_str)
