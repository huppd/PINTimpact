""" runner for mh TGV conv """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterTGV.xml')
ROOT = TREE.getroot()

RE = 10.

ma.setParameter(ROOT, 'withoutput', 1)
ma.setParameter(ROOT, 'initial guess', 'zero')
# ma.setParameter( ROOT, 'refinement level', 1 )

ma.setParameter(ROOT, 'Re', RE)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.exe_path)
os.system('make '+EXE+' -j4')


CASE_PATH = ['', '', '', '', '', '', '', '', '']

RUNS = range(1)

A2S = [1., 10., 0.1]
# A2S  = [ 1. ]
# A2S  = [ 10., 0.1 ]


CASE_PATH[0] = pp.DATA_PATH + '/MHTGV_conv'
pp.mkdir(CASE_PATH, 0)

for a2 in A2S:
    CASE_PATH[1] = '/a2_'+str(a2)
    pp.mkdir(CASE_PATH, 1)
    #
    pp.chdir(CASE_PATH, 1)
    #
    ma.setParameter(ROOT, 'alpha2', 2.*pi*a2*RE)
    # ma.setParameter( ROOT, 'nx', 64*+1 )
    # ma.setParameter( ROOT, 'ny', 64*+1 )
    # ma.setParameter( ROOT, 'nz', 5 )
    # ma.setParameter( ROOT, 'nf', nf )
    ma.setParameter(ROOT, 'npx', 1)
    ma.setParameter(ROOT, 'npy', 1)
    ma.setParameter(ROOT, 'npz', 1)
    ma.setParameter(ROOT, 'npf', 1)
    TREE.write('parameter3D.xml')
    nptot = 1
    mem = int(max(1024, 16*1024/nptot))
    for run in RUNS:
        print()
        print(CASE_PATH)
        exe_str = \
            pp.exe_pre(
                nptot,
                ' -N -R beta -R "span[ptile=4]" -R "rusage[mem=' +
                str(mem) + ']" ',
                run) + pp.exe_path+'/'+EXE
        print(exe_str)
        os.system(exe_str)
