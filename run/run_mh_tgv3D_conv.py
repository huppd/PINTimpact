""" runner for mh TGV conv """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterTGV3D.xml')
ROOT = TREE.getroot()

RE = 10.

ma.set_parameter(ROOT, 'withoutput', 1)
ma.set_parameter(ROOT, 'initial guess', 'zero')
# ma.set_parameter( ROOT, 'refinement level', 1 )

ma.set_parameter(ROOT, 'Re', RE)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*10

RUNS = range(1)

STS = [1., 10., 0.1]
# STS  = [ 1. ]
# STS  = [ 10., 0.1 ]


CASE_PATH[0] = pp.DATA_PATH + '/MHTGV3D_conv'
pp.mkdir(CASE_PATH, 0)

for st in STS:
    CASE_PATH[1] = '/a2_'+str(st)
    pp.mkdir(CASE_PATH, 1)
    #
    pp.chdir(CASE_PATH, 1)
    #
    ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*RE)
    # ma.set_parameter( ROOT, 'nx', 64*+1 )
    # ma.set_parameter( ROOT, 'ny', 64*+1 )
    # ma.set_parameter( ROOT, 'nz', 5 )
    # ma.set_parameter( ROOT, 'nf', nf )
    ma.set_parameter(ROOT, 'npx', 1)
    ma.set_parameter(ROOT, 'npy', 1)
    ma.set_parameter(ROOT, 'npz', 1)
    ma.set_parameter(ROOT, 'npf', 1)
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
                run) + pp.EXE_PATH+'/'+EXE
        print(exe_str)
        os.system(exe_str)
