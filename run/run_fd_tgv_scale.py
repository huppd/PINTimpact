""" runner for fd Taylor-Green vortex scaling """
import os
from math import pi
import xml.etree.ElementTree as ET

import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3DTime.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 0)
ma.set_parameter(ROOT, 'initial guess', 'exact')
# ma.set_parameter(ROOT, 'refinement level', 1)

# make executable ready
EXE = 'peri_navier3DTime'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*10

RUNS = range(10)

# res = [ 1,10,100]
STS = [0.1, 1., 10.]

RE = 10

# time first
NPX = [1, 1, 2, 2, 2, 4]
NPY = [1, 1, 1, 2, 2, 2]
NPF = [1, 3, 3, 3, 6, 6]
# space first
# NPX = [1, 2, 2, 2, 4, 4]
# NPY = [1, 1, 2, 2, 2, 4]
# NPF = [1, 1, 1, 3, 3, 3]

ma.set_parameter(ROOT, 'nx', 65)
ma.set_parameter(ROOT, 'ny', 65)
ma.set_parameter(ROOT, 'nz', 5)

CASE_PATH[0] = '/FD_scale'
pp.mkdir(CASE_PATH, 0)

for st in STS:
    CASE_PATH[1] = '/a2_'+str(st)
    pp.mkdir(CASE_PATH, 1)
    for i, npx in enumerate(NPX):
        CASE_PATH[2] = '/np_'+str(i)
        pp.mkdir(CASE_PATH, 2)
        #
        pp.chdir(CASE_PATH, 2)
        #
        ma.set_parameter(ROOT, 'Re', RE)
        ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*RE)
        ma.set_parameter(ROOT, 'nf', 72)
        ma.set_parameter(ROOT, 'npx', npx)
        ma.set_parameter(ROOT, 'npy', NPY[i])
        ma.set_parameter(ROOT, 'npz', 1)
        ma.set_parameter(ROOT, 'npf', NPF[i])
        TREE.write('parameter3D.xml')
        nptot = npx*NPY[i]*NPF[i]
        mem = int(max(1024, 29*1024/nptot))
        for run in RUNS:
            print()
            print(CASE_PATH)
            exeString = pp.exe_pre(nptot, ' -N -R beta -R "span[ptile=4]" ' +
                                   '-R "rusage[mem=' + str(mem) + ']" ', run) \
                + pp.EXE_PATH+'/'+EXE
            print(exeString)
            os.system(exeString)
