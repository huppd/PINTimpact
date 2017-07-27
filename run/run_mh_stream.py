""" runner for Rayleigh streaming """
from math import pi
import os
import xml.etree.ElementTree as ET
import manipulator as ma
import platform_paths as pp


# load parameter file
ma.set_ids('../XML/parameterStreaming2D.xml')
TREE = ET.parse('../XML/parameterStreaming2D.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'nx', 33)
ma.set_parameter(ROOT, 'ny', 33)

ma.set_parameter(ROOT, 'lx', 2.)
ma.set_parameter(ROOT, 'ly', 2.)


# case_consts = ' --nfe=17 --tolNOX=1.e-6 --tolBelos=1.e-4 ' + \
    # '--tolNF=1.e-6 --maxIter=10  --domain=2 '

RES = [1, 10, 100]
STS = [1, 10, 100, 1000]

# make executable ready
EXE = 'peri_navier2D'
os.chdir(pp.EXE_PATH)
os.system('make ' + EXE + ' -j4')

# work direcotries
CASE_PATH = ['']*3
CASE_PATH[0] = pp.DATA_PATH + '/streaming'
pp.mkdir(CASE_PATH, 0)

for st in STS:
    CASE_PATH[1] = '/alpha2_'+str(st)
    pp.mkdir(CASE_PATH, 1)
    for re in RES:
        CASE_PATH[2] = '/re_'+str(re)
        pp.mkdir(CASE_PATH, 2)
        pp.chdir(CASE_PATH, 2)
        #
        ma.set_parameter(ROOT, 'Re', re)
        ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*re)
        print(CASE_PATH)
        TREE.write('parameter3D.xml')
        nptot = 1
        memtot = int(1024.*max(16/nptot, 2))
        print()
        print(CASE_PATH)
        EXE_STRING = pp.exe_pre(nptot, ' -N -W 8:00 ' +
                                '-R "rusage[mem=' + str(memtot) +
                                ']" ') + pp.EXE_PATH + '/'+EXE
        print(EXE_STRING)
        os.system(EXE_STRING)
