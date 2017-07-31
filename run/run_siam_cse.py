""" runner script used for SHL presented at SIAM CSE """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterSHLabs.xml')
ROOT = TREE.getroot()

# ma.set_parameter( ROOT, 'withoutput', 0 )
# ma.set_parameter( ROOT, 'refinement level', 1 )

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


RUNS = range(1)

RES = [300]
STS = [1./60., 1./30., 1./10.]
# STS  = [ 1./15.]
NF = 1
NX = 1
NPX = 1
NPF = 1


CASE_PATH = ['']*10

CASE_PATH[0] = '/ultimate'
pp.mkdir(CASE_PATH, 0)

for re in RES:
    CASE_PATH[1] = '/re_'+str(re)
    pp.mkdir(CASE_PATH, 1)
    for st in STS:
        CASE_PATH[2] = '/a2_'+str(int(1./st))
        pp.mkdir(CASE_PATH, 2)
        # for NF in nfs:
        # CASE_PATH[3] = '/nf_'+str(NF)
        # mkdir( CASE_PATH, 3 )
        # for NX in NXs:
        # CASE_PATH[4] = '/NX_'+str(NX)
        # mkdir( CASE_PATH, 4 )
        # for NPX in NPXs:
        # CASE_PATH[5] = '/npx_'+str(NPX)
        # mkdir( CASE_PATH, 5 )
        # for NPF in range( 1, 2 ):
        # # for NPF in range( 1, NF+2 ):
        # # for NPF in range( NF+1, NF+2 ):
        # CASE_PATH[6] = '/npf_'+str(NPF)
        # mkdir( CASE_PATH, 6 )
        #
        pp.chdir(CASE_PATH, 2)
        #
        ma.set_parameter(ROOT, 'Re', re)
        ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*re)
        ma.set_parameter(ROOT, 'nx', 48*NX+1)
        ma.set_parameter(ROOT, 'ny', 96*NX+1)
        ma.set_parameter(ROOT, 'nz', 64*NX+1)
        ma.set_parameter(ROOT, 'nf', NF)
        ma.set_parameter(ROOT, 'NPX', 1)
        ma.set_parameter(ROOT, 'npy', 1)
        ma.set_parameter(ROOT, 'npz', 1)
        ma.set_parameter(ROOT, 'npf', NPF)
        TREE.write('parameter3D.xml')
        nptot = 1*NPF
        for run in RUNS:
            print
            print CASE_PATH
            exe_str = \
                pp.exe_pre(nptot,
                           ' -N -W 8:00 -R "rusage[mem=' +
                           str(max(1024*4, 1024)) + ']" ', run) +  \
                pp.EXE_PATH+'/'+EXE
            print exe_str
            os.system(exe_str)
