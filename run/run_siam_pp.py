""" runner for SIAM PP results """
import os
from math import pi
import xml.etree.ElementTree as ET

import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3D.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 0)
ma.set_parameter(ROOT, 'refinement level', 1)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*7

RUNS = range(1)

RES = [100]
STS = [0.1]
NFS = [1]
NXS = [1, 2, 4]
NPXS = [1, 2, 4]


CASE_PATH[0] = '/ultimate'
pp.mkdir(CASE_PATH, 0)

for re in RES:
    CASE_PATH[1] = '/re_'+str(re)
    pp.mkdir(CASE_PATH, 1)
    for st in STS:
        CASE_PATH[2] = '/a2_'+str(st)
        pp.mkdir(CASE_PATH, 2)
        for nf in NFS:
            CASE_PATH[3] = '/nf_'+str(nf)
            pp.mkdir(CASE_PATH, 3)
            for nx in NXS:
                CASE_PATH[4] = '/nx_'+str(nx)
                pp.mkdir(CASE_PATH, 4)
                for npx in NPXS:
                    CASE_PATH[5] = '/npx_'+str(npx)
                    pp.mkdir(CASE_PATH, 5)
                    for npf in range(1, nf+2):
                    # for npf in range(nf+1, nf+2):
                        CASE_PATH[6] = '/npf_'+str(npf)
                        pp.mkdir(CASE_PATH, 6)
                        #
                        pp.chdir(CASE_PATH, 6)
                        #
                        ma.set_parameter(ROOT, 'Re', re)
                        ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*re)
                        ma.set_parameter(ROOT, 'nx', 64*nx+1)
                        ma.set_parameter(ROOT, 'ny', 16*nx+1)
                        ma.set_parameter(ROOT, 'nz', 32*nx+1)
                        ma.set_parameter(ROOT, 'nf', nf)
                        ma.set_parameter(ROOT, 'npx', npx)
                        ma.set_parameter(ROOT, 'npy', max(npx/4, 1))
                        ma.set_parameter(ROOT, 'npz', max(npx/2, 1))
                        ma.set_parameter(ROOT, 'npf', npf)
                        TREE.write('parameter3D.xml')
                        nptot = npx*max(npx/4, 1)*max(npx/2, 1)*npf
                        for run in RUNS:
                            print()
                            print(CASE_PATH)
                            exe_str = pp.exe_pre(nptot, ' -N -R "select[' +
                                                 'model==Opteron8380"] ' +
                                                 '-R "rusage[mem=' +
                                                 str(max(1024*nx*nf*npf/nptot,
                                                         1024))+']" ',
                                                 run) + pp.EXE_PATH+'/' + EXE
                            print(exe_str)
                            os.system(exe_str)
