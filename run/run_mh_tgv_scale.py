""" runner for multiharmoic TGV scaling """
from math import pi
import os
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameterTGV.xml')
ROOT = TREE.getroot()

ma.setParameter(ROOT, 'withoutput', 0)
ma.setParameter(ROOT, 'initial guess', 'disturbed')
# ma.setParameter( ROOT, 'refinement level', 1 )

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*10

RUNS = range(10, 20)

A2S = [0.1, 1., 10.]

RE = 10

NPX = [1, 2, 2, 4, 4, 4, 4]
NPY = [1, 1, 2, 2, 4, 4, 4]
NPF = [1, 1, 1, 1, 1, 2, 3]

NXS = [65, 129, 257]
# NXS = [ 65 ]


for nx in NXS:
    CASE_PATH[0] = '/MHTGV_'+str(nx)
    pp.mkdir(CASE_PATH, 0)
    #
    for a2 in A2S:
        CASE_PATH[1] = '/a2_'+str(a2)
        pp.mkdir(CASE_PATH, 1)
        #
        for i in enumerate(NPX):
            CASE_PATH[2] = '/np_'+str(i)
            pp.mkdir(CASE_PATH, 2)
            #
            pp.chdir(CASE_PATH, 2)
            #
            ma.setParameter(ROOT, 'Re', RE)
            ma.setParameter(ROOT, 'alpha2', 2.*pi*a2*RE)
            ma.setParameter(ROOT, 'nx', nx)
            ma.setParameter(ROOT, 'ny', nx)
            ma.setParameter(ROOT, 'nz', 5)
            # ma.setParameter(ROOT, 'nf', nf)
            ma.setParameter(ROOT, 'npx', NPX[i])
            ma.setParameter(ROOT, 'npy', NPY[i])
            ma.setParameter(ROOT, 'npz', 1)
            ma.setParameter(ROOT, 'npf', NPF[i])
            TREE.write('parameter3D.xml')
            nptot = NPX[i]*NPY[i]*NPF[i]
            mem = int(max(1024, 16*1024/nptot))
            for run in RUNS:
                print()
                print CASE_PATH
                exe_str =  \
                    pp.exe_pre(
                        nptot,
                        ' -N -R beta -R "span[ptile=4]" -R "rusage[mem='
                        + str(mem) + ']" ', run) + pp.EXE_PATH+'/'+EXE
                print exe_str
                os.system(exe_str)
