""" running converferce for finite differences and Taylor-Green vortex """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3DTime.xml')
ROOT = TREE.getroot()

ma.setParameter(ROOT, 'withoutput', 1)
ma.setParameter(ROOT, 'initial guess', 'zero')
# ma.setParameter( ROOT, 'refinement level', 1 )

# make executable ready
EXE = 'peri_navier3DTime'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*4

RUNS = range(1)

RES = [10]
A2S = [0.1, 10., 1.]

NFS = [72]

ma.setParameter(ROOT, 'nx', 65)
ma.setParameter(ROOT, 'ny', 65)
ma.setParameter(ROOT, 'nz', 5)

CASE_PATH[0] = '/FDTGV_conv2'
pp.mkdir(CASE_PATH, 0)

for re in RES:
    CASE_PATH[1] = '/re_'+str(re)
    pp.mkdir(CASE_PATH, 1)
    for a2 in A2S:
        CASE_PATH[2] = '/a2_'+str(a2)
        pp.mkdir(CASE_PATH, 2)
        for nf in NFS:
            CASE_PATH[3] = '/nt_'+str(nf)
            pp.mkdir(CASE_PATH, 3)
            #
            pp.chdir(CASE_PATH, 3)
            #
            ma.setParameter(ROOT, 'Re', re)
            ma.setParameter(ROOT, 'alpha2', 2.*pi*a2*re)
            ma.setParameter(ROOT, 'nf', nf)
            ma.setParameter(ROOT, 'npx', 1)
            ma.setParameter(ROOT, 'npy', 1)
            ma.setParameter(ROOT, 'npz', 1)
            ma.setParameter(ROOT, 'npf', 12)
            TREE.write('parameter3D.xml')
            # nptot = npx[i]*npy[i]*npf[i]
            nptot = 12
            mem = int(max(1024, 60*1024/nptot))
            for run in RUNS:
                print()
                print(CASE_PATH)
                exeString = \
                    pp.exe_pre(nptot, ' -N  -R "rusage[mem=' +
                               str(mem) + ']" -W 6:00', run) + \
                    pp.EXE_PATH+'/'+EXE
                print(exeString)
                os.system(exeString)
