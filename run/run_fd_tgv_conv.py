""" running converferce for finite differences and Taylor-Green vortex """
import os
from math import pi
import xml.etree.ElementTree as ET
import platform_paths as pp
import manipulator as ma


# load parameter file
ma.set_ids('../XML/parameterTGVTime.xml')
TREE = ET.parse('../XML/parameterTGVTime.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 1)
ma.set_parameter(ROOT, 'initial guess', 'zero')
# ma.set_parameter( ROOT, 'refinement level', 1 )

# make executable ready
EXE = 'peri_navier3DTime'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*4

RUNS = range(1)

RES = [10]
STS = [0.1, 10., 1.]

NFS = [72]

ma.set_parameter(ROOT, 'nx', 65)
ma.set_parameter(ROOT, 'ny', 65)
ma.set_parameter(ROOT, 'nz', 5)

CASE_PATH[0] = pp.DATA_PATH + '/FDTGV_conv2'
pp.mkdir(CASE_PATH, 0)

for re in RES:
    CASE_PATH[1] = '/re_'+str(re)
    pp.mkdir(CASE_PATH, 1)
    for st in STS:
        CASE_PATH[2] = '/a2_'+str(st)
        pp.mkdir(CASE_PATH, 2)
        for nf in NFS:
            CASE_PATH[3] = '/nt_'+str(nf)
            pp.mkdir(CASE_PATH, 3)
            #
            pp.chdir(CASE_PATH, 3)
            #
            ma.set_parameter(ROOT, 'Re', re)
            ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*re)
            ma.set_parameter(ROOT, 'nf', nf)
            ma.set_parameter(ROOT, 'npx', 1)
            ma.set_parameter(ROOT, 'npy', 1)
            ma.set_parameter(ROOT, 'npz', 1)
            ma.set_parameter(ROOT, 'npf', 12)
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
