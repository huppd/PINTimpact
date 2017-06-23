import os
from math import pi
import platform_paths as pp
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3D.xml')
ROOT = TREE.getroot()

ma.setParameter(ROOT, 'withoutput', 0)
ma.setParameter(ROOT, 'refinement tol', 1.e-6)
npx = 8
npy = 2
npz = 4
ma.setParameter(ROOT, 'npx', npx)
ma.setParameter(ROOT, 'npy', npy)
ma.setParameter(ROOT, 'npz', npz)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*4

ns = [2]
res = [100, 100, 400]
STS = [0.05, 0.1, 0.2, 0.4]

res = [200]
STS = [0.05, 0.1]

CASE_PATH[0] = '/ultimate'
if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])

for n in ns:
    CASE_PATH[1] = ''  # '/n_'+str(n)
    #if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        #os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    for re in res:
        CASE_PATH[2] = '/re_'+str(re)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        for st in STS:
            CASE_PATH[3] = '/a2_'+str(st)
            if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
            os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
            os.system(' rm ./* -r -v  ')
            #
            ma.setParameter(ROOT, 'Re', re)
            ma.setParameter(ROOT, 'alpha2', 2.*pi*st*re)
            ma.setParameter(ROOT, 'nx', 128*2+1)
            ma.setParameter(ROOT, 'ny',  32*2+1)
            ma.setParameter(ROOT, 'nz',  64*2+1)
            TREE.write('parameter3D.xml')
            # os.system(pp.exe_pre(npx*npy*npz,' -W 48:00 ') + pp.EXE_PATH+EXE + ' > output ')
            print(pp.exe_pre(npx*npy*npz, ' -W 48:00 ') + pp.EXE_PATH+EXE + ' > output ')
