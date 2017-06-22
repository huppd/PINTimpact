import os
from math import pi
from platform_paths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3D.xml')
root = tree.getroot()

ma.setParameter(root, 'withoutput', 0)
ma.setParameter(root, 'refinement tol', 1.e-6)
npx = 8
npy = 2
npz = 4
ma.setParameter(root, 'npx', npx)
ma.setParameter(root, 'npy', npy)
ma.setParameter(root, 'npz', npz)

# make executable ready
exe = 'peri_navier3D'
os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')


CASE_PATH = ['']*4

ns = [2]
res = [100, 100, 400]
STS = [0.05, 0.1, 0.2, 0.4]

res = [200]
STS = [0.05, 0.1]

CASE_PATH[0] = '/ultimate'
if not os.path.exists(DATA_PATH+CASE_PATH[0]):
    os.mkdir(DATA_PATH+CASE_PATH[0])

for n in ns:
    CASE_PATH[1] = ''  # '/n_'+str(n)
    #if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        #os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    for re in res:
        CASE_PATH[2] = '/re_'+str(re)
        if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        for st in STS:
            CASE_PATH[3] = '/a2_'+str(st)
            if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
            os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
            os.system(' rm ./* -r -v  ')
            #
            ma.setParameter(root, 'Re', re)
            ma.setParameter(root, 'alpha2', 2.*pi*st*re)
            ma.setParameter(root, 'nx', 128*2+1)
            ma.setParameter(root, 'ny',  32*2+1)
            ma.setParameter(root, 'nz',  64*2+1)
            tree.write('parameter3D.xml')
            # os.system(exe_pre(npx*npy*npz,' -W 48:00 ') + EXE_PATH+exe + ' > output ')
            print(exe_pre(npx*npy*npz,' -W 48:00 ') + EXE_PATH+exe + ' > output ')
