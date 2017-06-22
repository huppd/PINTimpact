from math import pi
from platformPaths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3DTime.xml')
root = tree.getroot()

ma.setParameter(root, 'withoutput', 0)
ma.setParameter(root, 'initial guess', 'exact')
# ma.setParameter(root, 'refinement level', 1)

# make executable ready
exe = 'peri_navier3DTime'
os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')


CASE_PATH = ['']*10

runs = range(10)

# res = [ 1,10,100]
STS = [0.1, 1., 10.]

re = 10

# time first
npx = [1, 1, 2, 2, 2, 4]
npy = [1, 1, 1, 2, 2, 2]
npf = [1, 3, 3, 3, 6, 6]
# space first
# npx = [1, 2, 2, 2, 4, 4]
# npy = [1, 1, 2, 2, 2, 4]
# npf = [1, 1, 1, 3, 3, 3]

ma.setParameter(root, 'nx', 65)
ma.setParameter(root, 'ny', 65)
ma.setParameter(root, 'nz', 5)

CASE_PATH[0] = '/FD_scale'
mkdir(CASE_PATH, 0)

for st in STS:
    CASE_PATH[1] = '/a2_'+str(st)
    mkdir(CASE_PATH, 1)
    for i in range(len(npx)):
        CASE_PATH[2] = '/np_'+str(i)
        mkdir(CASE_PATH, 2)
        #
        chdir(CASE_PATH, 2)
        #
        ma.setParameter(root, 'Re', re)
        ma.setParameter(root, 'alpha2', 2.*pi*st*re)
        ma.setParameter(root, 'nf', 72)
        ma.setParameter(root, 'npx', npx[i])
        ma.setParameter(root, 'npy', npy[i])
        ma.setParameter(root, 'npz', 1  )
        ma.setParameter(root, 'npf', npf[i])
        tree.write('parameter3D.xml')
        nptot = npx[i]*npy[i]*npf[i]
        mem = int(max(1024, 29*1024/nptot))
        for run in runs:
            print()
            print(CASE_PATH)
            exeString = exe_pre(nptot, ' -N -R beta -R "span[ptile=4]" -R "rusage[mem=' +str(mem) + ']" ', run) + EXE_PATH+'/'+exe
            print(exeString )
            os.system(exeString)
