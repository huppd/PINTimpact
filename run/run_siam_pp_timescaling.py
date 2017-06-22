from numpy import linspace
import numpy as np
from math import pi
from platformPaths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3D.xml')
root = tree.getroot()

ma.setParameter(root, 'withoutput', 0)
ma.setParameter(root, 'refinement level', 1)

# make executable ready
exe = 'peri_navier3D'
os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')


CASE_PATH = ['']*10

runs = range(10, 12)

re = 100
st = 0.1
nfs = 2**np.arange(1, 5) - 1
nfs = [15]
nx = 2
npx = 4


CASE_PATH[0] = '/scale_time_re_'+str(re)+'_a2_'+str(st)+'_nx_'+str(nx)+'_npx_'+str(npx)
mkdir(CASE_PATH, 0)

for nf in nfs:
    CASE_PATH[3] = '/nf_'+str(nf)
    mkdir(CASE_PATH, 3)
    for npf in 2**np.arange(4,5):
        if(npf<=(nf+1) and nf/npf<=15):
            CASE_PATH[6] = '/npf_'+str(npf)
            mkdir(CASE_PATH, 6)
            #
            chdir(CASE_PATH, 6)
            #
            ma.setParameter(root, 'Re', re)
            ma.setParameter(root, 'alpha2', 2.*pi*st*re)
            ma.setParameter(root, 'nx', 64*nx+1)
            ma.setParameter(root, 'ny', 16*nx+1)
            ma.setParameter(root, 'nz', 32*nx+1)
            ma.setParameter(root, 'nf', nf)
            ma.setParameter(root, 'npx',     npx)
            ma.setParameter(root, 'npy',  max(npx/2, 1))
            ma.setParameter(root, 'npz',     npx)
            ma.setParameter(root, 'npf',     npf)
            tree.write('parameter3D.xml')
            nptot = npx*npx*max(npx/2, 1)*npf
            for run in runs:
                print()
                print(CASE_PATH)
                # print(exe_pre(nptot, ' -N -R "select[model==Opteron6174"] -R "rusage[mem='+str(1024*1)+']" ', run) + EXE_PATH+'/'+exe)
                # os.system(exe_pre(nptot, ' -N -R "select[model==Opteron6174"] -R "rusage[mem='+str(1024*1)+']" ', run) + EXE_PATH+'/'+exe)
                print(exe_pre(nptot, ' -N -R "select[model==Opteron6174"]  ', run) + EXE_PATH+'/'+exe)
                os.system(exe_pre(nptot, ' -N -R "select[model==Opteron6174"]  ', run) + EXE_PATH+'/'+exe)
