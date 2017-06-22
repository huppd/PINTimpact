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



runs = range(20, 25)

re = 100
st = 0.1
nf = 7
nx = 1
npx = 1


CASE_PATH = ['']*10

CASE_PATH[0] = '/time_re_'+str(re)+'_a2_'+str(st)+'_nf_'+str(nf)+'_nx_'+str(nx)+'_npx_'+str(npx)
mkdir(CASE_PATH, 0)

for npf in 2**np.arange(4):
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
    ma.setParameter(root, 'npx',     npx    )
    ma.setParameter(root, 'npy',  max(npx/2,1))
    ma.setParameter(root, 'npz',     npx    )
    ma.setParameter(root, 'npf',     npf    )
    tree.write('parameter3D.xml')
    nptot = npx*npx*max(npx/2,1)*npf
    for run in runs:
            print()
            print(CASE_PATH)
            print(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3)+']" ', run) + EXE_PATH+'/'+exe)
            os.system(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3)+']" ', run) + EXE_PATH+'/'+exe)

npx = [1, 2, 4, 4, 8, 8, 16, 16, 16]
npy = [1, 1, 1, 1, 1, 2,  2,  2, 4]
npz = [1, 1, 1, 2, 2, 4,  4,  8, 8]
npf = 1

nps = range(len(npx))

CASE_PATH = ['']*10

CASE_PATH[0] = '/space_re_'+str(re)+'_a2_'+str(st)+'_nf_'+str(nf)+'_nx_'+str(nx)+'_npf_'+str(npf)
mkdir(CASE_PATH, 0)

for i in nps:
    nptot = npx[i]*npy[i]*npz[i]*npf
    CASE_PATH[5] = '/npx_'+str(nptot)
    mkdir(CASE_PATH, 5)
    #
    chdir(CASE_PATH, 6)
    #
    ma.setParameter(root, 'Re', re)
    ma.setParameter(root, 'alpha2', 2.*pi*st*re)
    ma.setParameter(root, 'nx', 64*nx+1)
    ma.setParameter(root, 'ny', 16*nx+1)
    ma.setParameter(root, 'nz', 32*nx+1)
    ma.setParameter(root, 'nf', nf)
    ma.setParameter(root, 'npx',  npx[i])
    ma.setParameter(root, 'npy',  npy[i])
    ma.setParameter(root, 'npz',  npz[i])
    ma.setParameter(root, 'npf',  npf  )
    tree.write('parameter3D.xml')
    for run in runs:
        print()
        print(CASE_PATH)
        print(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + EXE_PATH+'/'+exe)
        os.system(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + EXE_PATH+'/'+exe)

npx = [1, 2, 4, 4, 8, 8, 8, 8, 8]
npy = [1, 1, 1, 1, 1, 1, 1, 2, 2]
npz = [1, 1, 1, 2, 2, 2, 2, 2, 4]
npf = [1, 1, 1, 1, 1, 4, 8, 8, 8]

nps = range(len(npx))

CASE_PATH = ['']*10

CASE_PATH[0] = '/spacetime_re_'+str(re)+'_a2_'+str(st)+'_nf_'+str(nf)+'_nx_'+str(nx)
mkdir(CASE_PATH, 0)

for i in nps:
    nptot = npx[i]*npy[i]*npz[i]*npf[i]
    CASE_PATH[5] = '/npxf_'+str(nptot)
    mkdir(CASE_PATH, 5)
    #
    chdir(CASE_PATH, 6)
    #
    ma.setParameter(root, 'Re', re)
    ma.setParameter(root, 'alpha2', 2.*pi*st*re)
    ma.setParameter(root, 'nx', 64*nx+1)
    ma.setParameter(root, 'ny', 16*nx+1)
    ma.setParameter(root, 'nz', 32*nx+1)
    ma.setParameter(root, 'nf', nf)
    ma.setParameter(root, 'npx',  npx[i])
    ma.setParameter(root, 'npy',  npy[i])
    ma.setParameter(root, 'npz',  npz[i])
    ma.setParameter(root, 'npf',  npf[i])
    tree.write('parameter3D.xml')
    for run in runs:
        print()
        print(CASE_PATH)
        print(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + EXE_PATH+'/'+exe)
        os.system(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + EXE_PATH+'/'+exe)


npx = [1, 1, 1, 1, 2, 4, 8, 8, 8]
npy = [1, 1, 1, 1, 1, 1, 1, 2, 2]
npz = [1, 1, 1, 1, 1, 2, 2, 2, 4]
npf = [1, 2, 4, 8, 8, 8, 8, 8, 8]

nps  = range(len(npx))

CASE_PATH = ['']*10
CASE_PATH[0] = '/timeSpace_re_'+str(re)+'_a2_'+str(st)+'_nf_'+str(nf)+'_nx_'+str(nx)
mkdir(CASE_PATH, 0)

for i in nps:
    nptot = npx[i]*npy[i]*npz[i]*npf[i]
    CASE_PATH[5] = '/npxf_'+str(nptot)
    mkdir(CASE_PATH, 5)
    #
    chdir(CASE_PATH, 6)
    #
    ma.setParameter(root, 'Re', re)
    ma.setParameter(root, 'alpha2', 2.*pi*st*re)
    ma.setParameter(root, 'nx', 64*nx+1)
    ma.setParameter(root, 'ny', 16*nx+1)
    ma.setParameter(root, 'nz', 32*nx+1)
    ma.setParameter(root, 'nf', nf)
    ma.setParameter(root, 'npx',  npx[i])
    ma.setParameter(root, 'npy',  npy[i])
    ma.setParameter(root, 'npz',  npz[i])
    ma.setParameter(root, 'npf',  npf[i])
    tree.write('parameter3D.xml')
    for run in runs:
        print()
        print(CASE_PATH)
        print(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + EXE_PATH+'/'+exe)
        os.system(exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + EXE_PATH+'/'+exe)
