import os
import numpy as np
from math import pi
import platform_paths as pp
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3D.xml')
ROOT = TREE.getroot()

ma.setParameter(ROOT, 'withoutput', 0)
ma.setParameter(ROOT, 'refinement level', 1)

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


RUNS = range(20, 25)

re = 100
st = 0.1
nf = 7
nx = 1
npx = 1


CASE_PATH = ['']*10

CASE_PATH[0] = '/time_re_' + str(re) + '_a2_' + str(st) + '_nf_' + str(nf) + \
    '_nx_' + str(nx) + '_npx_' + str(npx)

pp.mkdir(CASE_PATH, 0)

for npf in 2**np.arange(4):
    CASE_PATH[6] = '/npf_'+str(npf)
    pp.mkdir(CASE_PATH, 6)
    #
    pp.chdir(CASE_PATH, 6)
    #
    ma.setParameter(ROOT, 'Re', re)
    ma.setParameter(ROOT, 'alpha2', 2.*pi*st*re)
    ma.setParameter(ROOT, 'nx', 64*nx+1)
    ma.setParameter(ROOT, 'ny', 16*nx+1)
    ma.setParameter(ROOT, 'nz', 32*nx+1)
    ma.setParameter(ROOT, 'nf', nf)
    ma.setParameter(ROOT, 'npx', npx)
    ma.setParameter(ROOT, 'npy', max(npx/2, 1))
    ma.setParameter(ROOT, 'npz', npx)
    ma.setParameter(ROOT, 'npf', npf)
    TREE.write('parameter3D.xml')
    nptot = npx*npx*max(npx/2, 1)*npf
    for run in RUNS:
            print()
            print(CASE_PATH)
            print(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3)+']" ', run) + pp.EXE_PATH+'/'+EXE)
            os.system(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3)+']" ', run) + pp.EXE_PATH+'/'+EXE)

npx = [1, 2, 4, 4, 8, 8, 16, 16, 16]
npy = [1, 1, 1, 1, 1, 2, 2, 2, 4]
npz = [1, 1, 1, 2, 2, 4, 4, 8, 8]
npf = 1

nps = range(len(npx))

CASE_PATH = ['']*10

CASE_PATH[0] = '/space_re_'+str(re)+'_a2_'+str(st)+'_nf_'+str(nf)+'_nx_'+str(nx)+'_npf_'+str(npf)
pp.mkdir(CASE_PATH, 0)

for i in nps:
    nptot = npx[i]*npy[i]*npz[i]*npf
    CASE_PATH[5] = '/npx_'+str(nptot)
    pp.mkdir(CASE_PATH, 5)
    #
    pp.chdir(CASE_PATH, 6)
    #
    ma.setParameter(ROOT, 'Re', re)
    ma.setParameter(ROOT, 'alpha2', 2.*pi*st*re)
    ma.setParameter(ROOT, 'nx', 64*nx+1)
    ma.setParameter(ROOT, 'ny', 16*nx+1)
    ma.setParameter(ROOT, 'nz', 32*nx+1)
    ma.setParameter(ROOT, 'nf', nf)
    ma.setParameter(ROOT, 'npx', npx[i])
    ma.setParameter(ROOT, 'npy', npy[i])
    ma.setParameter(ROOT, 'npz', npz[i])
    ma.setParameter(ROOT, 'npf', npf)
    TREE.write('parameter3D.xml')
    for run in RUNS:
        print()
        print(CASE_PATH)
        print(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + pp.EXE_PATH+'/'+EXE)
        os.system(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + pp.EXE_PATH+'/'+EXE)

npx = [1, 2, 4, 4, 8, 8, 8, 8, 8]
npy = [1, 1, 1, 1, 1, 1, 1, 2, 2]
npz = [1, 1, 1, 2, 2, 2, 2, 2, 4]
npf = [1, 1, 1, 1, 1, 4, 8, 8, 8]

nps = range(len(npx))

CASE_PATH = ['']*10

CASE_PATH[0] = '/spacetime_re_'+str(re)+'_a2_'+str(st)+'_nf_'+str(nf)+'_nx_'+str(nx)
pp.mkdir(CASE_PATH, 0)

for i in nps:
    nptot = npx[i]*npy[i]*npz[i]*npf[i]
    CASE_PATH[5] = '/npxf_'+str(nptot)
    pp.mkdir(CASE_PATH, 5)
    #
    pp.chdir(CASE_PATH, 6)
    #
    ma.setParameter(ROOT, 'Re', re)
    ma.setParameter(ROOT, 'alpha2', 2.*pi*st*re)
    ma.setParameter(ROOT, 'nx', 64*nx+1)
    ma.setParameter(ROOT, 'ny', 16*nx+1)
    ma.setParameter(ROOT, 'nz', 32*nx+1)
    ma.setParameter(ROOT, 'nf', nf)
    ma.setParameter(ROOT, 'npx', npx[i])
    ma.setParameter(ROOT, 'npy', npy[i])
    ma.setParameter(ROOT, 'npz', npz[i])
    ma.setParameter(ROOT, 'npf', npf[i])
    TREE.write('parameter3D.xml')
    for run in RUNS:
        print()
        print(CASE_PATH)
        print(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + pp.EXE_PATH+'/'+EXE)
        os.system(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + pp.EXE_PATH+'/'+EXE)


npx = [1, 1, 1, 1, 2, 4, 8, 8, 8]
npy = [1, 1, 1, 1, 1, 1, 1, 2, 2]
npz = [1, 1, 1, 1, 1, 2, 2, 2, 4]
npf = [1, 2, 4, 8, 8, 8, 8, 8, 8]

nps = range(len(npx))

CASE_PATH = ['']*10
CASE_PATH[0] = '/timeSpace_re_' + str(re) + '_a2_' + str(st) + '_nf_' + \
    str(nf) + '_nx_'+str(nx)
pp.mkdir(CASE_PATH, 0)

for i in nps:
    nptot = npx[i]*npy[i]*npz[i]*npf[i]
    CASE_PATH[5] = '/npxf_'+str(nptot)
    pp.mkdir(CASE_PATH, 5)
    #
    pp.chdir(CASE_PATH, 6)
    #
    ma.setParameter(ROOT, 'Re', re)
    ma.setParameter(ROOT, 'alpha2', 2.*pi*st*re)
    ma.setParameter(ROOT, 'nx', 64*nx+1)
    ma.setParameter(ROOT, 'ny', 16*nx+1)
    ma.setParameter(ROOT, 'nz', 32*nx+1)
    ma.setParameter(ROOT, 'nf', nf)
    ma.setParameter(ROOT, 'npx', npx[i])
    ma.setParameter(ROOT, 'npy', npy[i])
    ma.setParameter(ROOT, 'npz', npz[i])
    ma.setParameter(ROOT, 'npf', npf[i])
    TREE.write('parameter3D.xml')
    for run in RUNS:
        print()
        print(CASE_PATH)
        print(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + pp.EXE_PATH+'/'+EXE)
        os.system(pp.exe_pre(nptot, ' -W 2:00 -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(1024*3/2)+']" ', run) + pp.EXE_PATH+'/'+EXE)
