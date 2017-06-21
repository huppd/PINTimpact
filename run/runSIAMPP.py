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


CASE_PATH = ['']*7

runs = range(1)

res = [100]
alpha2s = [0.1]
nfs = [1]
nxs = [1, 2, 4]
npxs = [1, 2, 4]


CASE_PATH[0] = '/ultimate'
mkdir(CASE_PATH, 0)

for re in res:
    CASE_PATH[1] = '/re_'+str(re)
    mkdir(CASE_PATH, 1)
    for alpha2 in alpha2s:
        CASE_PATH[2] = '/a2_'+str(alpha2)
        mkdir(CASE_PATH, 2)
        for nf in nfs:
            CASE_PATH[3] = '/nf_'+str(nf)
            mkdir(CASE_PATH, 3)
            for nx in nxs:
                CASE_PATH[4] = '/nx_'+str(nx)
                mkdir(CASE_PATH, 4)
                for npx in npxs:
                    CASE_PATH[5] = '/npx_'+str(npx)
                    mkdir(CASE_PATH, 5)
                    for npf in range(1, nf+2):
                    # for npf in range(nf+1, nf+2):
                        CASE_PATH[6] = '/npf_'+str(npf)
                        mkdir(CASE_PATH, 6)
                        #
                        chdir(CASE_PATH, 6)
                        #
                        ma.setParameter(root, 'Re', re)
                        ma.setParameter(root, 'alpha2', 2.*pi*alpha2*re)
                        ma.setParameter(root, 'nx', 64*nx+1)
                        ma.setParameter(root, 'ny', 16*nx+1)
                        ma.setParameter(root, 'nz', 32*nx+1)
                        ma.setParameter(root, 'nf', nf)
                        ma.setParameter(root, 'npx',     npx)
                        ma.setParameter(root, 'npy', max(npx/4,1))
                        ma.setParameter(root, 'npz', max(npx/2,1))
                        ma.setParameter(root, 'npf',     npf)
                        tree.write('parameter3D.xml')
                        nptot = npx*max(npx/4,1)*max(npx/2,1)*npf
                        for run in runs:
                            print()
                            print(CASE_PATH)
                            print(exe_pre(nptot, ' -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(max(1024*1*nx*nf*npf/nptot,1024))+']" ', run) + EXE_PATH+'/'+exe)
                            os.system(exe_pre(nptot, ' -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(max(1024*1*nx*nf*npf/nptot,1024))+']" ', run) + EXE_PATH+'/'+exe)
