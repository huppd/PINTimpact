import os
import numpy as np
from pylab import pi
from platform_paths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3D.xml')
root = tree.getroot()

ma.setParameter(root, 'withoutput', 0)
ma.setParameter(root, 'refinement level', 1)
# ma.setParameter(root, 'Verbosity', 0)


# make executable ready
exe = 'peri_navier3D'
os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')
exe = 'peri_navier3D'


runs = range(1)


CASE_PATH = ['']*5


nps = [4, 8, 12]
nxs = [1, 2, 3, 4]
nfs = [1, 2, 3, 4]

nps = [4]
nxs = [1]
nfs = [1]

re = 200
alpha2 = 0.2


CASE_PATH[0] = '/strong'
if not os.path.exists(DATA_PATH+CASE_PATH[0]):
    os.mkdir(DATA_PATH+CASE_PATH[0])
for nf in nfs:
    CASE_PATH[1] = '/nf_'+str(nf)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    for n in nxs:
        CASE_PATH[2] = '/n_'+str(n)
        if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        for np in nps:
            CASE_PATH[3] = '/np_'+str(np)
            if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]  ):
                os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
            for npf in range(1, nf+2):
                CASE_PATH[4] = '/npf_'+str(npf)
                if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4]):
                    os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                os.system(' rm ./* -r -v  ')
                ma.setParameter(root, 'Re', re)
                ma.setParameter(root, 'alpha2', 2.*pi*alpha2*re)
                ma.setParameter(root, 'nx', 64*n+1)
                ma.setParameter(root, 'ny', 16*n+1)
                ma.setParameter(root, 'nz', 32*n+1)
                ma.setParameter(root, 'nf', nf)
                ma.setParameter(root, 'npx',     np     )
                ma.setParameter(root, 'npy', max(np/4,1))
                ma.setParameter(root, 'npz', max(np/2,1))
                ma.setParameter(root, 'npf',     npf    )
                tree.write('parameter3D.xml')
                #os.system(exe_pre(npx*npy*npt,' -R lustre ')+EXE_PATH+exe+case_para+case_consts)
                for run in runs:
                    print(exe_pre(npf*np*max(np/2,1)*max(np/4,1),' -R "select[model==Opteron8380"] ')+EXE_PATH+exe)
                    # os.system(exe_pre(npf*np*max(np/2,1)*max(np/4,1),' -R "select[model==Opteron6174"] -R "rusage[mem=8192]" -W 24:00',run)+EXE_PATH+exe)
                    os.system(exe_pre(npf*np*max(np/2,1)*max(np/4,1),' -R "select[model==Opteron6174"] -R "rusage[mem=8192]" -W 4:00',run)+EXE_PATH+exe)
