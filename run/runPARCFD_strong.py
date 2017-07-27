import os
import xml.etree.ElementTree as ET
from pylab import pi

import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3D.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 0)
ma.set_parameter(ROOT, 'refinement level', 1)
# ma.set_parameter(ROOT, 'Verbosity', 0)


# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

RUNS = range(1)


CASE_PATH = ['']*5


nps = [4, 8, 12]
nxs = [1, 2, 3, 4]
nfs = [1, 2, 3, 4]

nps = [4]
nxs = [1]
nfs = [1]

re = 200
st = 0.2


CASE_PATH[0] = '/strong'
if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])
for nf in nfs:
    CASE_PATH[1] = '/nf_'+str(nf)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    for n in nxs:
        CASE_PATH[2] = '/n_'+str(n)
        if not os.path.exists(pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1] +
                              CASE_PATH[2]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        for np in nps:
            CASE_PATH[3] = '/np_'+str(np)
            if not os.path.exists(pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1] +
                                  CASE_PATH[2] + CASE_PATH[3]):
                os.mkdir(pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1] +
                         CASE_PATH[2] + CASE_PATH[3])
            for npf in range(1, nf+2):
                CASE_PATH[4] = '/npf_'+str(npf)
                if not os.path.exists(pp.DATA_PATH + CASE_PATH[0] +
                                      CASE_PATH[1] + CASE_PATH[2] +
                                      CASE_PATH[3] + CASE_PATH[4]):
                    os.mkdir(pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1] +
                             CASE_PATH[2] + CASE_PATH[3] + CASE_PATH[4])
                os.chdir(pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1] +
                         CASE_PATH[2] + CASE_PATH[3] + CASE_PATH[4])
                os.system(' rm ./* -r -v  ')
                ma.set_parameter(ROOT, 'Re', re)
                ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*re)
                ma.set_parameter(ROOT, 'nx', 64*n+1)
                ma.set_parameter(ROOT, 'ny', 16*n+1)
                ma.set_parameter(ROOT, 'nz', 32*n+1)
                ma.set_parameter(ROOT, 'nf', nf)
                ma.set_parameter(ROOT, 'npx', np)
                ma.set_parameter(ROOT, 'npy', max(np/4, 1))
                ma.set_parameter(ROOT, 'npz', max(np/2, 1))
                ma.set_parameter(ROOT, 'npf', npf)
                TREE.write('parameter3D.xml')
                for run in RUNS:
                    exe_str = pp.exe_pre(npf*np*max(np/2, 1)*max(np/4, 1),
                                         ' -R "select[model==Opteron6174"] ' +
                                         ' -R "rusage[mem=8192]" -W 4:00',
                                         run) + pp.EXE_PATH + EXE
                    print(exe_str)
                    os.system(exe_str)
