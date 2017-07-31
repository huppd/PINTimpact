import os
from math import pi
import xml.etree.ElementTree as ET

import platform_paths as pp
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3D.xml')
ROOT = TREE.getroot()

ma.set_parameter(ROOT, 'withoutput', 1)
npx = 8
npy = 2
npz = 2
ma.set_parameter(ROOT, 'npx', npx)
ma.set_parameter(ROOT, 'npy', npy)
ma.set_parameter(ROOT, 'npz', npz)

ma.set_parameter(ROOT, 'Convergence Tolerance', 1.e-1)

EXE = 'peri_navier3D'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*4

refinement = [12, 1]
ns = [1, 7]

ma.set_parameter(ROOT, 'alpha2', 2*pi*0.2*200)

CASE_PATH[0] = '/case_study_small'
if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])
for i in range(0, 1):
    CASE_PATH[1] = '/ref_'+str(i)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.system(' rm ./* -r -v  ')
    ma.set_parameter(ROOT, 'refinement level', refinement[i])
    ma.set_parameter(ROOT, 'nx', 256+1)
    ma.set_parameter(ROOT, 'ny', 64+1)
    ma.set_parameter(ROOT, 'nz', 128+1)
    ma.set_parameter(ROOT, 'nf', ns[i])
    TREE.write('parameter3D.xml')
    os.system(pp.exe_pre(npx*npy*npz, ' -W 48:00 -N -R "rusage[mem=8192]"') +
              pp.EXE_PATH+EXE + ' > output ')
