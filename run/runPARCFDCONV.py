import os
from math import pi
import xml.etree.ElementTree as ET

import platform_paths as pp
import manipulator as ma

TREE = ET.parse('../XML/parameter3D.xml')

ROOT = TREE.getroot()

EXE = 'peri_navier3D'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*3

lambdas = [4]

ma.setParameter(ROOT, 'alpha2', 2*pi*0.2*200)
ma.setParameter(ROOT, 'Tolerance', 1.e-9)

npx = 8
npy = 2
npz = 4
ma.setParameter(ROOT, 'npx', npx)
ma.setParameter(ROOT, 'npy', npy)
ma.setParameter(ROOT, 'npz', npz)

CASE_PATH[0] = '/case_study'
if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])
for lam in lambdas:
    CASE_PATH[1] = '/n_'+str(lam)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.system(' rm ./* -r -v  ')
    #ma.setParameter(ROOT, 'refinement', lam)
    ma.setParameter(ROOT, 'nx', 256+1)
    ma.setParameter(ROOT, 'ny', 64+1)
    ma.setParameter(ROOT, 'nz', 128+1)
    ma.setParameter(ROOT, 'nf', 1)
    TREE.write('parameter3D.xml')
    os.system(pp.exe_pre(npx*npy*npz, ' -W 21:00 ') + pp.EXE_PATH + EXE +
              ' > output ')
