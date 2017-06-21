import os
from math import pi
from platform_paths import *
import xml.etree.ElementTree as ET
import manipulator as ma

tree = ET.parse('../XML/parameter3D.xml')

root = tree.getroot()

exe = 'peri_navier3D'


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')



CASE_PATH = ['']*3

lambdas = [4]

ma.setParameter(root, 'alpha2', 2*pi*0.2*200)
ma.setParameter(root, 'Tolerance', 1.e-9)

npx = 8
npy = 2
npz = 4
ma.setParameter(root, 'npx', npx)
ma.setParameter(root, 'npy', npy)
ma.setParameter(root, 'npz', npz)

CASE_PATH[0] = '/case_study'
if not os.path.exists(DATA_PATH+CASE_PATH[0]):
    os.mkdir(DATA_PATH+CASE_PATH[0])
for lam in lambdas:
    CASE_PATH[1] = '/n_'+str(lam)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.system(' rm ./* -r -v  ')
    #ma.setParameter(root, 'refinement', lam)
    ma.setParameter(root, 'nx', 256+1)
    ma.setParameter(root, 'ny',  64+1)
    ma.setParameter(root, 'nz', 128+1)
    ma.setParameter(root, 'nf',  1)
    tree.write('parameter3D.xml')
    os.system(exe_pre(npx*npy*npz, ' -W 21:00 ') + EXE_PATH+exe + ' > output ')
