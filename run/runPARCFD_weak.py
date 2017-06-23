import os
from pylab import pi, log2
import platform_paths as pp
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
TREE = ET.parse('../XML/parameter3D.xml')
ROOT = TREE.getroot()

ma.setParameter(ROOT, 'withoutput', 0)


# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

runs = range(0, 4)


CASE_PATH = ['']*5


nxs = [1, 2, 3]
nps = [1, 2, 3]
#nxs = [2]
#nps = [1]
re = 200
st = 0.2


CASE_PATH[0] = '/weak'
if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])
#os.chdir(pp.DATA_PATH+CASE_PATH[0])
#os.system(' rm ./* -r -v  ')
#for st in sts:
    #CASE_PATH[1] = '/a2_'+str(st)
    #if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        #os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
for nx in nxs:
    CASE_PATH[2] = '/nx_'+str(nx)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
    for np in nps:
        CASE_PATH[3] = '/np_'+str(np)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]  ):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
        os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
        os.system(' rm ./* -r -v  ')
        case_para = ' --nx='+str(128*np*nx/2+1)+' --ny='+str(32*np*nx/2+1)+' --nz='+str(64*np*nx/2+1)+' --npx='+str(4*np)+' --npy='+str(np)+' --npz='+str(np*2)+' --maxGrids='+str(int(3+log2(np*nx)))+' --re='+str(re)+' --alpha2='+str(st)+' --lx='+str(8)+'. --lz='+str(4)+'. '
        ma.setParameter(ROOT, 'Re', re)
        ma.setParameter(ROOT, 'alpha2', 2.*pi*st*re)
        ma.setParameter(ROOT, 'nx', 128*nx*np+1)
        ma.setParameter(ROOT, 'ny',  32*nx*np/2+1)
        ma.setParameter(ROOT, 'nz',  64*nx*np+1)
        ma.setParameter(ROOT, 'npx', 4*np)
        ma.setParameter(ROOT, 'npy',   np)
        ma.setParameter(ROOT, 'npz', 2*np)
        TREE.write('parameter3D.xml')
        for run in runs:
            print(pp.exe_pre(8*np*np*np, ' -R "select[model==Opteron6174"] ')+pp.EXE_PATH+EXE+case_para +' > output ')
            os.system(pp.exe_pre(8*np*np*np, ' -R "select[model==Opteron6174"] -W 4:00 ', run)+pp.EXE_PATH+EXE)
