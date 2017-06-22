import os
from pylab import log2
from platform_paths import *


exe = 'peri_navier3D'

runs = 4

os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')


CASE_PATH = ['']*5


case_consts = [' --baseflow=1 --force=1  --domain=0 --nf=4  --tolNOX=1.e-6  --tolInnerBelos=1.e-2 --maxIter=2 --ly=2. ']

nxs = [1, 2, 3]
nps = [1, 2, 3]
#nxs = [2]
#nps = [1]
re = 100
st = 576


CASE_PATH[0] = '/weak'
if not os.path.exists(DATA_PATH+CASE_PATH[0]):
    os.mkdir(DATA_PATH+CASE_PATH[0])
#os.chdir(DATA_PATH+CASE_PATH[0])
#os.system(' rm ./* -r -v  ')
#for st in sts:
    #CASE_PATH[1] = '/a2_'+str(st)
    #if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        #os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    #for re in res:
        #CASE_PATH[2] = '/re_'+str(re)
        #if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
            #os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
for n in ns:
    CASE_PATH[3] = '/n2_'+str(n)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]  ):
        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
    os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
    case_para = ' --nx='+str(192*n+1)+' --ny='+str(48+1)+' --nz='+str(96*n+1)+' --npx='+str(n)+' --npy='+str(2)+' --npz='+str(n)+' --maxGrids='+str(int(3+log2(1)))+' --re='+str(re)+' --alpha2='+str(st)+' --lx='+str(8*n)+'. --lz='+str(4*n)+'. '
    #os.system(exe_pre(npx*npy*npt,' -R lustre ')+EXE_PATH+exe+case_para+case_consts)
    print(exe_pre(2*n**2,' -R "select[model==Opteron6174"] ')+EXE_PATH+exe+case_para+case_consts[0] +' > output ')
    for run in range(runs):
        os.system(exe_pre(2*n**2,'',run)+EXE_PATH+exe+case_para+case_consts[0] +'  ')
