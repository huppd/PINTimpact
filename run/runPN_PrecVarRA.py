import os
from pylab import pi
from platform_paths import *


exe = 'peri_navier'


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')

CASE_PATH = ['']*5

npx = 1
npy = 1
npt = 1


case_consts = ' --npx='+str(npx)+' --npy='+str(npy)+' --npf='+str(npt)+' --tolNOX=1.e-6  --tolBelos=1.e-3 --tolInnerBelos=1.e-3 --maxIter=20  --lx=2. --ly=2.  --initZero=0 --numCycles=2 '#--domain=1 --linSolver=GCRODR '

precTypes = [0, 1, 2]
precTypes = [2]
ns = [6]
#ns = [ 5, 6]
res = [1, 25, 50, 75, 100]
STS = [1, 9, 25, 64, 100]

#precTypes = [2, 3]
#precTypes = [2]
#precTypes = [3]
#ns = [5]
#res       = [100]
#STS   = [100]


#ns  = [6, 7]

for st in STS:
    CASE_PATH[0] = '/a2_'+str(st)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]):
        os.mkdir(DATA_PATH+CASE_PATH[0])
    for re in res:
        CASE_PATH[1] = '/re_'+str(re)
        if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        for n in ns:
            CASE_PATH[2] = '/n2_'+str(n)
            if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            for precType in precTypes:
                CASE_PATH[3] = '/precType_'+str(precType)
                if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]  ):
                    os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                os.system(' rm ./* -r -v  ')
                case_para = ' --nx='+str(2**n+1)+' --ny='+str(2**n+1)+' --nf='+str(2*(n-1))+' --withprec='+str(precType)+'  --re='+str(re)+' --alpha2='+str(st)+' --maxGrids='+str(n-2)+' ' 
                #os.system(exe_pre(npx*npy*npt,' -R lustre ')+EXE_PATH+exe+case_para+case_consts)
                print(exe_pre(npx*npy*npt)+EXE_PATH+exe+case_para+case_consts +' > output ')
                os.system(exe_pre(npx*npy*npt)+EXE_PATH+exe+case_para+case_consts +' > output ')
