import os
from pylab import pi
from platform_paths import *


exe = 'peri_navier'


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')

CASE_PATH = ['']*5
npx = 4
npy = 4

case_consts = ' --linSolName="GCRODR"  --nf=8   --npx='+str(npx)+' --npy='+str(npy)+' --tolNOX=1.e-4  --tolBelos=1.e-5  --maxIter=10  --lx=2. --ly=2.  '

precTypes = [99]
ns = [4, 5, 6, 7]
res = [10, 100, 200]
STS = [10, 100, 200]
fixTypes = [1, 2, 4, 6, 9, 10]

#ns  = [6]
#precTypes = [0]
#res = [400]
res = [100]
#STS = [25**2]
STS = [125]
fixTypes = [2]

for precType in precTypes:
    CASE_PATH[0] = '/precType_'+str(precType)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]):
        os.mkdir(DATA_PATH+CASE_PATH[0])
    for n in ns:
        CASE_PATH[1] = '/n2_'+str(n)
        if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        for re in res:
            CASE_PATH[2] = '/re_'+str(re)
            if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            for st in STS:
                CASE_PATH[3] = '/alpha2_'+str(st)
                if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                    os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                for fixType in fixTypes:
                    CASE_PATH[4] = '/fixType_'+str(fixType)
                    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --nx='+str(2**n+1)+' --ny='+str(2**n+1)+' --re='+str(re)+' --alpha2='+str(st)+' --fixType='+str(fixType)+' '
                    print case_consts + case_para
                    os.system(exe_pre(npx*npy)+EXE_PATH+exe+case_para+case_consts)
