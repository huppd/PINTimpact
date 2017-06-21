import os
from pylab import pi
from platform_paths import *


exe = 'peri_navier'


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')

CASE_PATH = ['']*5

case_consts = ' --linSolName="TFQMR" --flow=1 --domain=1 --force=3 --radius=0.1 --rotation=1 --nfe=8   --npx=4 --npy=1 --tolNOX=1.e-6 --tolBelos=1.e-2  --maxIter=10  --lx=2. --ly=1.  '

precTypes = [0, 10]
ns = [4, 5, 6, 7]
res = [10, 100, 200]
alpha2s = [10, 100, 200]
fixTypes = [1, 2, 4, 6, 9, 10]

ns = [4]
res = [11]
alpha2s = [12]

for precType in precTypes:
    CASE_PATH[0] = 'precType_'+str(precType)
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
            for alpha2 in alpha2s:
                CASE_PATH[3] = '/alpha2_'+str(alpha2)
                if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                    os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                for fixType in fixTypes:
                    CASE_PATH[4] = '/fixType_'+str(fixType)
                    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --nx='+str(2*2**n+1)+' --ny='+str(2**n+1)+' --re='+str(re)+' --alpha2='+str(alpha2)+' --fixType='+str(fixType)+' '
                    print case_consts + case_para
                    os.system(exe_pre+EXE_PATH+exe+case_para+case_consts)
