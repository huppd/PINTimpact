import os
from math import pi
from platform_paths import *


exe = 'peri_navier4'


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')

CASE_PATH = ['']*5

npx = 4
npy = 1
npt = 8

case_consts = ' --linSolName="GCRODR" --flow=1 --domain=1 --force=1   --npx='+str(npx)+' --npy='+str(npy)+' --npt='+str(npt)+' --tolNOX=1.e-6  --tolBelos=1.e-2  --maxIter=40  --lx=8. --ly=2. --amp=0.2 --radius=0.2 --xm='+str(1./8.)+' '

precTypes = [0]
#ns = [4, 5, 6]
ns = [5, 6, 7]
res = [10, 100, 200]
alpha2s = [10, 100, 200]
fixTypes = [1]


#ns  = [7]
#precTypes = [0,1,2]
res = [150]
alpha2s = [12]
#alpha2s = [2.*pi*0.1*200, 2.*pi*0.2*200, 2.*pi*0.3*200]
alpha2s = [2.*pi*0.2*res[0]]
fixTypes = [1]

for precType in precTypes:
    CASE_PATH[0] = '/convshow'+str(precType)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]):
        os.mkdir(DATA_PATH+CASE_PATH[0])
    for fixType in fixTypes:
        CASE_PATH[1] = '/fixType_'+str(fixType)
        if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        for re in res:
            CASE_PATH[2] = '/re_'+str(re)
            if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            for alpha2 in alpha2s:
                CASE_PATH[3] = '/alpha2_'+str(int(alpha2))
                if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                    os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                for n in ns:
                    CASE_PATH[4] = '/n2_'+str(n)
                    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --nx='+str(193)+' --ny='+str(49)+' --nt='+str(2**n)+' --re='+str(re)+' --alpha2='+str(alpha2)+' --fixType='+str(fixType)+' '
                    print case_consts + case_para
                    os.system(exe_pre(npx*npy*npt,' -R lustre ')+EXE_PATH+exe+case_para+case_consts)
                    #os.system(exe_pre(npx*npy*npt)+EXE_PATH+exe+case_para+case_consts)
