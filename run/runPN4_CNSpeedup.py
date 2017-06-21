import os
from platform_paths import *


exe = 'peri_navier4'
runs = 3


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')

CASE_PATH = ['']*5

case_consts = ' --piccard --linSolName="GMRES" --flow=5 --domain=2 --force=0 --tolNOX=1.e-1 --tolBelos=1.e-2   --lx=2. --ly=2.  --fixType=1   '

ns = [6]
res = [10, 100, 200]
alpha2s = [10, 100, 200]
fixTypes = [1, 2, 4, 6, 9, 10]

ns = [4, 5, 6]
ns = [4, 5, 6]
ns = [4]
precTypes = [0]
res = [1]
alpha2s = [1]
fixTypes = [1]
npxs = [1, 1, 1, 2, 2, 2, 4, 4]
npys = [1, 1, 2, 2, 2, 4, 4, 4]
npts = [1, 2, 2, 2, 4, 4, 4, 8]

npxs = [1, 1, 1, 2, 2]
npys = [1, 1, 2, 2, 2]
npts = [1, 2, 2, 2, 4]


for precType in precTypes:
    CASE_PATH[0] = '/speedup'
    if not os.path.exists(DATA_PATH+CASE_PATH[0]):
        os.mkdir(DATA_PATH+CASE_PATH[0])
    for n in ns:
        CASE_PATH[1] = '/ns_'+str(n)
        if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        for re in res:
            #CASE_PATH[2] = '/re_'+str(re)
            CASE_PATH[2] = ''
            if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            for alpha2 in alpha2s:
                #CASE_PATH[3] = '/alpha2_'+str(alpha2)
                CASE_PATH[3] = ''
                if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                    os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                for i in range(len(npxs)):
                    CASE_PATH[4] = '/np_'+str(npxs[i]*npys[i]*npts[i])
                    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --nx='+str(2**n+1)+' --ny='+str(2**n+1)+' --nt='+str(2**n)+' --re='+str(re)+' --alpha2='+str(alpha2)+' --npx='+str(npxs[i])+' --npy='+str(npys[i])+' --npt='+str(npts[i])
                    print exe_pre(npxs[i]*npys[i],' -R "select[model==Opteron8384"] ')+EXE_PATH+exe+case_para+case_consts
                    for run in range(runs):
                        os.system(exe_pre(npxs[i]*npys[i]*npts[i],' -R "select[model==Opteron8380"] ',run=run)+EXE_PATH+exe+case_para+case_consts)
#os.system(exe_pre(npxs[i]*npys[i],run=run)+EXE_PATH+exe+case_para+case_consts)
#os.system(exe_pre(npxs[i]*npys[i],run=run)+EXE_PATH+exe+case_para+case_consts)
