import os
from pylab import pi
import platform_paths as pp


EXE = 'peri_navier4'
RUNS = 4


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

CASE_PATH = ['']*5

case_consts = ' --linSolName="GCRODR" --piccard --flow=1 --domain=1 --force=1 --radius=0.2 --amp=0.2  --tolNOX=1.e-2 --tolBelos=1.e-1  --maxIter=20  --ly=2. --xm='+str(0.25) + '  '


ns = [3, 4, 5]

precTypes = [0]
res = [150]
STS = [2.*pi*0.2*res[0]]
fixTypes = [1]
npxs = [1, 1, 1, 2, 2, 4, 4]
npys = [1, 1, 2, 2, 2, 4, 4]
npts = [1, 2, 2, 2, 4, 4, 8]

npxs = [1, 1, 1, 2, 2, 4]
npys = [1, 1, 2, 2, 2, 4]
npts = [1, 2, 2, 2, 4, 4]

npxs = [1, 1, 2, 2, 4]
npys = [1, 1, 1, 1, 1]
npts = [1, 2, 2, 4, 4]

#npxs = [4]
#npys = [4]
#npts = [4]


for precType in precTypes:
    CASE_PATH[0] = '/weakscaling'
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0])
    for n in ns:
        CASE_PATH[1] = '/ns_'+str(n)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        for re in res:
            #CASE_PATH[2] = '/re_'+str(re)
            CASE_PATH[2] = ''
            if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            for st in STS:
                #CASE_PATH[3] = '/alpha2_'+str(st)
                CASE_PATH[3] = ''
                if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                    os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                for i in range(len(npxs)):
                    CASE_PATH[4] = '/np_'+str(npxs[i]*npys[i]*npts[i])
                    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --lx='+str(4*npxs[i])+' --nx='+str(npxs[i]*2*2**n+1)+' --ny='+str(2**n+1)+' --nt='+str(npts[i]*2**n)+' --re='+str(re)+' --alpha2='+str(st)+' --npx='+str(npxs[i])+' --npy='+str(npys[i])+' --npt='+str(npts[i])
                    print pp.exe_pre(npxs[i]*npys[i],' -R "select[model==Opteron8384"] ')+pp.EXE_PATH+EXE+case_para+case_consts
                    for run in range(RUNS, 10):
                        os.system(pp.exe_pre(npxs[i]*npys[i]*npts[i], ' -R "select[model==Opteron8380"] ', run=run)+pp.EXE_PATH+EXE+case_para+case_consts)
