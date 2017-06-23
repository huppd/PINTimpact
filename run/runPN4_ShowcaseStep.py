import os
from pylab import pi
import platform_paths as pp


EXE = 'peri_navier4'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

CASE_PATH = ['']*5

npx = 4
npy = 1
npt = 8

case_consts = ' --linSolName="GCRODR" --flow=1 --domain=1 --force=1   --npx='+str(npx)+' --npy='+str(npy)+' --npt='+str(npt)+' --tolNOX=1.e-6  --tolBelos=1.e-2    --lx=8. --ly=2. --amp=0.2 --radius=0.2 --xm='+str(1./8.)+' '

precTypes = [0]
ns = [5, 6]
res = [10, 100, 200]
STS = [10, 100, 200]
fixTypes = [1]


#ns  = [7]
#precTypes = [0,1,2]
res = [150]
re = 150
steps = range(20, 23, 2)
STS = [12]
STS = [2.*pi*0.1*200, 2.*pi*0.2*200, 2.*pi*0.3*200]
STS = [2.*pi*0.2*150]
fixTypes = [1]
fixType = 1

for precType in precTypes:
    CASE_PATH[0] = '/convshowstep'
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0])
    for st in STS:
        CASE_PATH[1] = '/alpha2_'+str(int(st))
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        for n in ns:
            CASE_PATH[2] = '/n2_'+str(n)
            if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            for step in steps:
                CASE_PATH[3] = '/step_'+str(step)
                if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4]):
                    os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                os.system(' rm ./* -r -v  ')
                case_para = ' --maxIter='+str(step)+' --precType='+str(precType)+' --nx='+str(193)+' --ny='+str(49)+' --nt='+str(2**n)+' --re='+str(re)+' --alpha2='+str(st)+' --fixType='+str(fixType)+' '
                print case_consts + case_para
                os.system(pp.exe_pre(npx*npy*npt,' -R lustre ')+pp.EXE_PATH+EXE+case_para+case_consts)
                #os.system(pp.exe_pre(npx*npy*npt)+pp.EXE_PATH+EXE+case_para+case_consts)
