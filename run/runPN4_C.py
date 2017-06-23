import os
import platform_paths as pp


EXE = 'peri_navier4'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*5

npx = 1
npy = 2
npt = 8

case_consts = ' --linSolName="GCRODR" --flow=5 --domain=2 --force=0   --npx='+str(npx)+' --npy='+str(npy)+' --npt='+str(npt)+' --tolNOX=1.e-6  --tolBelos=1.e-2  --maxIter=20  --lx=2. --ly=2.  '

precTypes = [0]
ns = [4, 5, 6, 7, 8, 9, 10, 11]
res = [10, 100, 200]
STS = [10, 100, 200]
fixTypes = [1]


#ns  = [7]
#precTypes = [0,1,2]
res = [10]
#res = [200]
STS = [12]
#STS = [15**2]
fixTypes = [1]

for precType in precTypes:
    CASE_PATH[0] = '/precType_2'+str(precType)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0])
    for fixType in fixTypes:
        CASE_PATH[1] = '/fixType_'+str(fixType)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        for re in res:
            CASE_PATH[2] = '/re_'+str(re)
            if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            for st in STS:
                CASE_PATH[3] = '/alpha2_'+str(st)
                if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
                    os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
                for n in ns:
                    CASE_PATH[4] = '/n2_'+str(n)
                    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --nx='+str(2**5+1)+' --ny='+str(2**5+1)+' --nt='+str(2**n)+' --re='+str(re)+' --alpha2='+str(st)+' --fixType='+str(fixType)+' '
                    print case_consts + case_para
                    os.system(pp.exe_pre(npx*npy*npt, ' -R lustre ')+pp.EXE_PATH+EXE+case_para+case_consts)
