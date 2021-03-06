import os
import platform_paths as pp


EXE = 'peri_navier4'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

CASE_PATH = ['']*5

npx = 8
npy = 1
npt = 4

case_consts = ' --linSolName="GMRES" --piccard --flow=1 --domain=1 --force=1 --radius=0.1 --amp=0.1  --npx='+str(npx)+' --npy='+str(npy)+' --npt='+str(npt)+' --tolNOX=1.e-2 --tolBelos=1.e-1  --maxIter=20  --lx=4. --ly=2. --xm='+str(0.25) + '  '

precTypes = [0, 10]
ns = [4]
res = [10, 100, 200]
STS = [10, 100, 200]
fixTypes = [1, 2, 4, 6, 9, 10]


ns = [5]
precTypes = [0, 1, 2]
res = [200]
STS = [251]
fixTypes = [1]

for precType in precTypes:
    CASE_PATH[0] = 'precType_'+str(precType)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0])
    for n in ns:
        CASE_PATH[1] = '/n2_'+str(n)
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
                for fixType in fixTypes:
                    CASE_PATH[4] = '/fixType_'+str(fixType)
                    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --nx='+str(2*2**n+1)+' --ny='+str(2**n+1)+' --nt='+str(2**(n-1)+1)+' --re='+str(re)+' --alpha2='+str(st)+' --fixType='+str(fixType)+' '
                    print case_consts + case_para
                    os.system(pp.exe_pre(npx*npy*npt, ' -R lustre ')+pp.EXE_PATH+EXE+case_para+case_consts)
