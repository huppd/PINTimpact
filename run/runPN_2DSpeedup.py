import os
import platform_paths as pp


EXE = 'peri_navier'
runs = 4


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

CASE_PATH = ['']*5

case_consts = ' --linSolName="GMRES" --flow=1 --domain=1 --force=4 --radius=0.1 --rotation=2 --nf=2 --tolNOX=1.e-1 --tolNF=1.e-1 --tolBelos=1.e-6  --maxIter=2  --lx=8. --ly=2.  --fixType=1 --xm='+str(1./8.) + '  '

precTypes = [0, 10]
ns = [6]
res = [10, 100, 200]
STS = [10, 100, 200]
fixTypes = [1, 2, 4, 6, 9, 10]

ns = [4, 5, 6]
ns = [4, 5, 6, 7]
precTypes = [0]
res = [100]
STS = [10**2]
fixTypes = [1]
npxs = [1, 2, 4, 4, 8, 16, 16, 32]
npys = [1, 1, 1, 2, 2,  2,  4,  4]

npxs = [1, 2, 4, 4, 8]
npys = [1, 1, 1, 2, 2]
#npxs = [16]
#npys = [4]

for precType in precTypes:
    CASE_PATH[0] = '/speedup3'
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
                    CASE_PATH[4] = '/np_'+str(npxs[i]*npys[i])
                    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4] ):
                        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]+CASE_PATH[4])
                    os.system(' rm ./* -r -v  ')
                    case_para = ' --precType='+str(precType)+' --nx='+str(4*2**n+1)+' --ny='+str(2**n+1)+' --re='+str(re)+' --alpha2='+str(st)+' --npx='+str(npxs[i])+' --npy='+str(npys[i])
                    print pp.exe_pre(npxs[i]*npys[i], ' -R "select[model==Opteron8384"] ')+pp.EXE_PATH+EXE+case_para+case_consts
                    for run in range(runs):
                        os.system(pp.exe_pre(npxs[i]*npys[i], ' -R "select[model==Opteron8380"] ', run=run)+pp.EXE_PATH+EXE+case_para+case_consts)
