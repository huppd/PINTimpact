import os
import platform_paths as pp


EXE = 'peri_navier3D'

RUNS = range(0, 4)

os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*5

npx = 1
npy = 1
npt = 1


case_consts = ['--withoutput=0 --baseflow=1 --force=1  --domain=0 --nf=4  ' +
               '--tolNOX=1.e-6  --tolInnerBelos=1.e-1 --maxIter=2 --lx=8. ' +
               '--ly=2. --lz=4.']

nps = [4, 8, 12]
nxs = [1, 2, 3]
#ns = [ 2, 3, 4, 6, 8,12]
#ns = [1, 2]
#nps = [8]
#nxs = [2]
#nps = [4, 8]
re = 100
st = 576


CASE_PATH[0] = '/strong2'
if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])
#os.chdir(pp.DATA_PATH+CASE_PATH[0])
#os.system(' rm ./* -r -v  ')
#for st in sts:
    #CASE_PATH[1] = '/a2_'+str(st)
    #if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        #os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
for n in nxs:
    CASE_PATH[2] = '/n_'+str(n)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
    for np in nps:
        CASE_PATH[3] = '/np_'+str(np)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
        os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]+CASE_PATH[3])
        os.system(' rm ./* -r -v  ')
        case_para = ' --nx='+str(128*n+1)+' --ny='+str(32*n+1)+' --nz='+str(64*n+1)+' --npx='+str(np)+' --npy='+str(max(np/4,1))+' --npz='+str(max(np/2,1))+' --maxGrids='+str(int(2+n))+' --re='+str(re)+' --alpha2='+str(st)+' '
        for run in RUNS:
            exe_str = pp.exe_pre(np*max(np/2, 1)*max(np/4, 1),
                                 ' -R "select[model==Opteron6174"] ' +
                                 ' -R "rusage[mem=8192]" -W 2:00', run) + \
                pp.EXE_PATH + EXE + case_para + case_consts[0] + '  '
            print(exe_str)
            os.system(exe_str)
