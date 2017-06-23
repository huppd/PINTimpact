import os
from platform_paths import EXE_PATH, DATA_PATH, exe_pre


EXE = 'peri_burgers'


os.chdir(EXE_PATH)
os.system('make -j4')

CASE_PATH = ['']*3

CASE_PATH[0] = 'burgy1D/'

itMs = [1, 2, 4, 6]
case_consts = ' --dim=1 --nx=97 --ny=7 --npx=4 --npy=1 --nfs=1 --nfe=17' + \
    '--tolNOX=1.e-4 --tol=1.e-6  --maxI=20  --linesearch="Polynomial" '

if not os.path.exists(DATA_PATH+CASE_PATH[0]):
    os.mkdir(DATA_PATH+CASE_PATH[0])

for rex in [3, 4, 5, 6]:
    CASE_PATH[1] = 're_1e'+str(rex)+'/'
    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    for itM in itMs:
        CASE_PATH[2] = 'itM_'+str(itM)
        if not os.path.exists(DATA_PATH + CASE_PATH[0] + CASE_PATH[1] +
                              CASE_PATH[2]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        print DATA_PATH + CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]
        os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        os.system(' rm -v ./* ')
        case_para = ' --re=1e'+str(rex)+' --iterM='+str(itM)+' '
        print case_consts + case_para
        os.system(exe_pre+EXE_PATH+EXE+case_para+case_consts)
