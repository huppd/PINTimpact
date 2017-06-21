import os
from platform_paths import *


exe = 'peri_burgers'


os.chdir(EXE_PATH)
os.system('make -j4')

CASE_PATH = ['']*2

itMs = [1, 2, 4, 6]
case_consts = ' --nf=16 --tolNOX=1.e-6 --tol=1.e-1  --maxI=40 --re=1e3 --linesearch="Polynomial" '

for itM in itMs:
    CASE_PATH[0] = 'itM_'+str(itM)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]):
        os.mkdir(DATA_PATH+CASE_PATH[0])
    print DATA_PATH + CASE_PATH[0]
    os.chdir(DATA_PATH+CASE_PATH[0])
    os.system(' rm -rv ./* ')
    if itM == 1:
        precs = [0, 3, 4, 5, 6, 7, 8, 9]
    elif itM == 2:
        precs = [0, 4, 6, 8, 9]
    elif itM == 4:
        precs = [0, 4, 6, 8, 9]
    elif itM == 6:
        precs = [0, 6, 8, 9]
    for prec in precs:
        CASE_PATH[1] = '/precType_'+str(prec)
        if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        print DATA_PATH + CASE_PATH[0]+CASE_PATH[1]
        os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        os.system(' rm -vr ./* ')
        case_para = ' --iterM='+str(itM)+' --precType='+str(prec)+' '
        print case_consts + case_para
        os.system(exe_pre+EXE_PATH+exe+case_para+case_consts)
