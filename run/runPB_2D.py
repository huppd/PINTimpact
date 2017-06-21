import os
from platform_paths import *


exe = 'peri_burgers'


os.chdir(EXE_PATH)
os.system('make -j4')

CASE_PATH = ['']*2

itMs = [1, 2, 4, 6]
case_consts = ' --nfe=17 --tolNOX=1.e-1 --tol=1.e-6  --maxI=40 --re=1e5 --linesearch="Polynomial" '

for itM in itMs:
    CASE_PATH[0] = 'itM_'+str(itM)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]):
        os.mkdir(DATA_PATH+CASE_PATH[0])
    print DATA_PATH + CASE_PATH[0]
    os.chdir(DATA_PATH+CASE_PATH[0])
    os.system(' rm -v ./* ')
    case_para = ' --iterM='+str(itM)+' '
    print case_consts + case_para
    os.system(exe_pre+EXE_PATH+exe+case_para+case_consts)
