import os
import platform_paths as pp


EXE = 'peri_burgers'


os.chdir(pp.EXE_PATH)
os.system('make -j4')

CASE_PATH = ['']*3

CASE_PATH[0] = 'varNX/'

itMs = [1, 2, 4, 6]
case_consts = ' --dim=1 --ny=7 --npx=4 --npy=1 --nfs=1 --nfe=17 --tolNOX=1.e-1 --tol=1.e-6  --iterM=1 --maxI=20  --linesearch="Polynomial" '

if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])

for nx in range(4, 8):
    CASE_PATH[1] = 'nx_22'+str(nx)+'/'
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    for itM in itMs:
        CASE_PATH[2] = 'itM_'+str(itM)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        print pp.DATA_PATH + CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]
        os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
        os.system(' rm -v ./* ')
        case_para = ' --re=1e'+str(rex)+'  --nx='+str(2**nx+1)+' '
        print case_consts + case_para
        os.system(pp.exe_pre+pp.EXE_PATH+EXE+case_para+case_consts)
