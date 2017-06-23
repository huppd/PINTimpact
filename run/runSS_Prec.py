import os
import numpy as np
import platform_paths as pp


EXE = 'stat_stokes'
EXE = 'peri_stokes'

sol = 'GMRES'
#sol='GCRODR'

os.chdir(pp.EXE_PATH)
os.system('make -j4')

woms = np.array([0.01, 0.05, 0.1, 0.5, 1., 5, 10., 50, 100., 225])
woms = 10**np.linspace(-1, 2, 5)

oms = woms

CASE_PATH = ['']*3

for i in range(5, 7):
    n = str(2**i+1)
    case_consts = ' --flow=3 --nx='+n+' --ny='+n+' '
    #
    CASE_PATH[0] = 'discr'+str(i)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0])
    print pp.DATA_PATH + CASE_PATH[0]
    for prec in [0, 2, 3, 4, 5, 6]:
        CASE_PATH[1] = '/prec'+str(prec)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        print pp.DATA_PATH + CASE_PATH[0]+CASE_PATH[1]
        for j, om in enumerate(oms):
            CASE_PATH[2] = '/case'+str(j)
            print pp.DATA_PATH + CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]
            if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
                    os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            os.system(' rm -v ./* ')
            case_para = ' --px='+str(1)+' --omega='+str(om)+' --solver1='+sol+' --prec='+str(prec)+' '
            print case_consts + case_para
            os.system(pp.exe_pre()+pp.EXE_PATH+EXE+case_para+case_consts)
