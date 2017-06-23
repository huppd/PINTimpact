import os
import numpy as np
import platform_paths as pp


EXE = 'stat_stokes'
EXE = 'peri_stokes'

i = 6
n = str(2**i+1)
woms = np.array([0.01, 0.05, 0.1, 0.5, 1., 5, 10., 50, 100., 225])
woms = 10**np.linspace(-2, 3, 5)
oms = woms

case_consts = ' --domain=2 --flow=5 --nx='+n+' --ny='+n+' '

os.chdir(pp.EXE_PATH)
os.system('make -j2')


for sol in ["GMRES", "GCRODR"]:
    CASE_PATH0 = sol
    if not os.path.exists(pp.DATA_PATH+CASE_PATH0):
        os.mkdir(pp.DATA_PATH+CASE_PATH0)
    print pp.DATA_PATH + CASE_PATH0
    i = 0
    for om in oms:
        CASE_PATH1 = '/case'+str(i)
        print pp.DATA_PATH + CASE_PATH0 + CASE_PATH1
        if not os.path.exists(pp.DATA_PATH+CASE_PATH0+CASE_PATH1):
            os.mkdir(pp.DATA_PATH+CASE_PATH0+CASE_PATH1)
        os.chdir(pp.DATA_PATH+CASE_PATH0+CASE_PATH1)
        case_para = ' --omega='+str(om) + ' --solver1='+sol+' '
        print case_consts + case_para
        os.system(pp.exe_pre+pp.EXE_PATH+EXE+case_para+case_consts)
        i += 1
