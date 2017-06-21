import os
from numpy import linspace
import numpy as np
from platform_paths import *


#EXE_PATH = '/home/huppd/workspace/pimpact-repo/release/src/src_c/'
exe = 'stat_stokes'
exe = 'peri_stokes'

i = 6
n = str(2**i+1)
woms = np.array([0.01, 0.05, 0.1, 0.5, 1., 5, 10., 50, 100., 225])
woms = 10**np.linspace(-2, 3, 5)
oms = woms

case_consts = ' --domain=2 --flow=5 --nx='+n+' --ny='+n+' '

os.chdir(EXE_PATH)
os.system('make -j2')


for sol in ["GMRES","GCRODR"]:
    CASE_PATH0 = sol
    if not os.path.exists(DATA_PATH+CASE_PATH0):
        os.mkdir(DATA_PATH+CASE_PATH0)
    print DATA_PATH + CASE_PATH0
    i = 0
    for om in oms:
        CASE_PATH1 = '/case'+str(i)
        print DATA_PATH + CASE_PATH0 + CASE_PATH1
        if not os.path.exists(DATA_PATH+CASE_PATH0+CASE_PATH1):
            os.mkdir(DATA_PATH+CASE_PATH0+CASE_PATH1)
        os.chdir(DATA_PATH+CASE_PATH0+CASE_PATH1)
        case_para = ' --omega='+str(om)+ ' --solver1='+sol+' '
        print case_consts + case_para
        os.system(exe_pre+EXE_PATH+exe+case_para+case_consts)
        i += 1
