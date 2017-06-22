import os
from IMPACT_loadfields import pulseC, pulseS
from numpy import linspace
import numpy as np
from platform_paths import DATA_PATH, EXE_PATH


exe = 'stat_stokes'
exe = 'peri_stokes'


i = 6
n = str(2**i+1)
woms = np.array([0.01, 0.05, 0.1, 0.5, 1., 5, 10., 50, 100., 225])
woms = 10**np.linspace(-1, 2, 5)
oms = woms
case_consts = ' --flow=5 --nx='+n+' --ny='+n+' '


os.chdir(EXE_PATH)
os.system('make -j2')
i = 0
for om in oms:
    CASE_PATH = 'case'+str(i)
    px = 1./max(max(pulseC(linspace(0, 1., 10000), 1., om, 1.)), max(pulseS(linspace(0, 1., 10000), 1., om, 1.)))
    case_para = ' --omega='+str(om)+' --px='+str(px)
    print case_consts + case_para
    if not os.path.exists(DATA_PATH+CASE_PATH):
        os.mkdir(DATA_PATH+CASE_PATH)
    os.chdir(DATA_PATH+CASE_PATH)
    os.system('mpirun -np 4 '+EXE_PATH+exe+case_para+case_consts)
    i += 1
