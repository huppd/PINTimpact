""" runner for periodic stokes """
import os
from velocity_profiles import pulseC, pulseS
from numpy import linspace
import numpy as np
from platform_paths import DATA_PATH, EXE_PATH


EXE = 'peri_stokes'


i = 6
N = str(2**i+1)
WOMS = np.array([0.01, 0.05, 0.1, 0.5, 1., 5, 10., 50, 100., 225])
WOMS = 10**np.linspace(-1, 2, 5)
OMS = WOMS
CASE_CONSTS = ' --flow=5 --nx='+N+' --ny='+N+' '


os.chdir(EXE_PATH)
os.system('make -j2')
i = 0
for om in OMS:
    CASE_PATH = 'case'+str(i)
    px = 1./max(max(pulseC(linspace(0, 1., 10000), 1., om, 1.)),
                max(pulseS(linspace(0, 1., 10000), 1., om, 1.)))
    case_para = ' --omega='+str(om)+' --px='+str(px)
    print CASE_CONSTS + case_para
    if not os.path.exists(DATA_PATH+CASE_PATH):
        os.mkdir(DATA_PATH+CASE_PATH)
    os.chdir(DATA_PATH+CASE_PATH)
    os.system('mpirun -np 4 '+EXE_PATH+EXE+case_para+CASE_CONSTS)
    i += 1
