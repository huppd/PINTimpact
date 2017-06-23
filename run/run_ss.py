""" runner for static stokes """
import os
import platform_paths as pp

EXE = 'stat_stokes'

CASE_PATHS = ['case1', 'case2', 'case3']

CASE_PARAS = []
for i in [5, 6, 7]:
    n = 2**i+1
    CASE_PARAS.append(' --nx='+str(n)+' --ny='+str(n))

for i in range(3):
    os.chdir(pp.DATA_PATH+CASE_PATHS[i])
    os.system('mpirun -np 4 '+pp.EXE_PATH+EXE+CASE_PARAS[i])
