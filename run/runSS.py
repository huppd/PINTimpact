import os
import platform_paths as pp

exe = 'stat_stokes'

CASE_PATHs = ['case1', 'case2', 'case3']

case_paras = []
for i in [5, 6, 7]:
    n = 2**i+1
    case_paras.append(' --nx='+str(n)+' --ny='+str(n))

for i in range(3):
    os.chdir(pp.DATA_PATH+CASE_PATHs[i])
    os.system('mpirun -np 4 '+pp.EXE_PATH+exe+case_paras[i])
