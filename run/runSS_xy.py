import subprocess
# import signal
import os
from IMPACT_loadfields import pulseC, pulseS
from numpy import linspace
import platform_paths as pp


EXE = 'stat_stokes'
EXE = 'peri_stokes'


CASE_PATHs = ['casex', 'casey']

case_paras = [' --flow=1', ' --flow=2']
case_paras = [' --flow=3', ' --flow=4']

i = 4
n = str(2**i+1)
om = 100.
px = 1./max(max(pulseC(linspace(0, 1., 10000), 1., om, 1.)),
            max(pulseS(linspace(0, 1., 10000), 1., om, 1.)))
case_consts = ' --omega='+str(om)+' --px='+str(px) + ' --nx='+n+' --ny='+n+' '

os.chdir(pp.EXE_PATH)
subprocess.call('make -j2', shell=True)
print(case_consts)
for i in range(2):
    if not os.path.exists(pp.DATA_PATH+CASE_PATHs[i]):
        os.mkdir(pp.DATA_PATH+CASE_PATHs[i])
    os.chdir(pp.DATA_PATH+CASE_PATHs[i])
    #os.system('/usr/bin/mpirun -np 4 '+pp.EXE_PATH+EXE+case_paras[i]+case_consts ,shell=True)
    #subprocess.call('/usr/bin/mpirun -np 4 '+pp.EXE_PATH+EXE+case_paras[i]+case_consts ,shell=True)
    mpicmd = 'mpirun -np 4 '+pp.EXE_PATH+EXE+case_paras[i]+case_consts
    #subprocess.call(mpicmd, shell=True)
    #subprocess.Popen(mpicmd, shell=True)
    subprocess.check_call(mpicmd, shell=True)
    #mpi = subprocess.Popen(mpicmd,
                   #shell=True,
                   #stdout=subprocess.PIPE,
                   #stderr=subprocess.PIPE)
    #mpi.wait()
    #stdout, stderr = mpi.communicate()
    #print(stdout)
