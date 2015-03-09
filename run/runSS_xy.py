import subprocess
import signal
import os
from IMPACT_loadfields import pulseC, pulseS
from numpy import linspace


exe_path = '/home/huppd/workspace/pimpact-repo/release/src/src_c/'
exe = 'stat_stokes'
exe = 'peri_stokes'

data_path = '/home/huppd/workspace/pimpact-repo/data/'

case_paths = ['casex','casey']

case_paras = [' --flow=1',' --flow=2']
case_paras = [' --flow=3',' --flow=4']

i = 4
n = str(2**i+1)
om = 100.
px = 1./max(max( pulseC(linspace(0,1.,10000),1.,om,1.) ),max( pulseS(linspace(0,1.,10000),1.,om,1.) ) )
case_consts = ' --omega='+str(om)+' --px='+str(px)+ ' --nx='+n+' --ny='+n+' '

os.chdir(exe_path)
subprocess.call('make -j2',shell=True)
print( case_consts )
for i in range(2):
	if not os.path.exists( data_path+case_paths[i] ):
		os.mkdir( data_path+case_paths[i] )
	os.chdir( data_path+case_paths[i] )
	#os.system('/usr/bin/mpirun -np 4 '+exe_path+exe+case_paras[i]+case_consts ,shell=True)
	#subprocess.call('/usr/bin/mpirun -np 4 '+exe_path+exe+case_paras[i]+case_consts ,shell=True )
	mpicmd =  'mpirun -np 4 '+exe_path+exe+case_paras[i]+case_consts 
	#subprocess.call(mpicmd, shell=True)
	#subprocess.Popen(mpicmd, shell=True)
	subprocess.check_call(mpicmd, shell=True)
	#mpi = subprocess.Popen(mpicmd, 
                       #shell=True, 
                       #stdout=subprocess.PIPE, 
                       #stderr=subprocess.PIPE)
	#mpi.wait()
	#stdout, stderr = mpi.communicate()
	#print( stdout )
