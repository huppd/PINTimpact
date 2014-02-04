import os

exe_path = '/home/huppd/workspace/PIMPACT/RELEASE/src/src_c/'
exe = 'stat_stokes'

data_path = '/home/huppd/workspace/PIMPACT/data/'

case_paths = ['case1','case2','case3']

case_paras = []
for i in [5,6,7]:
	n = 2**i+1
	case_paras.append( ' --nx='+str(n)+' --ny='+str(n) )

for i in range(3):
	os.chdir( data_path+case_paths[i] )
	os.system('mpirun -np 4 '+exe_path+exe+case_paras[i] )
