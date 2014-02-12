import os
from numpy import linspace
import numpy as np
from platform_paths import *


#exe_path = '/home/huppd/workspace/pimpact-repo/release/src/src_c/'
exe = 'stat_stokes'
exe = 'peri_stokes'

#data_path = '/home/huppd/workspace/pimpact-repo/data/'

i = 6
n = str(2**i+1)
woms = np.array([0.01,0.05,0.1,0.5,1.,5,10.,50,100.,225])
woms = 10**np.linspace(-2,3,20)
oms  = woms

case_consts = ' --domain=2 --flow=5 --nx='+n+' --ny='+n+' '

os.chdir(exe_path)
os.system('make -j2')


for sol in ["GMRES","GCRODR"]:
	case_path0 = sol
	if not os.path.exists( data_path+case_path0 ):
		os.mkdir( data_path+case_path0 )
	print data_path + case_path0
	i = 0
	for om in oms:
		case_path1 = '/case'+str(i)
		print data_path + case_path0 + case_path1
		if not os.path.exists( data_path+case_path0+case_path1 ):
			os.mkdir( data_path+case_path0+case_path1 )
		os.chdir( data_path+case_path0+case_path1 )
		case_para = ' --omega='+str(om)+ ' --solver1='+sol+' '
		print case_consts + case_para
		os.system(exe_pre+exe_path+exe+case_para+case_consts )
		i += 1
