import os
from numpy import linspace
import numpy as np
from platform_paths import *
from IMPACT_loadfields import pulseC, pulseS


exe = 'stat_stokes'
exe = 'peri_stokes'


i = 6
n = str(2**i+1)

woms = np.array([0.01,0.05,0.1,0.5,1.,5,10.,50,100.,225])
woms = 10**np.linspace(-2,2,5)

oms  = woms

case_consts = ' --flow=3 --nx='+n+' --ny='+n+' '

os.chdir(exe_path)
os.system('make -j2')

sol='GMRES'
#sol='GCRODR'

for prec in range(0,4):
	case_path0 = 'prec'+str(prec)
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
		os.system(' rm -v ./* ')
		px = 1./max(max( pulseC(linspace(0,1.,10000),1.,om,1.) ),max( pulseS(linspace(0,1.,10000),1.,om,1.) ) )
		case_para = ' --px='+str(px)+' --omega='+str(om)+ ' --solver1='+sol+' --prec='+str(prec)+' '
		print case_consts + case_para
		os.system(exe_pre+exe_path+exe+case_para+case_consts )
		i += 1
