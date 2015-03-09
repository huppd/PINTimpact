import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


exe = 'peri_burgers'


os.chdir( exe_path )
os.system( 'make -j4' )

case_path = ['varNX/','','']

itMs = [1,2,4,6]
case_consts = ' --dim=1 --ny=7 --npx=4 --npy=1 --nfs=1 --nfe=17 --tolNOX=1.e-1 --tol=1.e-6  --iterM=1 --maxI=20  --linesearch="Polynomial" '

if not os.path.exists( data_path+case_path[0] ):
  os.mkdir( data_path+case_path[0] )
  
for nx in range(4,8):
	case_path[1] = 'nx_22'+str(nx)+'/'
	if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		os.mkdir( data_path+case_path[0]+case_path[1] )
	for itM in itMs:
		case_path[2] = 'itM_'+str(itM)
		if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
			os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
		print data_path + case_path[0]+case_path[1]+case_path[2]
		os.chdir( data_path+case_path[0]+case_path[1]+case_path[2] )
		os.system(' rm -v ./* ')
		case_para = ' --re=1e'+str(rex)+'  --nx='+str(2**nx+1)+' '
		print case_consts + case_para
		os.system(exe_pre+exe_path+exe+case_para+case_consts )
