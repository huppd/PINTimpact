import os
from numpy import linspace
import numpy as np
from platform_paths import *
from IMPACT_loadfields import pulseC, pulseS


exe = 'stat_stokes'
exe = 'peri_stokes'

sol='GMRES'
#sol='GCRODR'

os.chdir(exe_path)
os.system('make -j4')

woms = np.array([0.01,0.05,0.1,0.5,1.,5,10.,50,100.,225])
woms = 10**np.linspace(-1,3,20)

oms  = woms

case_path = ['','','']

for i in range(4,8):
  n = str(2**i+1)
  case_consts = ' --flow=3 --nx='+n+' --ny='+n+' '
  
  case_path[0] = 'discr'+str(i)
  if not os.path.exists( data_path+case_path[0] ):
	  os.mkdir( data_path+case_path[0] )
  print data_path + case_path[0]
  for prec in [0,2,3]:
  #for prec in [0]:
	  case_path[1] = '/prec'+str(prec)
	  if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		  os.mkdir( data_path+case_path[0]+case_path[1] )
	  print data_path + case_path[0]+case_path[1]
	  j = 0
	  for om in oms:
		  case_path[2] = '/case'+str(j)
		  print data_path + case_path[0]+case_path[1]+case_path[2]
		  if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
			  os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
		  os.chdir( data_path+case_path[0]+case_path[1]+case_path[2])
		  os.system(' rm -v ./* ')
		  #px = 1./max(max( pulseC(linspace(0,1.,10000),1.,om,1.) ),max( pulseS(linspace(0,1.,10000),1.,om,1.) ) )
		  case_para = ' --px='+str(1)+' --omega='+str(om)+ ' --solver1='+sol+' --prec='+str(prec)+' '
		  print case_consts + case_para
		  os.system(exe_pre+exe_path+exe+case_para+case_consts )
		  j += 1