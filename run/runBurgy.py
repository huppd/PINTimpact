import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *
from IMPACT_loadfields import pulseC, pulseS


exe = 'peri_burgers'


os.chdir( exe_path )
os.system( 'make -j4' )

case_path = ['','']

case_consts = ' --nf=4 --tolNOX=1.e-1 --tol=1.e-1  --maxI=10 '

for itM in range(1,10):
  case_path[0] = 'itM_'+str(itM)
  if not os.path.exists( data_path+case_path[0] ):
	  os.mkdir( data_path+case_path[0] )
  print data_path + case_path[0]
  os.chdir( data_path+case_path[0] )
  os.system(' rm -v ./* ')
  case_para = ' --iterM='+str(itM)+' '
  print case_consts + case_para
  os.system(exe_pre+exe_path+exe+case_para+case_consts )
