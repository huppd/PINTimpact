import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *
from IMPACT_loadfields import pulseC, pulseS


exe = 'peri_burgers'

sol='GMRES'
#sol='GCRODR'

os.chdir(exe_path)
os.system('make -j4')

woms = 4.*pi

oms  = woms

case_path = ['','','']

for i in range(7,8):
  n = str(2**i+1)
  case_consts = ' --nx='+n+' --ny=7 '
  
  case_path[0] = 'discr'+str(i)
  if not os.path.exists( data_path+case_path[0] ):
	  os.mkdir( data_path+case_path[0] )
  print data_path + case_path[0]
  for nf in range(16,17,1):
    case_path[1] = '/nf'+str(nf)
    if not os.path.exists( data_path+case_path[0]+case_path[1] ):
      os.mkdir( data_path+case_path[0]+case_path[1] )
    print data_path + case_path[0]+case_path[1]
    print data_path + case_path[0]+case_path[1]+case_path[2]
    if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
      os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
    os.chdir( data_path+case_path[0]+case_path[1]+case_path[2])
    os.system(' rm -v ./* ')
    case_para = ' --nf='+str(nf)+' --alpha2='+str(4.*pi)+ '  --tol=0.1 --re=1000 --npx=4 --npy=1 --domain=4 '
    print case_consts + case_para
    os.system(exe_pre+exe_path+exe+case_para+case_consts )
