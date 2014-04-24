import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


exe = 'peri_navier'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','']

case_consts = ' --nfe=17   --maxIter=10 --flow=7  --domain=2 --lx=2. --ly=2.  '

res     = [10, 100, 1000]
alpha2s = [10, 100, 1000,10000]

for alpha2 in alpha2s:
  case_path[0] = 'alpha2_'+str(alpha2)
  if not os.path.exists( data_path+case_path[0] ):
	  os.mkdir( data_path+case_path[0] )
  print data_path + case_path[0]
  os.chdir( data_path+case_path[0] )
  os.system(' rm -rv ./* ')
  for re in res:
    case_path[1] = '/re_'+str(re)
    if not os.path.exists( data_path+case_path[0]+case_path[1] ):
	    os.mkdir( data_path+case_path[0]+case_path[1] )
    print data_path + case_path[0]+case_path[1]
    os.chdir( data_path+case_path[0]+case_path[1] )
    os.system(' rm -vr ./* ')
    case_para = ' --alpha2='+str(alpha2)+' --re='+str(re)+' '
    print case_consts + case_para
    os.system(exe_pre+exe_path+exe+case_para+case_consts )