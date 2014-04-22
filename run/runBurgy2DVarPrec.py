import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


exe = 'peri_burgers'


os.chdir( exe_path )
os.system( 'make -j4' )

case_path = ['','']

itMs = [1,2,4,6]
case_consts = ' --nf=16 --tolNOX=1.e-6 --tol=1.e-1  --maxI=40 --re=1e3 --linesearch="Polynomial" '

for itM in itMs:
  case_path[0] = 'itM_'+str(itM)
  if not os.path.exists( data_path+case_path[0] ):
	  os.mkdir( data_path+case_path[0] )
  print data_path + case_path[0]
  os.chdir( data_path+case_path[0] )
  os.system(' rm -rv ./* ')
  if itM==1:
    precs = [0,3,4,5,6,7,8,9]
  elif itM==2:
    precs = [0,4,6,8,9]
  elif itM==4:
    precs = [0,4,6,8,9]
  elif itM==6:
    precs = [0,6,8,9]
  for prec in precs:
    case_path[1] = '/precType_'+str(prec)
    if not os.path.exists( data_path+case_path[0]+case_path[1] ):
	    os.mkdir( data_path+case_path[0]+case_path[1] )
    print data_path + case_path[0]+case_path[1]
    os.chdir( data_path+case_path[0]+case_path[1] )
    os.system(' rm -vr ./* ')
    case_para = ' --iterM='+str(itM)+' --precType='+str(prec)+' '
    print case_consts + case_para
    os.system(exe_pre+exe_path+exe+case_para+case_consts )