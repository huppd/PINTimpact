import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


exe = 'peri_navier'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','']

case_consts = ' --flow=1 --domain=1 --force=3 --rad=0.2 --rot=10 --nfe=17  --nx=33 --ny=64 --tolNOX=1.e-6 --tolBelos=1.e-4 --tolNF=1.e-6 --maxIter=10  --domain=2 --lx=2. --ly=4.  '

flows   = [ 1 ]
res     = [ 1, 10, 100, 500 ]
alpha2s = [ 1, 10, 100, 200, 800 ]

for flow in flows:
  case_path[0] = 'disc'
  if not os.path.exists( data_path+case_path[0] ):
    os.mkdir( data_path+case_path[0] )
  print data_path + case_path[0]
  for alpha2 in alpha2s:
    case_path[1] = '/alpha2_'+str(alpha2)
    if not os.path.exists( data_path+case_path[0]+case_path[1] ):
      os.mkdir( data_path+case_path[0]+case_path[1] )
    print data_path + case_path[0] + case_path[1]
    os.chdir( data_path+case_path[0]+case_path[1] )
    os.system(' rm -rv ./* ')
    for re in res:
      case_path[2] = '/re_'+str(re)
      if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
	os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
      print data_path + case_path[0]+case_path[1]+case_path[2]
      os.chdir( data_path+case_path[0]+case_path[1]+case_path[2] )
      os.system(' rm -vr ./* ')
      case_para = ' --alpha2='+str(alpha2)+' --re='+str(re)+' '#+'  --flow='+str(flow)+' '
      print case_consts + case_para
      os.system( exe_pre+exe_path+exe+case_para+case_consts )
