import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


data_path = '/cluster/home04/math/huppd/trunk/data2/'
exe = 'peri_navier'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']
npx= 1 
npy= 1
npt= 1


case_consts = ' --npx='+str(npx)+' --npy='+str(npy)+' --npf='+str(npt)+' --tolNOX=1.e-6  --tolBelos=1.e-2 --tolInnerBelos=1.e-2 --maxIter=20  --lx=2. --ly=2.  --initZero=0 --numCycles=2 '#--domain=1 --linSolver=GCRODR '

precTypes = [ 0, 1, 2 ]
precTypes = [ 2 ]
ns        = [ 4, 5, 6, 7 ]
#ns        = [  5, 6 ]
res       = [ 1, 25, 50, 75, 100 ]
res = np.arange(1,220,20)
alpha2s   = [ 1, 9, 25, 64, 100 ]
alpha2s   = np.arange(1,16,2)**2

#precTypes = [ 2, 3 ]
#precTypes = [ 2 ]
#precTypes = [ 3 ]
#ns = [ 5 ]
#res       = [ 100 ]
#alpha2s   = [ 100 ]

#data_path = '/cluster/scratch_xp/public/huppd'

#ns  = [ 6, 7 ]

for alpha2 in alpha2s:
	case_path[0] = '/a2_'+str(alpha2)
	if not os.path.exists( data_path+case_path[0] ):
		os.mkdir( data_path+case_path[0] )
	for re in res:
		case_path[1] = '/re_'+str(re)
		if not os.path.exists( data_path+case_path[0]+case_path[1] ):
			os.mkdir( data_path+case_path[0]+case_path[1] )
		for n in ns:
			case_path[2] = '/n2_'+str(n)
			if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
				os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
			for precType in precTypes:
				case_path[3] = '/precType_'+str(precType)
				if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]   ):
					os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
				os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
				os.system(' rm ./* -r -v  ')
				case_para = ' --nx='+str(2**n+1)+' --ny='+str(2**n+1)+' --nf='+str(2*(n-1))+' --withprec='+str(precType)+'  --re='+str(re)+' --alpha2='+str(alpha2)+' --maxGrids='+str(n-2)+' ' 
				#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
				print( exe_pre(npx*npy*npt)+exe_path+exe+case_para+case_consts +' > output ' )
				os.system( exe_pre(npx*npy*npt)+exe_path+exe+case_para+case_consts +' > output ' )
