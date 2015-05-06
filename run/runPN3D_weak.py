import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


#data_path = '/cluster/home04/math/huppd/trunk/data/'
#data_path = '/cluster/scratch_xp/public/huppd'
#data_path = '/cluster/scratch_xp/public/huppd'
exe = 'peri_navier3D'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']
npx= 1 
npy= 1
npt= 1


case_consts = [' --baseflow=1 --flow=17 --alpha2=9  --re=1500  --domain=1 --nf=4  --tolNOX=1.e-6 ',
' --nf=4   --baseflow=22 --flow=0 --re=300 --alpha2=9 --tolNOX=1.e-6  --tolBelos=1.e-1 --tolInnerBelos=1.e-3 --maxIter=20  --lx=30. --ly=30. --ly=30 --initZero=0 --numCycles=2 --domain=0 ']

ns        = [ 1 ]
ns        = [ 1, 2 ]
#res       = [ 1, 25, 50, 75, 100 ]
#alpha2s   = [ 1, 9, 25, 64, 100 ]

#precTypes = [ 2, 3 ]
#precTypes = [ 2 ]
#precTypes = [ 3 ]
#ns = [ 5 ]
#res       = [ 1600 ]
#res       = [ 400 ]
#alpha2s   = [ 125 ]


#ns  = [ 6, 7 ]

#for alpha2 in alpha2s:
	#case_path[0] = '/a2_'+str(alpha2)
	#if not os.path.exists( data_path+case_path[0] ):
		#os.mkdir( data_path+case_path[0] )
for case in range(2):
	case_path[1] = '/case_'+str(case)
	if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		os.mkdir( data_path+case_path[0]+case_path[1] )
	for n in ns:
		case_path[2] = '/n2_'+str(n)
		if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
			os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
		#for precType in precTypes:
			#case_path[3] = '/precType_'+str(precType)
			#if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]   ):
				#os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.system(' rm ./* -r -v  ')
		case_para = ' --nx='+str(32*n+1)+' --ny='+str(16*n+1)+' --nz='+str(32*n+1)+' --npx='+str(n)+' --npy='+str( max(n/2,1) )+' --npz='+str(n)+' --maxGrids='+str(4+n)+' '
		#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
		print( exe_pre(n*n*max(n/2,1))+exe_path+exe+case_para+case_consts[case] +' > output ' )
		os.system( exe_pre(n*n*max(n/2,1))+exe_path+exe+case_para+case_consts[case] +'  ' )
