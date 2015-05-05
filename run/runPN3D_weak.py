import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


data_path = '/cluster/home04/math/huppd/trunk/data/'
exe = 'peri_navier3D'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']
npx= 1 
npy= 1
npt= 1


case_consts = '  --maxGrids=10 --baseflow=1 --flow=17 --re=1600 --alpha2=125 --tolNOX=1.e-6  --tolBelos=1.e-1 --tolInnerBelos=1.e-3 --maxIter=20  --lx=4. --ly=2. --ly=4 --initZero=0 --numCycles=2 --domain=1 '
case_consts = ' --maxGrids=10  --baseflow=22 --flow=23 --re=400 --alpha2=125 --tolNOX=1.e-6  --tolBelos=1.e-1 --tolInnerBelos=1.e-3 --maxIter=20  --lx=30. --ly=30. --ly=30 --initZero=0 --numCycles=2 --domain=0 --dim=3 '

precTypes = [ 2 ]
ns        = [ 6 ]
ns        = [ 1, 2, 3, 4 ]
res       = [ 1, 25, 50, 75, 100 ]
alpha2s   = [ 1, 9, 25, 64, 100 ]

#precTypes = [ 2, 3 ]
#precTypes = [ 2 ]
#precTypes = [ 3 ]
#ns = [ 5 ]
res       = [ 1600 ]
res       = [ 400 ]
alpha2s   = [ 125 ]

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
				case_para = ' --nx='+str(64*n+1)+' --ny='+str(32*n+1)+' --nz='+str(64*n+1)+' --nf=8 --withprec=2  --npx='+str(n)+' --npy='+str( n )+' --npz='+str(n)+' '
				#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
				print( exe_pre(n*n*n)+exe_path+exe+case_para+case_consts +' > output ' )
				os.system( exe_pre(n*n*n)+exe_path+exe+case_para+case_consts +' > output ' )
