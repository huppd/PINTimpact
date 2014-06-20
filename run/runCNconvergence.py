import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


exe = 'peri_navier4'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']
npx= 2 
npy= 2
npt= 2

case_consts = ' --linSolName="GMRES" --flow=5 --domain=2 --force=0   --npx='+str(npx)+' --npy='+str(npy)+' --npt='+str(npt)+' --tolNOX=1.e-6  --tolBelos=1.e-2  --maxIter=20  --lx=2. --ly=2.  '

precTypes = [ 0 ]
ns        = [ 4, 5, 6, 7, 8, 9   ]
res       = [ 10, 100, 200 ]
alpha2s   = [ 10, 100, 200 ]
fixTypes  = [ 1 ]

data_path = '/cluster/scratch_xl/public/huppd'

#ns  = [ 7 ]
#precTypes = [ 0,1,2 ]
res = [ 10  ]
#res = [ 200 ]
alpha2s = [ 20 ] 
#alpha2s = [ 15**2 ] 
fixTypes  = [ 1 ]

for precType in precTypes:
	case_path[0] = '/1precType_'+str(precType)
	if not os.path.exists( data_path+case_path[0] ):
		os.mkdir( data_path+case_path[0] )
	for fixType in fixTypes:
		case_path[1] = '/fixType_'+str(fixType)
		if not os.path.exists( data_path+case_path[0]+case_path[1] ):
			os.mkdir( data_path+case_path[0]+case_path[1] )
		for re in res:
			case_path[2] = '/re_'+str(re)
			if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
				os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
			for alpha2 in alpha2s:
				case_path[3] = '/alpha2_'+str(alpha2)
				if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] ):
					os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
				for n in ns:
					case_path[4] = '/n2_'+str(n)
					if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4]  ):
						os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4] )
					os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4] )
					os.system(' rm ./* -r -v  ')
					case_para = ' --precType='+str(precType)+' --nx='+str(2**5+1)+' --ny='+str(2**5+1)+' --nt='+str(2**n)+' --re='+str(re)+' --alpha2='+str(alpha2)+' --fixType='+str(fixType)+' '
					print case_consts + case_para
					os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
					#os.system( exe_pre(npx*npy*npt)+exe_path+exe+case_para+case_consts )
