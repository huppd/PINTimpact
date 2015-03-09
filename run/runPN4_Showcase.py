import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


exe = 'peri_navier4'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']
npx= 4 
npy= 2
npt= 4

case_consts = ' --linSolName="GMRES" --piccard --flow=1 --domain=1 --force=1 --radius=0.1 --amp=0.1  --npx='+str(npx)+' --npy='+str(npy)+' --npt='+str(npt)+' --tolNOX=1.e-2 --tolBelos=1.e-1  --maxIter=20  --lx=4. --ly=2. --xm='+str(0.25) + '  '

precTypes = [ 0, 10 ]
ns        = [ 4 ]
res       = [ 10, 100, 200 ]
alpha2s   = [ 10, 100, 200 ]
fixTypes  = [ 1, 2, 4, 6, 9, 10]

data_path = '/cluster/scratch_xp/public/huppd/'

ns  = [ 5 ]
precTypes = [ 0,1,2 ]
res = [ 200  ]
alpha2s = [ 251 ] 
fixTypes  = [ 1 ]

for precType in precTypes:
	case_path[0] = 'precType_'+str(precType)
	if not os.path.exists( data_path+case_path[0] ):
		os.mkdir( data_path+case_path[0] )
	for n in ns:
		case_path[1] = '/n2_'+str(n)
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
				for fixType in fixTypes:
					case_path[4] = '/fixType_'+str(fixType)
					if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4]  ):
						os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4] )
					os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4] )
					os.system(' rm ./* -r -v  ')
					case_para = ' --precType='+str(precType)+' --nx='+str(2*2**n+1)+' --ny='+str(2**n+1)+' --nt='+str(2**(n-1)+1)+' --re='+str(re)+' --alpha2='+str(alpha2)+' --fixType='+str(fixType)+' '
					print case_consts + case_para
					os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
