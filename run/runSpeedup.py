import os
from numpy import linspace
import numpy as np
from pylab import pi
from platform_paths import *


exe = 'peri_navier'
runs = 10  


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']

case_consts = ' --linSolName="GMRES" --flow=1 --domain=1 --force=4 --radius=0.1 --rotation=2 --nf=2 --tolNOX=1.e-1 --tolNF=1.e-1 --tolBelos=1.e-6  --maxIter=2  --lx=8. --ly=2.  --fixType=1 --xm='+str(1./8.) + '  '

precTypes = [ 0, 10 ]
ns        = [ 6 ]
res       = [ 10, 100, 200 ]
alpha2s   = [ 10, 100, 200 ]
fixTypes  = [ 1, 2, 4, 6, 9, 10]

ns  = [ 4, 5, 6 ]
ns  = [ 4, 5, 6, 7 ]
precTypes = [ 0 ]
res = [ 100 ]
alpha2s = [ 10**2 ] 
fixTypes  = [ 1 ]
npxs = [ 1, 2, 4, 4, 8, 16, 16, 32 ]
npys = [ 1, 1, 1, 2, 2,  2,  4,  4 ]

npxs = [ 1, 2, 4, 4, 8 ]
npys = [ 1, 1, 1, 2, 2 ]
#npxs = [16]
#npys = [4]

for precType in precTypes:
	case_path[0] = '/speedup2'
	if not os.path.exists( data_path+case_path[0] ):
		os.mkdir( data_path+case_path[0] )
	for n in ns:
		case_path[1] = '/ns_'+str(n)
		if not os.path.exists( data_path+case_path[0]+case_path[1] ):
			os.mkdir( data_path+case_path[0]+case_path[1] )
		for re in res:
			#case_path[2] = '/re_'+str(re)
			case_path[2] = ''
			if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
				os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
			for alpha2 in alpha2s:
				#case_path[3] = '/alpha2_'+str(alpha2)
				case_path[3] = ''
				if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] ):
					os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
				for i in range(len(npxs)):
					case_path[4] = '/np_'+str(npxs[i]*npys[i])
					if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4]  ):
						os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4] )
					os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]+case_path[4] )
					os.system(' rm ./* -r -v  ')
					case_para = ' --precType='+str(precType)+' --nx='+str(4*2**n+1)+' --ny='+str(2**n+1)+' --re='+str(re)+' --alpha2='+str(alpha2)+' --npx='+str(npxs[i])+' --npy='+str(npys[i])
					print exe_pre(npxs[i]*npys[i],' -R "select[model==Opteron8384"] ')+exe_path+exe+case_para+case_consts
					for run in range(runs):
						os.system( exe_pre(npxs[i]*npys[i],' -R "select[model==Opteron8384"] ',run=run)+exe_path+exe+case_para+case_consts )
		  #os.system( exe_pre(npxs[i]*npys[i],run=run )+exe_path+exe+case_para+case_consts )
	  	#os.system( exe_pre(npxs[i]*npys[i],run=run)+exe_path+exe+case_para+case_consts )
