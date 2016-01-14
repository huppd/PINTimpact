import os
from numpy import linspace
import numpy as np
from math import pi
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


case_consts = [' --baseflow=1 --force=1  --domain=0 --nf=4  --tolNOX=1.e-6  --tolInnerBelos=1.e-1 --maxIter=5 --lx=8. --ly=2. --lz=4. --refinement=3 ']

ns        = [ 4 ]
res       = [ 100, 200, 400 ]
alpha2s   = [ 0 ]

#res = [100]
#alpha2s = [576]

n = ns[0]

case_path[0] = '/case_study'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
#os.chdir( data_path+case_path[0] )
for alpha2 in alpha2s:
	case_path[1] = '/a2_'+str(alpha2)
	if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		os.mkdir( data_path+case_path[0]+case_path[1] )
	for re in res:
		case_path[2] = '/re_'+str(re)
		if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
			os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
		#for n in ns:
			#case_path[3] = '/n2_'+str(n)
			#if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]   ):
				#os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.system(' rm ./* -r -v  ')
		case_para = ' --nx='+str(32+1)+' --ny='+str(8+1)+' --nz='+str(16+1)+' --npx='+str(n)+' --npy='+str( max(n/4,1) )+' --npz='+str(max(n/2,1))+' --maxGrids='+str(4)+' --re='+str(re)+' --alpha2='+str(alpha2)+' ' 
		#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
		print( exe_pre(n*max(n/2,1)*max(n/4,1),' -W 21:00 ')+exe_path+exe+case_para+case_consts[0] +' > output ' )
		os.system( exe_pre(int(n*max(n/2,1)*max(n/4,1)),' -W 21:00 ' )+exe_path+exe+case_para+case_consts[0] +' > output ' )

