import os
from numpy import linspace
import numpy as np
from pylab import pi, log2
from platform_paths import *


#data_path = '/cluster/home04/math/huppd/trunk/data/'
#data_path = '/cluster/scratch_xp/public/huppd'
#data_path = '/cluster/scratch_xp/public/huppd'
exe = 'peri_navier3D'

runs = range(0,4)

os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']
npx= 1 
npy= 1
npt= 1


case_consts = ['--withoutput=0 --baseflow=1 --force=1  --domain=0 --nf=4  --tolNOX=1.e-6  --tolInnerBelos=1.e-1 --maxIter=2 --lx=8. --ly=2. --lz=4.']

nps       = [ 4, 8, 12 ]
nxs       = [ 1, 2, 3 ]
#ns       = [  2, 3, 4, 6, 8,12 ]
#ns       = [ 1, 2 ]
#nps = [8]
#nxs = [2]
#nps       = [ 4, 8 ]
re       = 100
alpha2   = 576


case_path[0] = '/strong2'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
#os.chdir( data_path+case_path[0] )
#os.system(' rm ./* -r -v  ')
#for alpha2 in alpha2s:
	#case_path[1] = '/a2_'+str(alpha2)
	#if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		#os.mkdir( data_path+case_path[0]+case_path[1] )
for n in nxs:
	case_path[2] = '/n_'+str(n)
	if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
		os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
	for np in nps:
		case_path[3] = '/np_'+str(np)
		if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]   ):
			os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.system(' rm ./* -r -v  ')
		case_para = ' --nx='+str(128*n+1)+' --ny='+str(32*n+1)+' --nz='+str(64*n+1)+' --npx='+str(np)+' --npy='+str( max(np/4,1) )+' --npz='+str(max(np/2,1))+' --maxGrids='+str(int(2+n))+' --re='+str(re)+' --alpha2='+str(alpha2)+' ' 
		#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
		print( exe_pre(np*max(np/2,1)*max(np/4,1),' -R "select[model==Opteron8380"] ')+exe_path+exe+case_para+case_consts[0] +'  ' )
		for run in runs:
			os.system( exe_pre(np*max(np/2,1)*max(np/4,1),' -R "select[model==Opteron6174"] -R "rusage[mem=8192]" -W 2:00',run)+exe_path+exe+case_para+case_consts[0] +'  ' )
