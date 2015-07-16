import os
from numpy import linspace
#import numpy as nump
from pylab import pi, log2
from platform_paths import *


#data_path = '/cluster/home04/math/huppd/trunk/data/'
#data_path = '/cluster/scratch_xp/public/huppd'
#data_path = '/cluster/scratch_xp/public/huppd'
exe = 'peri_navier3D'

runs = 4

os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )

case_path = ['','','','','']


case_consts = [' --baseflow=1 --force=1  --domain=0 --nf=4  --tolNOX=1.e-6  --tolInnerBelos=1.e-2 --maxIter=2 --ly=2. ']

nxs = [ 1, 2, 3 ]
nps = [ 1, 2, 3 ]
#nxs = [ 2 ]
#nps = [ 1 ]
re       = 100
alpha2   = 576


case_path[0] = '/weak'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
#os.chdir( data_path+case_path[0] )
#os.system(' rm ./* -r -v  ')
#for alpha2 in alpha2s:
	#case_path[1] = '/a2_'+str(alpha2)
	#if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		#os.mkdir( data_path+case_path[0]+case_path[1] )
for nx in nxs:
	case_path[2] = '/nx_'+str(nx)
	if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
		os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
	for np in nps:
		case_path[3] = '/np_'+str(np)
		if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]   ):
			os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
		os.system(' rm ./* -r -v  ')
		case_para = ' --nx='+str(128*np*nx/2+1)+' --ny='+str(32*np*nx/2+1)+' --nz='+str(64*np*nx/2+1)+' --npx='+str(4*np)+' --npy='+str(np)+' --npz='+str(np*2)+' --maxGrids='+str(int(3+log2(np*nx)))+' --re='+str(re)+' --alpha2='+str(alpha2)+' --lx='+str(8)+'. --lz='+str(4)+'. '
		#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
		print( exe_pre(8*np*np*np,' -R "select[model==Opteron6174"] ')+exe_path+exe+case_para+case_consts[0] +' > output ' )
		for run in range(runs):
			os.system( exe_pre(8*np*np*np,' -R "select[model==Opteron6174"] -W 4:00 ',run)+exe_path+exe+case_para+case_consts[0] +'  ' )
