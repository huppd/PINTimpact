import os
from numpy import linspace
import numpy as np
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
npx= 1 
npy= 1
npt= 1


case_consts = [' --baseflow=1 --force=1  --domain=1 --nf=8  --tolNOX=1.e-6  --tolInnerBelos=1.e-1 --maxIter=2 --ly=2. ']

ns       = [ 1, 2, 4 ]
re       = 100
alpha2   = 576


case_path[0] = '/weakCheap'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
os.chdir( data_path+case_path[0] )
os.system(' rm ./* -r -v  ')
#for alpha2 in alpha2s:
	#case_path[1] = '/a2_'+str(alpha2)
	#if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		#os.mkdir( data_path+case_path[0]+case_path[1] )
	#for re in res:
		#case_path[2] = '/re_'+str(re)
		#if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
			#os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
for n in ns:
	case_path[3] = '/n2_'+str(n)
	if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3]   ):
		os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
	os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
	case_para = ' --nx='+str(192*n+1)+' --ny='+str(48+1)+' --nz='+str(96*n+1)+' --npx='+str(n)+' --npy='+str( 2 )+' --npz='+str(n)+' --maxGrids='+str(int(3+log2(1)))+' --re='+str(re)+' --alpha2='+str(alpha2)+' --lx='+str(8*n)+'. --lz='+str(4*n)'. '
	#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
	print( exe_pre(2*n**2,' -R "select[model==Opteron6174"] ')+exe_path+exe+case_para+case_consts[0] +' > output ' )
	for run in range(runs):
		os.system( exe_pre(2*n**2,' -R "select[model==Opteron6174"] ',run)+exe_path+exe+case_para+case_consts[0] +'  ' )
