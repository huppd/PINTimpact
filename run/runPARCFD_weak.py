import os
from numpy import linspace
#import numpy as nump
from pylab import pi, log2
from platform_paths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3D.xml')
root = tree.getroot()

ma.setParameter( root, 'withoutput', 0 )


# make executable ready
exe = 'peri_navier3D'
os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )
exe = 'peri_navier3D'

runs = range(0,4)


case_path = ['','','','','']



nxs = [ 1, 2, 3 ]
nps = [ 1, 2, 3 ]
#nxs = [ 2 ]
#nps = [ 1 ]
re       = 200
alpha2   = 0.2 


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
		ma.setParameter( root, 'Re', re )
		ma.setParameter( root, 'alpha2', 2.*pi*alpha2*re )
		ma.setParameter( root, 'nx', 128*nx*np+1 )
		ma.setParameter( root, 'ny',  32*nx*np/2+1 )
		ma.setParameter( root, 'nz',  64*nx*np+1 )
		ma.setParameter( root, 'npx', 4*np )
		ma.setParameter( root, 'npy',   np )
		ma.setParameter( root, 'npz', 2*np )
		tree.write('parameter3D.xml')
		#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
		for run in runs:
			print( exe_pre(8*np*np*np,' -R "select[model==Opteron6174"] ')+exe_path+exe+case_para +' > output ' )
			os.system( exe_pre(8*np*np*np,' -R "select[model==Opteron6174"] -W 4:00 ',run)+exe_path+exe+  )
