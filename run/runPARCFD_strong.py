import os
from numpy import linspace
import numpy as np
from pylab import pi, log2
from platform_paths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3D.xml')
root = tree.getroot()

ma.setParameter( root, 'withoutput', 0 )
ma.setParameter( root, 'refinement level', 1 )


# make executable ready
exe = 'peri_navier3D'
os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )
exe = 'peri_navier3D'


runs = range(0,4)


case_path = ['','','','','']



nps       = [ 4, 8, 12 ]
nxs       = [ 1, 2, 3, 4 ]
#ns       = [  2, 3, 4, 6, 8,12 ]
#ns       = [ 1, 2 ]
#nps = [8]
#nxs = [2]
#nps       = [ 4, 8 ]
re       = 200
alpha2   = 0.2 


case_path[0] = '/strong'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
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
		ma.setParameter( root, 'Re', re )
		ma.setParameter( root, 'alpha2', 2.*pi*alpha2*re )
		ma.setParameter( root, 'nx', 128*n+1 )
		ma.setParameter( root, 'ny',  32*n+1 )
		ma.setParameter( root, 'nz',  64*n+1 )
		ma.setParameter( root, 'npx',     np      )
		ma.setParameter( root, 'npy', max(np/4,1) )
		ma.setParameter( root, 'npz', max(np/2,1) )
		tree.write('parameter3D.xml')
		#os.system( exe_pre(npx*npy*npt,' -R lustre ')+exe_path+exe+case_para+case_consts )
		for run in runs:
			print( exe_pre(np*max(np/2,1)*max(np/4,1),' -R "select[model==Opteron8380"] ')+exe_path+exe )
			os.system( exe_pre(np*max(np/2,1)*max(np/4,1),' -R "select[model==Opteron6174"] -R "rusage[mem=8192]" -W 2:00',run)+exe_path+exe )
