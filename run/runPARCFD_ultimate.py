import os
from numpy import linspace
import numpy as np
from math import pi
from platform_paths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3D.xml')
root = tree.getroot()

ma.setParameter( root, 'withoutput', 0 )
npx = 8
npy = 1
npz = 2
ma.setParameter( root, 'npx', npx )
ma.setParameter( root, 'npy', npy )
ma.setParameter( root, 'npz', npz )

# make executable ready
exe = 'peri_navier3D'
os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )


case_path = ['','','','']

ns        = [ 2 ]
res       = [ 10, 100, 200 ]
alpha2s   = [ 0.05, 0.1, 0.2, 0.4 ]


case_path[0] = '/ultimate'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )

for n in ns:
	case_path[1] = '' #'/n_'+str(n)
	#if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		#os.mkdir( data_path+case_path[0]+case_path[1] )
	for re in res:
		case_path[2] = '/re_'+str(re)
		if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2] ):
			os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2] )
		for alpha2 in alpha2s:
			case_path[3] = '/a2_'+str(alpha2)
			if not os.path.exists( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] ):
				os.mkdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
			os.chdir( data_path+case_path[0]+case_path[1]+case_path[2]+case_path[3] )
			os.system(' rm ./* -r -v  ')
			#
			ma.setParameter( root, 'Re', re )
			ma.setParameter( root, 'alpha2', 2.*pi*alpha2*re )
			ma.setParameter( root, 'nx', 128*2+1 )
			ma.setParameter( root, 'ny',  32*2+1 )
			ma.setParameter( root, 'nz',  64*2+1 )
			tree.write('parameter3D.xml')
			os.system( exe_pre(npx*npy*npz,' -W 21:00 ') + exe_path+exe + ' > output ' )

