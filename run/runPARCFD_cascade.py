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

exe = 'peri_navier3D'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )



case_path = ['','','','']

refinement = [ 1, 4 ]
ns = [4,1 ]

ma.setParameter( root, 'alpha2', 2*pi*0.2*200 )

case_path[0] = '/case_study'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
for i in range(2):
	case_path[1] = '/n_'+str(i)
	if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		os.mkdir( data_path+case_path[0]+case_path[1] )
	os.chdir( data_path+case_path[0]+case_path[1] )
	os.system(' rm ./* -r -v  ')
	ma.setParameter( root, 'refinement', refinement[i] )
	ma.setParameter( root, 'nx', 32*ns[i]+1 )
	ma.setParameter( root, 'ny',  8*ns[i]+1 )
	ma.setParameter( root, 'nz', 16*ns[i]+1 )
	ma.setParameter( root, 'nf',  4 )
	tree.write('parameter3D.xml')
	os.system( exe_pre(npx*npy*npz,' -W 21:00 ') + exe_path+exe + ' > output ' )
