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

ma.setParameter( root, 'withoutput', 1 )
npx = 8
npy = 2
npz = 4
ma.setParameter( root, 'npx', npx )
ma.setParameter( root, 'npy', npy )
ma.setParameter( root, 'npz', npz )

exe = 'peri_navier3D'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )



case_path = ['','','','']

refinement = [ 4,1 ]
ns         = [ 1,4 ]

ma.setParameter( root, 'alpha2', 2*pi*0.2*200 )

case_path[0] = '/case_study'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
for i in range(2):
	case_path[1] = '/ref_'+str(i)
	if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		os.mkdir( data_path+case_path[0]+case_path[1] )
	os.chdir( data_path+case_path[0]+case_path[1] )
	os.system(' rm ./* -r -v  ')
	ma.setParameter( root, 'refinement level', refinement[i] )
	ma.setParameter( root, 'nx', 256+1 )
	ma.setParameter( root, 'ny',  64+1 )
	ma.setParameter( root, 'nz', 128+1 )
	ma.setParameter( root, 'nf',  ns[i] )
	tree.write('parameter3D.xml')
	os.system( exe_pre(npx*npy*npz,' -W 21:00 ') + exe_path+exe + ' > output ' )
