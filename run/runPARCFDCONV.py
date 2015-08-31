import os
from numpy import linspace
import numpy as np
from math import pi
from platform_paths import *
import xml.etree.ElementTree as ET
import manipulator as ma

tree = ET.parse('../XML/parameter3D.xml')

root = tree.getroot()

exe = 'peri_navier3D'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )



case_path = ['','','','']

lambdas = [ 4 ]

ma.setParameter( root, 'alpha2', 2*pi*0.2*200 )
ma.setParameter( root, 'Tolerance', 1.e-9 )

npx = 8
npy = 2
npz = 4
ma.setParameter( root, 'npx', npx )
ma.setParameter( root, 'npy', npy )
ma.setParameter( root, 'npz', npz )

case_path[0] = '/case_study'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
for lam in lambdas:
	case_path[1] = '/n_'+str(lam)
	if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		os.mkdir( data_path+case_path[0]+case_path[1] )
	os.chdir( data_path+case_path[0]+case_path[1] )
	os.system(' rm ./* -r -v  ')
	#ma.setParameter( root, 'refinement', lam )
	ma.setParameter( root, 'nx', 256+1 )
	ma.setParameter( root, 'ny',  64+1 )
	ma.setParameter( root, 'nz', 128+1 )
	ma.setParameter( root, 'nf',  1 )
	tree.write('parameter3D.xml')
	os.system( exe_pre(npx*npy*npz,' -W 21:00 ') + exe_path+exe + ' > output ' )

