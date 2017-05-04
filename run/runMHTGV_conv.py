from numpy import linspace
import numpy as np
from math import pi
from platformPaths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameterTGV.xml')
root = tree.getroot()

re = 10.

ma.setParameter( root, 'withoutput', 1 )
ma.setParameter( root, 'initial guess', 'zero' )
# ma.setParameter( root, 'refinement level', 1 )

ma.setParameter( root, 'Re', re )

# make executable ready
exe = 'peri_navier3D'
os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )


case_path = ['','','','','','','','','']

runs = range( 1 )

a2s  = [ 1., 10., 0.1 ]
# a2s  = [ 1. ]
# a2s  = [ 10., 0.1 ]


case_path[0] = '/MHTGV_conv'
mkdir( case_path, 0 )

for a2 in a2s:
	case_path[1] = '/a2_'+str(a2)
	mkdir( case_path, 1 )
	#
	chdir( case_path, 1 )
	#
	ma.setParameter( root, 'alpha2', 2.*pi*a2*re )
	# ma.setParameter( root, 'nx', 64*+1 )
	# ma.setParameter( root, 'ny', 64*+1 )
	# ma.setParameter( root, 'nz', 5 )
	# ma.setParameter( root, 'nf', nf )
	ma.setParameter( root, 'npx', 1 )
	ma.setParameter( root, 'npy', 1 )
	ma.setParameter( root, 'npz', 1   )
	ma.setParameter( root, 'npf', 1 )
	tree.write( 'parameter3D.xml' )
	nptot = 1 
	mem = int( max( 1024, 16*1024/nptot ) )
	for run in runs:
		print()
		print( case_path )
		# print(     exe_pre( nptot, ' -N -R "select[model==Opteron6174"] -R "rusage[mem='+str(1024*1)+']" ', run ) + exe_path+'/'+exe  )
		# os.system( exe_pre( nptot, ' -N -R "select[model==Opteron6174"] -R "rusage[mem='+str(1024*1)+']" ', run ) + exe_path+'/'+exe  )
		exeString = exe_pre( nptot, ' -N -R beta -R "span[ptile=4]" -R "rusage[mem=' +str(mem) + ']" ', run ) + exe_path+'/'+exe
		print( exeString  )
		os.system( exeString )
