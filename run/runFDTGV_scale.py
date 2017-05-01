from numpy import linspace
import numpy as np
from math import pi
from platformPaths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3DTime.xml')
root = tree.getroot()

ma.setParameter( root, 'withoutput', 0 )
ma.setParameter( root, 'initial guess', 'disturbed' )
# ma.setParameter( root, 'refinement level', 1 )

# make executable ready
exe = 'peri_navier3DTime'
os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )


case_path = ['','','','','','','','','']

runs = range( 0, 5 )

res = [  1,10,100 ]
a2  = [0.1, 1., 10. ]

re = 1
a2 = 1.

npx = [1,2,2,4,4,4,4]
npy = [1,1,2,2,4,4,4]
npf = [1,1,1,1,1,2,3]

ma.setParameter( root, 'nx', 64*+1 )
ma.setParameter( root, 'ny', 64*+1 )
ma.setParameter( root, 'nz', 7 )

case_path[0] = '/MHTGV_re_'+str(re)+'_a2_'+str(a2)
mkdir( case_path, 0 )

for i in range(len(npx)):
	case_path[1] = '/np_'+str(i)
	mkdir( case_path, 1 )
	#
	chdir( case_path, 1 )
	#
	ma.setParameter( root, 'Re', re )
	ma.setParameter( root, 'alpha2', 2.*pi*a2*re )
	# ma.setParameter( root, 'nf', nf )
	ma.setParameter( root, 'npx', npx[i] )
	ma.setParameter( root, 'npy', npy[i] )
	ma.setParameter( root, 'npz', 1   )
	ma.setParameter( root, 'npf', npf[i] )
	tree.write( 'parameter3D.xml' )
	nptot = npx[i]*npy[i]*npf[i]
	mem = int( max( 1024, 16*1024/nptot ) )
	for run in runs:
		print()
		print( case_path )
		exeString = exe_pre( nptot, ' -N -R beta -R "span[ptile=4]" -R "rusage[mem=' +str(mem) + ']" ', run ) + exe_path+'/'+exe
		print( exeString  )
		os.system( exeString )
