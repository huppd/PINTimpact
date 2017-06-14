from numpy import linspace
import numpy as np
from math import pi
from platformPaths import *
import xml.etree.ElementTree as ET
import manipulator as ma


# load parameter file
tree = ET.parse('../XML/parameter3DTime.xml')
root = tree.getroot()

ma.setParameter( root, 'withoutput', 1 )
ma.setParameter( root, 'initial guess', 'exact' )
# ma.setParameter( root, 'refinement level', 1 )

# make executable ready
exe = 'peri_navier3DTime'
os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )


case_path = ['','','','','','','','','']

runs = range( 1 )

res = [ 1, 10 ]
res = [ 10 ]
# a2s = [ 1., 10., 100. ]
# a2s = [ 1., 10., 0.1 ]
a2s = [ 0.1, 1., 10. ]
# a2s = [ 1. ]

# nfs = [ 32, 48 , 64]#, 96, 128 ]
nfs = [ 96, 128, 192 ]

ma.setParameter( root, 'nx', 65 )
ma.setParameter( root, 'ny', 65 )
ma.setParameter( root, 'nz', 5 )

case_path[0] = '/FDTGV_dis'
mkdir( case_path, 0 )

for re in res:
	case_path[1] = '/re_'+str(re)
	mkdir( case_path, 1 )
	for a2 in a2s:
		case_path[2] = '/a2_'+str(a2)
		mkdir( case_path, 2 )
		for nf in nfs:
			case_path[3] = '/nt_'+str(nf)
			mkdir( case_path, 3 )
			#
			chdir( case_path, 3 )
			#
			ma.setParameter( root, 'Re', re )
			ma.setParameter( root, 'alpha2', 2.*pi*a2*re )
			ma.setParameter( root, 'nf', nf )
			ma.setParameter( root, 'npx', 1 )
			ma.setParameter( root, 'npy', 1 )
			ma.setParameter( root, 'npz', 1 )
			ma.setParameter( root, 'npf', 12 )
			tree.write( 'parameter3D.xml' )
			# nptot = npx[i]*npy[i]*npf[i]
			nptot = 12 
			mem = int( max( 1024, 60*1024/nptot ) )
			for run in runs:
				print()
				print( case_path )
				# exeString = exe_pre( nptot, ' -N -R beta -R "span[ptile=4]" -R "rusage[mem=' +str(mem) + ']" -W 1:00', run ) + exe_path+'/'+exe
				exeString = exe_pre( nptot, ' -N  -R "rusage[mem=' +str(mem) + ']" -W 6:00', run ) + exe_path+'/'+exe
				print( exeString  )
				os.system( exeString )
