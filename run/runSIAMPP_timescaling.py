from numpy import linspace
import numpy as np
from math import pi
from platformPaths import *
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


case_path = ['','','','','','','','','']

runs = range( 1 )

res       = [ 100 ]
alpha2s   = [ 0.1 ]
nfs       = range( 1, 16 )
nxs       = [ 1, 2, 4 ]
npxs      = [ 4 ]


case_path[0] = '/scale_time'
mkdir( case_path, 0 )

for re in res:
	case_path[1] = '/re_'+str(re)
	mkdir( case_path, 1 )
	for alpha2 in alpha2s:
		case_path[2] = '/a2_'+str(alpha2)
		mkdir( case_path, 2 )
		for nf in nfs:
			case_path[3] = '/nf_'+str(nf)
			mkdir( case_path, 3 )
			for nx in nxs:
				case_path[4] = '/nx_'+str(nx)
				mkdir( case_path, 4 )
				for npx in npxs:
					case_path[5] = '/npx_'+str(npx)
					mkdir( case_path, 5 )
					for npf in range( 1, nf+2 ):
					# for npf in range( nf+1, nf+2 ):
						case_path[6] = '/npf_'+str(npf)
						mkdir( case_path, 6 )
						#
						chdir( case_path, 6 )
						#
						ma.setParameter( root, 'Re', re )
						ma.setParameter( root, 'alpha2', 2.*pi*alpha2*re )
						ma.setParameter( root, 'nx', 64*nx+1 )
						ma.setParameter( root, 'ny', 16*nx+1 )
						ma.setParameter( root, 'nz', 32*nx+1 )
						ma.setParameter( root, 'nf', nf )
						ma.setParameter( root, 'npx',     npx      )
						ma.setParameter( root, 'npy',  max(npx/2,1))
						ma.setParameter( root, 'npz',     npx      )
						ma.setParameter( root, 'npf',     npf      )
						tree.write( 'parameter3D.xml' )
						nptot = npx*npx*max(npx/2,1)*npf
						for run in runs:
							print()
							print( case_path )
							print(     exe_pre( nptot, ' -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(max(1024*1*nx*nf*npf/nptot,1024))+']" ', run ) + exe_path+'/'+exe  )
							os.system( exe_pre( nptot, ' -N -R "select[model==Opteron8380"] -R "rusage[mem='+str(max(1024*1*nx*nf*npf/nptot,1024))+']" ', run ) + exe_path+'/'+exe  )

