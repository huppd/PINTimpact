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

runs = range( 10, 12 )

re  = 100 
a2  = 0.1 
nfs = 2**np.arange(1,5) - 1
nfs = [15]
nx  = 2
npx = 4


case_path[0] = '/scale_time_re_'+str(re)+'_a2_'+str(a2)+'_nx_'+str(nx)+'_npx_'+str(npx)
mkdir( case_path, 0 )

for nf in nfs:
	case_path[3] = '/nf_'+str(nf)
	mkdir( case_path, 3 )
	for npf in 2**np.arange(4,5):
		if( npf<=(nf+1) and nf/npf<=15 ):
			case_path[6] = '/npf_'+str(npf)
			mkdir( case_path, 6 )
			#
			chdir( case_path, 6 )
			#
			ma.setParameter( root, 'Re', re )
			ma.setParameter( root, 'alpha2', 2.*pi*a2*re )
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
				# print(     exe_pre( nptot, ' -N -R "select[model==Opteron6174"] -R "rusage[mem='+str(1024*1)+']" ', run ) + exe_path+'/'+exe  )
				# os.system( exe_pre( nptot, ' -N -R "select[model==Opteron6174"] -R "rusage[mem='+str(1024*1)+']" ', run ) + exe_path+'/'+exe  )
				print(     exe_pre( nptot, ' -N -R "select[model==Opteron6174"]  ', run ) + exe_path+'/'+exe  )
				os.system( exe_pre( nptot, ' -N -R "select[model==Opteron6174"]  ', run ) + exe_path+'/'+exe  )