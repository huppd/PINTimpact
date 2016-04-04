import os

# exection path
exe_path = "/cluster/home04/math/huppd/PImpact/release/src_c"

# data patch
data_path = "/cluster/home04/math/huppd/data"

# exection command
def exe_pre( np, ops='', run=0 ):
	# return( "bsub -n "+str(np)+" "+ops+" -o output"+str(run)+" " )
	return( "bsub -n "+str(np)+ ' '+ops+' '+' -oo output'+str(run)+' ' )

def mkdir( path, npa ):
	fullPath = data_path
	for i in range( npa+1 ):
		fullPath += path[i]
	if not os.path.exists( fullPath ):
		os.mkdir( fullPath )

def chdir( path, npa ):
	fullPath = data_path
	for i in range( npa+1 ):
		fullPath += path[i]
	os.chdir( fullPath )
	os.system( ' rm ./* -r -v  ' )
