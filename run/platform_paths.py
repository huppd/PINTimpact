# exection path
exe_path = "/cluster/home04/math/huppd/PImpact/release/src_c"

# data patch
data_path = "/cluster/home04/math/huppd/data"

# exection command
def exe_pre( np, ops ):
	return( "bsub -n "+str(np)+" "+ops )
