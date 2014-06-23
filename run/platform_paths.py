exe_path = '/home/huppd/workspace/pimpact-repo/release/src/src_c/'
exe_path = '/home/huppd/workspace/pimpact-repo/debug/src/src_c/'
prof_path = '/home/huppd/workspace/pimpact-repo/profile/src/src_c/'
data_path = '/home/huppd/workspace/pimpact-repo/data/'
def exe_pre(n=4,bla=0,run=0):
	return 'mpirun -np '+str(n)+' '
prof_pre = ' valgrind --tool=callgrind --dump-instr=yes --cache-sim=yes --branch-sim=yes '
