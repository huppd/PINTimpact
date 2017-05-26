""" platfomrpaths and useful stuff """
import os

# exection path
exe_path = os.path.expanduser("~/Pimpact/release/src_c")  # euler
# exe_path = os.path.expanduser( "~/PImpact/release/src_c" ) #hpc

# data patch
data_path = os.path.expanduser("~/data")


# exection command
def exe_pre(np, ops='', run=0):
    # return( "bsub -n "+str(np)+" "+ops+" -o output"+str(run)+" " )
    # euler
    return "bsub -n "+str(np)+' '+ops+' '+' -oo output'+str(run)+' mpirun '
    # return( "mpirun -n "+str(np)+' ' ) # hpc


def getPath(path, npa):
    full_path = data_path
    for i in range(npa+1):
        full_path += path[i]
    return full_path


def mkdir(path, npa):
    full_path = getPath(path, npa)
    if not os.path.exists(full_path):
        os.mkdir(full_path)


def chdir(path, npa):
    os.chdir(getPath(path, npa))
    # os.system( ' rm ./* -r -v  ' )
