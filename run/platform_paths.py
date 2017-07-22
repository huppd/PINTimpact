""" platfomrpaths and useful stuff """
import os

# exection path
EXE_PATH = os.path.expanduser("~/Pimpact/release/src_c")  # euler
# EXE_PATH = os.path.expanduser("~/PImpact/release/src_c")  # hpc

# data patch
# DATA_PATH = os.path.expanduser("~/data") # hpc
DATA_PATH = os.path.expanduser("/cluster/scratch/huppd/data") # euler


# exection command
def exe_pre(npro, ops='', run=0):
    """ return pre execution string """
    # euler
    return "bsub -n "+str(npro)+' '+ops+' '+' -oo output'+str(run)+' mpirun '
    # return "mpirun -n "+str(npro)+' '  # hpc


def get_path(path, npa):
    """ creates path string from list """
    return ''.join(path[:npa+1])


def mkdir(path, npa):
    """ makes director from path list """
    full_path = get_path(path, npa)
    if not os.path.exists(full_path):
        os.mkdir(full_path)


def chdir(path, npa):
    """ changed director path list """
    os.chdir(get_path(path, npa))
    # os.system( ' rm ./* -r -v  ' )
