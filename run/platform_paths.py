""" platfomrpaths and useful stuff """
import os
import time

HPC = True
# HPC = False

# exection path
EXE_PATH = os.path.expanduser("~/PImpact/release/src_c")  # hpc
# EXE_PATH = os.path.expanduser("~/PImpact/debug/src_c")  # hpc

# data patch
if HPC:
    DATA_PATH = os.path.expanduser("~/data")  # hpc
else:
    DATA_PATH = "/cluster/scratch/huppd/data"  # euler


# exection command
def exe_pre(npro, ops='', run=0):
    """ return pre execution string """
    # euler
    if HPC:
        return "mpirun -n "+str(npro)+' '  # hpc
    else:
        time.sleep(60)
        return "bsub -n " + str(npro) + ' ' + ops + ' ' + ' -oo output' + \
            str(run) + ' mpirun '



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
