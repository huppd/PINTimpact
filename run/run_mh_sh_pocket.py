""" runner for swept Hiemenz flow """
from math import pi
import os
import xml.etree.ElementTree as ET
import manipulator as ma
import platform_paths as pp


# load parameter file
ma.set_ids('../XML/parameterSHLabs.xml')
TREE = ET.parse('../XML/parameterSHLabs.xml')
ROOT = TREE.getroot()

# make executable ready
EXE = 'peri_navier3D'
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


ma.set_parameter(ROOT, 'withoutput', 1)
ma.set_parameter(ROOT, 'refinement step', 1)
ma.set_parameter(ROOT, 'refinement tol', 1.e-4)


STS = [1./60., 1./30., 1./10.]
STS = [1./60., 1./10.]
STS = [1./30.]
st = STS[0]

NF = 1

NYS = [129, 193, 257]
NYS = [513]
# NYS = [129, 145, 193]
# LYS = [1., 1.5, 2.25]
DIRS = ['Y', 'Z']
DIRS = ['Y']
DIRS = ['zero', 'base']
DIRS = ['base']

BACKS = [True, False]
BACKS = [False]

CASE_PATH = ['']*4

CASE_PATH[0] = pp.DATA_PATH + '/IWantItAll'
pp.mkdir(CASE_PATH, 0)
pp.chdir(CASE_PATH, 0)


for di in DIRS:
    CASE_PATH[1] = '/init_' + di
    pp.mkdir(CASE_PATH, 1)
    pp.chdir(CASE_PATH, 1)
    for back in BACKS:
        CASE_PATH[2] = '/backtrack_' + str(back)
        pp.mkdir(CASE_PATH, 2)
        pp.chdir(CASE_PATH, 2)
        # for i, ys in enumerate(LYS):
        for NY in NYS:
            #
            NPX = 1
            NPY = 6
            NPZ = 8
            NPF = 1
            #
            LXO = 22.5
            LYO = 600.
            LZO = 150.
            NXO = 97
            NYO = 1537
            NZO = 513
            #
            NX = 65
            # NY = 129
            NZ = 193
            #
            LX = round(1.5*LXO/(NXO-1)*(NX-1), 1)
            LZ = round(1.5*LZO/(NZO-1)*(NZ-1), 1)
            LY = round(1.5*LYO/(NYO-1)*(NY-1), 1)
            CASE_PATH[3] = '/NY_' + str(NY)
            pp.mkdir(CASE_PATH, 3)
            pp.chdir(CASE_PATH, 3)
            #
            print('LX', LX)
            print('LY', LY)
            print('LZ', LZ)
            #
            ma.set_parameter(ROOT, 'initial guess', di)
            ma.set_parameter(ROOT, 'Re', 300.)
            ma.set_parameter(ROOT, 'alpha2', 2.*pi*st*300.)
            ma.set_parameter(ROOT, 'lx', LX)
            ma.set_parameter(ROOT, 'ly', LY)
            ma.set_parameter(ROOT, 'lz', LZ)
            ma.set_parameter(ROOT, 'origin z', LZ/2.)
            ma.set_parameter(ROOT, 'nx', NX)
            ma.set_parameter(ROOT, 'ny', NY)
            ma.set_parameter(ROOT, 'nz', NZ)
            ma.set_parameter(ROOT, 'nf', NF)
            ma.set_parameter(ROOT, 'max refinement', 5-NF)
            ma.set_parameter(ROOT, 'npx', NPX)
            ma.set_parameter(ROOT, 'npy', NPY)
            ma.set_parameter(ROOT, 'npz', NPZ)
            ma.set_parameter(ROOT, 'npf', NPF)
            if back:
                ma.set_insublist(ROOT, 'Line Search', 'Method', 'Backtrack')
            else:
                ma.set_insublist(ROOT, 'Line Search', 'Method', 'Full Step')
            # ma.set_insublist(ROOT, 'Coarse Grid Solver', 'numIters', 4)
            # ma.set_insublist(ROOT, 'Coarse Grid Solver', 'Jacobi',
                    # True)
            TREE.write('parameter3D.xml')
            nptot = NPX*NPY*NPZ*NPF
            memtot = int(1024.*max(30/nptot, 2))
            print()
            print(CASE_PATH)
            EXE_STRING = pp.exe_pre(nptot, ' -N -W 24:00 ' +
                                    '-R "rusage[mem=' + str(memtot) +
                                    ']" ') + pp.EXE_PATH + '/'+EXE
            print(EXE_STRING)
            os.system(EXE_STRING)
