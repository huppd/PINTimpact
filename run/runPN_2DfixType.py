import os
from pylab import pi
from platform_paths import *


exe = 'peri_navier'


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')

CASE_PATH = ['']*3

case_consts = ' --nfe=17  --nx=33 --ny=33 --tolNOX=1.e-6 --tolBelos=1.e-4 --tolNF=1.e-6 --maxIter=10  --domain=2 --lx=2. --ly=2.  '

flow = 5
re = 100
st = 1000
fixTypes = [1, 2, 3, 4, 5, 6, 7]

for fixType in fixTypes:
    CASE_PATH[0] = 'fixType_'+str(fixType)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]):
        os.mkdir(DATA_PATH+CASE_PATH[0])
    print DATA_PATH + CASE_PATH[0]
    os.chdir(DATA_PATH+CASE_PATH[0])
    #for st in STS:
        #CASE_PATH[1] = '/alpha2_'+str(st)
        #if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
          #os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        #print DATA_PATH + CASE_PATH[0] + CASE_PATH[1]
        #os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        #os.system(' rm -rv ./* ')
        #for re in res:
          #CASE_PATH[2] = '/re_'+str(re)
          #if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
            #os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
          #print DATA_PATH + CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]
          #os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
          #os.system(' rm -vr ./* ')
    case_para = ' --alpha2='+str(st)+' --re='+str(re)+' '+'  --flow='+str(flow)+' --fixType='+str(fixType)+' '
    print case_consts + case_para
    os.system(exe_pre+EXE_PATH+exe+case_para+case_consts)
