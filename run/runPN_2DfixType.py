import os
import platform_paths as pp


EXE = 'peri_navier'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

CASE_PATH = ['']*3

case_consts = ' --nfe=17  --nx=33 --ny=33 --tolNOX=1.e-6 --tolBelos=1.e-4 --tolNF=1.e-6 --maxIter=10  --domain=2 --lx=2. --ly=2.  '

flow = 5
re = 100
st = 1000
fixTypes = [1, 2, 3, 4, 5, 6, 7]

for fixType in fixTypes:
    CASE_PATH[0] = 'fixType_'+str(fixType)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0])
    print pp.DATA_PATH + CASE_PATH[0]
    os.chdir(pp.DATA_PATH+CASE_PATH[0])
    #for st in STS:
        #CASE_PATH[1] = '/alpha2_'+str(st)
        #if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
          #os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        #print pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1]
        #os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        #os.system(' rm -rv ./* ')
        #for re in res:
          #CASE_PATH[2] = '/re_'+str(re)
          #if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]):
            #os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
          #print pp.DATA_PATH + CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]
          #os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
          #os.system(' rm -vr ./* ')
    case_para = ' --alpha2='+str(st)+' --re='+str(re)+' '+'  --flow='+str(flow)+' --fixType='+str(fixType)+' '
    print case_consts + case_para
    os.system(pp.exe_pre+pp.EXE_PATH+EXE+case_para+case_consts)
