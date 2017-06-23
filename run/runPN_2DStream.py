""" runner for Rayleigh streaming """
import os
import platform_paths as pp


EXE = 'peri_navier'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')

CASE_PATH = ['']*3

case_consts = ' --nfe=17  --nx=33 --ny=33 --tolNOX=1.e-6 --tolBelos=1.e-4 --tolNF=1.e-6 --maxIter=10  --domain=2 --lx=2. --ly=2.  '

FLOWS = [5, 6, 7]
RES = [1, 10, 100]
STS = [1, 10, 100, 1000]

for flow in FLOWS:
    CASE_PATH[0] = 'flow_'+str(flow)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0])
    print pp.DATA_PATH + CASE_PATH[0]
    for st in STS:
        CASE_PATH[1] = '/alpha2_'+str(st)
        if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
            os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        print pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1]
        os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
        os.system(' rm -rv ./* ')
        for re in RES:
            CASE_PATH[2] = '/re_'+str(re)
            if not os.path.exists(pp.DATA_PATH + CASE_PATH[0] + CASE_PATH[1] +
                                  CASE_PATH[2]):
                os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            print pp.DATA_PATH + CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2]
            os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]+CASE_PATH[2])
            os.system(' rm -vr ./* ')
            case_para = ' --alpha2='+str(st)+' --re='+str(re)+' '+'  --flow='+str(flow)+' '
            print case_consts + case_para
            os.system(pp.exe_pre+pp.EXE_PATH+EXE+case_para+case_consts)
