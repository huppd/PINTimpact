import os
import time
import platform_paths as pp


rep = 10
EXE = 'peri_burgers_timing'

os.chdir(pp.EXE_PATH)

os.system('make '+EXE+' -j4')

#case_consts = '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '
case_consts = ''  # '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '

os.chdir(pp.DATA_PATH)
f = open('timing.txt', 'a')

runtime = 1e22
for i in range(rep):
    start = time.time()
    os.system(pp.exe_pre()+pp.EXE_PATH+EXE+case_consts)
    rt = time.time()-start
    runtime = min(rt, runtime)

f.write('year:'+str(time.gmtime().tm_year) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

f.close()

EXE = 'peri_navier_timing'

os.chdir(pp.EXE_PATH)

os.system('make '+EXE+' -j4')

case_consts = ''  # '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '

os.chdir(pp.DATA_PATH)
f = open('timingNaveier.txt', 'a')

runtime = 1e22
for i in range(rep):
    start = time.time()
    os.system(pp.exe_pre()+pp.EXE_PATH+EXE+case_consts)
    rt = time.time()-start
    runtime = min(rt, runtime)

f.write('year:'+str(time.gmtime().tm_year) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

f.close()

EXE = 'peri_navier4'
os.system('make '+EXE+' -j4')
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')
print pp.exe_pre()+pp.EXE_PATH+EXE
os.system(pp.exe_pre()+pp.EXE_PATH+'./'+EXE+' --nx=33 --ny=33 --nt=32 --domain=1 --flow=1 --force=1 --precType=0 --maxIter=1 --rightPrec --piccard --tolNOX=0.1 --tolBelos=0.1 --tolSchur=0.1 --tolPrec=0.1 --re=200 --alpha2=250 ')
