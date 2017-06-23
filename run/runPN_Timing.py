import os
import time
import platform_paths as pp


rep = 4
EXE = 'peri_navier'

os.chdir(pp.EXE_PATH)

os.system('make '+EXE+' -j4')

case_consts = ' --nx=65 --ny=65 --nf=8 --npx='+str(1)+' --npy='+str(1)+' --npf='+str(1)+' --tolNOX=1.e-6  --tolBelos=1.e-4  --maxIter=5  --lx=2. --ly=2.  --initZero=0 --tolInnerBelos=1.e-6 --numCycles=2 --maxGrids=4  --withprec=2 --re=10 --alpha2=10  '

os.chdir(pp.DATA_PATH)
f = open('timingNavier.txt', 'a')

runtime = 1e22
for i in range(rep):
    start = time.time()
    os.system(pp.exe_pre()+pp.EXE_PATH+EXE+case_consts)
    rt = time.time()-start
    runtime = min(rt, runtime)

f.write('year:'+str(time.gmtime().tm_year) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

f.close()

EXE = 'peri_navier'
os.system('make '+EXE+' -j4')
os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')
print pp.exe_pre()+pp.EXE_PATH+EXE
os.system(pp.exe_pre()+pp.EXE_PATH+'./'+EXE+case_consts)
