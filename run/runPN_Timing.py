import os
import time
from platform_paths import *


data_path = '/home/huppd/pimpact-repo/trunk/timing/'

rep = 4 
exe = 'peri_navier'

os.chdir( exe_path )

os.system( 'make '+exe+' -j4' )

case_consts = ' --nx=65 --ny=65 --nf=8 --npx='+str(1)+' --npy='+str(1)+' --npf='+str(1)+' --tolNOX=1.e-6  --tolBelos=1.e-4  --maxIter=5  --lx=2. --ly=2.  --initZero=0 --tolInnerBelos=1.e-6 --numCycles=2 --maxGrids=4  --withprec=2 --re=10 --alpha2=10  '

os.chdir( data_path )
f = open( 'timingNavier.txt','a' )

runtime = 1e22
for i in range(rep):
	start = time.time()
	os.system(exe_pre()+exe_path+exe+case_consts )
	rt = time.time()-start
	runtime = min(rt,runtime)

f.write('year:'+str(time.gmtime().tm_year ) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

f.close()

exe = 'peri_navier'
os.system( 'make '+exe+' -j4' )
os.chdir( prof_path )
os.system( 'make '+exe+' -j4' )
print exe_pre()+prof_pre+exe
os.system(exe_pre()+prof_pre+'./'+exe+case_consts )

