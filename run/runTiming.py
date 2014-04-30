import os
import time
from platform_paths import *


exe = 'peri_burgers_timing'

os.chdir( exe_path )

os.system( 'make '+exe+' -j4' )

#case_consts = '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '
case_consts = ''# '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '

os.chdir( data_path )
f = open( 'timing.txt','a' )

runtime = 1e22
for i in range(10):
  start = time.time()
  os.system(exe_pre+exe_path+exe+case_consts )
  rt = time.time()-start
  runtime = min(rt,runtime)

f.write('year:'+str(time.gmtime().tm_year ) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

f.close()

exe = 'peri_navier_timing'

os.chdir( exe_path )

os.system( 'make '+exe+' -j4' )

#case_consts = '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '
case_consts = ''# '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '

os.chdir( data_path )
f = open( 'timingNaveier.txt','a' )

runtime = 1e22
for i in range(10):
  start = time.time()
  os.system(exe_pre+exe_path+exe+case_consts )
  rt = time.time()-start
  runtime = min(rt,runtime)

f.write('year:'+str(time.gmtime().tm_year ) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

f.close()
