import os
import time
from platform_paths import *


rep = 10
exe = 'peri_burgers_timing'

os.chdir( exe_path )

#os.system( 'make '+exe+' -j4' )

##case_consts = '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '
#case_consts = ''# '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '

#os.chdir( data_path )
#f = open( 'timing.txt','a' )

#runtime = 1e22
#for i in range(rep):
  #start = time.time()
  #os.system(exe_pre()+exe_path+exe+case_consts )
  #rt = time.time()-start
  #runtime = min(rt,runtime)

#f.write('year:'+str(time.gmtime().tm_year ) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

#f.close()

#exe = 'peri_navier_timing'

#os.chdir( exe_path )

#os.system( 'make '+exe+' -j4' )

##case_consts = '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '
#case_consts = ''# '  --tolNOX=1.e-1 --tol=1.e-6  --maxI=30 --nx=33 --ny=33 --nf=16 --npx=2 --npy=2 '

#os.chdir( data_path )
#f = open( 'timingNaveier.txt','a' )

#runtime = 1e22
#for i in range(rep):
  #start = time.time()
  #os.system(exe_pre()+exe_path+exe+case_consts )
  #rt = time.time()-start
  #runtime = min(rt,runtime)

#f.write('year:'+str(time.gmtime().tm_year ) + ', month: '+str(time.gmtime().tm_mon)+', day: '+str(time.gmtime().tm_mday)+',\truntime: '+str(runtime)+'\n')

#f.close()

exe = 'peri_navier4'
os.system( 'make '+exe+' -j4' )
os.chdir( prof_path )
os.system( 'make '+exe+' -j4' )
print exe_pre()+prof_pre+exe
os.system(exe_pre()+prof_pre+'./'+exe+' --nx=33 --ny=33 --nt=32 --domain=1 --flow=1 --force=1 --precType=0 --maxIter=1 --rightPrec --piccard --tolNOX=0.5 --tolBelos=0.5 --tolSchur=0.1 --tolPrec=0.1 --re=200 --alpha2=250 ')

