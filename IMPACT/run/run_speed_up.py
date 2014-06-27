import os
import shutil
import sys
from numpy import array
from seandre import sar



## brutus stuff
#os.system('module purge')
#os.system('module load pgi')
#os.system('module load open_mpi')
#os.system('module load h5_par')


# gets path of the project
os.chdir('..')
path = os.getcwd()


#os.chdir( path+'/work' )
fac = 16 


# defines amount of processors in x and y direction
#mx = [1,1,2,2,4,4,8,8]
#my = [1,2,2,4,4,8,8,16]
if fac==1:
	mx = array([1,1,2,2,4,8, 8])
	my = array([1,2,2,4,4,8,16])
elif fac==2:
	mx = array([1,1,1,2,2, 4, 8])
	my = array([1,2,4,4,8,16,16])
elif fac==4:
	mx = array([1,1,1,1,2, 4, 4])
	my = array([1,2,4,8,8,16,32])
elif fac==8:
	mx = array([1,1,1,1, 1, 1, 2])
	my = array([1,2,4,8,16,64,64])
elif fac==16:
	mx = array([1,1,1,1, 1, 1,  1])
	my = array([1,2,4,8,16,64,128])
	mx = array([1])
	my = array([1])
elif fac==54:
	mx = array([1,1,1,1, 1, 1,  1])
	my = array([1,2,4,8,16,64,128])
#mx = mx/fac
#mx[mx<1] = 1
#my = my*fac

# set-up
os.chdir( path+'/src' )
sar( '  Re = '  , '50', 'usr_config.f90')
sar( '  freq = ', '0.2', 'usr_config.f90')
sar( '  dtime_max = ', '1./freq/100' )
sar( '  dtime0 = ', '1./freq/100' )
sar( '  M1 = '  , '81', 'usr_config.f90')
sar( '  M2 = '  , str(fac*160+1), 'usr_config.f90')
sar( '  L2 = '  , str(fac*12), 'usr_config.f90')
sar( '  dtime_max = ', '1./freq/100', 'usr_config.f90' )
sar( '  dtime0 = ', '1./freq/100', 'usr_config.f90' )
sar( '  CFL = ', '0.75', 'usr_config.f90' ) 
sar( '  periodic_tol = ' , '1.e-6', 'usr_config.f90')
sar( '  time_end = '     , '400./freq', 'usr_config.f90')
sar( '  write_stout_yes = '  , '.FALSE.', 'usr_config.f90')
sar( '  log_iteration_yes = ', '.FALSE.', 'usr_config.f90')
sar( '  dtime_out_vect = '   , '400./freq', 'usr_config.f90')

# loop over the different amounts of proccesors
for m in range(len(mx)):
#for m in range(2):
	# changes usr_config file to according amount of processors
	os.chdir( path+'/src' )
	sar( '  NB1 = ', '  NB1 = '+str(mx[m]), 'usr_config.f90')
	sar( '  NB2 = ', '  NB2 = '+str(my[m]), 'usr_config.f90')
	# compiles changed code
	#if(	os.system( 'make brutus' ) != 0 ):
	if(	os.system( 'make hpc2' ) != 0 ):
		sys.exit()	
	p = mx[m]*my[m]
	print
	print '______________________________________________________________'
	print 'p: ',p
	print '______________________________________________________________'
	print
	# creates and changes to working director for according processor amount
	proc_path = path+'/work/proc_'+str(p)
	try:
		os.mkdir( proc_path )
	except:
		pass
	os.chdir( proc_path )
	# moves binary to working director
	shutil.copy2(path+'/prog/impact.exe',proc_path+'/impact.exe')
	# loops over diverent runs for good time measurement
	for i in range(10):
		try:
			os.mkdir( proc_path+'/run'+str(i) )
		except:
			pass
		
		os.chdir( proc_path+'/run'+str(i) )
		os.system('bsub -N -n '+str(p)+' -R "select[model==Opteron8380]" -W 24:00 mpirun '+proc_path+'/impact.exe')

#os.chdir( path+'/run' )
