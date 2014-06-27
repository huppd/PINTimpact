import os
import shutil
import sys
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


# defines amount of processors in x and y direction
mx = 4
my = 1
# set-up
os.chdir( path+'/src' )
sar( '  restart = ', '0', 'usr_config.f90')
sar( '  Re = '  , '200', 'usr_config.f90')
sar( '  freq = ', '0.2', 'usr_config.f90')
sar( '  M1 = '  , '  M1 = 49', 'usr_config.f90')
sar( '  M2 = '  , '  M2 = 193', 'usr_config.f90')
sar( '  L2 = '  , '12', 'usr_config.f90')
sar( '  NB1 = ', str(mx), 'usr_config.f90')
sar( '  NB2 = ', str(my), 'usr_config.f90')
sar( '  dtime_max = ', '1./freq/250', 'usr_config.f90' )
sar( '  dtime0 = ', '1./freq/250', 'usr_config.f90' )
sar( '  CFL = ', '0.999', 'usr_config.f90' ) 
sar( '  periodic_tol ='  , ' 1.e-11', 'usr_config.f90')
sar( '  time_end = '     , '200./freq', 'usr_config.f90')
sar( '  write_stout_yes = '  , '.TRUE.', 'usr_config.f90')
sar( '  log_iteration_yes = '  , '.FALSE.', 'usr_config.f90')
sar( '  dtime_out_vect = '  , '1./freq', 'usr_config.f90')
sar( '  Int_dtime = ', '900000', 'usr_config.f90')

if(	os.system( 'make hpc2' ) != 0 ):
	sys.exit()	

# creates and changes to working director for according processor amount
#proc_path = path+'/work/Re200/full'
proc_path = path+'/work/Re200/star'
try:
	os.mkdir( proc_path )
except:
	pass
os.chdir( proc_path )
# moves binary to working director
shutil.copy2(path+'/prog/impact.exe',proc_path+'/impact.exe')
		
os.chdir( proc_path )

# loops over diverent runs for good time measurement
#for i in range(1):
	#try:
		#os.mkdir( proc_path+'/run'+str(i) )
	#except:
		#pass
	
	#os.chdir( proc_path+'/run'+str(i) )
os.system('bsub -N -n '+str(mx*my)+' -R "select[model==Opteron8380]" -W 14:00 mpirun '+proc_path+'/impact.exe')
#os.system('bsub -N -n '+str(mx*my)+'  -W 12:00 mpirun '+proc_path+'/impact.exe')
###############3
#for i in range(10)
	#try:
		#os.mkdir( proc_path )
	#except:
		#pass
	#os.system('bsub -n 4 mpirun '+proc_path+'/impact.exe')

#os.chdir( path+'/run' )
