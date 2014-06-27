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


# defines frequenzies
Res = ['50','100','200']
freqs = ['0.1','0.2','0.35']
dis = 1

# set-up
os.chdir( path+'/src' )
sar( '  dtime_max = ', '1./freq/100' )
sar( '  dtime0 = ', '1./freq/100' )
sar( '  M1 = '  , str(dis*49), 'usr_config.f90')
sar( '  M2 = '  , str(dis*193), 'usr_config.f90')
sar( '  NB1 = ', str(dis*2), 'usr_config.f90')
sar( '  NB2 = ', str(dis*8), 'usr_config.f90')
sar( '  L2 = '  , '12', 'usr_config.f90')
sar( '  dtime_max = ', '1./freq/100', 'usr_config.f90' )
sar( '  dtime0 = ', '1./freq/100', 'usr_config.f90' )
sar( '  CFL = ', '0.75', 'usr_config.f90' ) 
sar( '  periodic_tol ='  , ' 1.e-26', 'usr_config.f90')
sar( '  time_end = '     , '200./freq', 'usr_config.f90')
sar( '  write_stout_yes = '  , '.TRUE.', 'usr_config.f90')
sar( '  log_iteration_yes = '  , '.TRUE.', 'usr_config.f90')
sar( '  dtime_out_vect = '  , '1./freq', 'usr_config.f90')

for Re in Res:
	os.chdir( path+'/src' )
	sar( '  Re = ', Re, 'usr_config.f90')
	try:
		os.mkdir( path+'/work/Re_'+Re )
	except:
		pass
	for freq in freqs:
		os.chdir( path+'/src' )
		sar( '  freq = ', freq, 'usr_config.f90')
		if(	os.system( 'make hpc2' ) != 0 ):
			sys.exit()	
		
		# creates and changes to working director for according processor amount
		proc_path = path+'/work/Re_'+Re+'/freq_'+freq
		try:
			os.mkdir( proc_path )
		except:
			pass
		os.chdir( proc_path )
		# moves binary to working director
		shutil.copy2(path+'/prog/impact.exe',proc_path+'/impact.exe')
		
		os.chdir( proc_path )
		os.system('bsub -N -W 8:00 -n '+str(16*dis*dis)+' mpirun '+proc_path+'/impact.exe')

#os.chdir( path+'/run' )
