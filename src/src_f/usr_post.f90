!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  subroutine postprocess
  ! (basic subroutine)

  

  use mod_dims
  use mod_vars
  use mod_lib
  use mod_inout
  use usr_vars
  use usr_func
  
  implicit none
  
  
  !--- data post-processing ---
  ! note: - use the following subroutine for post-processing
  !       - the content is just an example/template, feel free to edit and modify it
  !       - currently only subroutines for reading 3D data plus ghost-cell update are available
  !       - available subroutines for 2D can only read (without ghost-cell update)
  !
  integer                ::  m
  character(len=5)       ::  count_char
  character(len=1)       ::  conc_number
  
  
                         open(21,file='post_TKE_mod_restart'       //restart_char//'.txt',status='UNKNOWN')
                         open(22,file='post_TKE_spectrum_restart'  //restart_char//'.txt',status='UNKNOWN')
                         open(23,file='post_turb_stat_restart'     //restart_char//'.txt',status='UNKNOWN')
  if (concentration_yes) open(24,file='post_conc_stat_restart'     //restart_char//'.txt',status='UNKNOWN')
  
  
  do timestep = 1, 50
     call num_to_string(5,timestep,count_char)
     call read2_hdf('velX_med_'//count_char,'velX' ,S1p,S2p,S3p,N1p,N2p,N3p,0,work1)
     call read2_hdf('velY_med_'//count_char,'velY' ,S1p,S2p,S3p,N1p,N2p,N3p,0,work2)
     call read2_hdf('velZ_med_'//count_char,'velZ' ,S1p,S2p,S3p,N1p,N2p,N3p,0,work3)
     do m = 1, n_conc
        write(conc_number,'(i1.1)') m
        call read2_hdf('conc'//conc_number//'_med_'//count_char,'conc1',S1p,S2p,S3p,N1p,N2p,N3p,0,conc(b1L,b2L,b3L,m))
     end do
     call compute_stats
  end do
  
  call close_stats
  
  
  end subroutine postprocess
