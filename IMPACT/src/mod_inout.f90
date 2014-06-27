!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* Apr 2012                                                                                                  *
!*************************************************************************************************************

MODULE mod_inout
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_lib
  USE HDF5
  
  
  PRIVATE
  
  PUBLIC write_fields
  PUBLIC write_restart, read_restart
  PUBLIC write_hdf, read_hdf, read2_hdf, write_hdf_velall
  PUBLIC write_2D_hdf, read2_2D_hdf
  PUBLIC write_1D_hdf, read2_1D_hdf
  PUBLIC write_part_hdf, read_part_hdf
  PUBLIC write_part_serial, read_part_serial
  PUBLIC write_stats_hdf_4D, read_stats_hdf_4D
  
  PUBLIC write_hdf_infoREAL, read_hdf_infoREAL
  PUBLIC write_hdf_infoINT , read_hdf_infoINT
  PUBLIC write_hdf_infoLOG , read_hdf_infoLOG
  
  PUBLIC filespace_props
  PUBLIC lambda
  
  
  
  INCLUDE 'mpif.h'
  
  INTEGER(HID_T)                 ::  file_id, plist_id, dset_id
  INTEGER(HID_T)                 ::  filespace, memspace
  INTEGER(HID_T)                 ::  memtypeREAL, memtypeINT
  
  INTEGER(HSIZE_T )              ::  dims_file  (1:3) ! wird global eingefuehrt, um Probleme bei Uebergabe als Argument nach "CALL filespace_props" zu vermeiden.
  INTEGER(HSSIZE_T)              ::  offset_file(1:3)
  
  INTEGER                        ::  i, j, k ! TEST!!! wo brauchts das global?
  
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  SUBROUTINE write_fields
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  CHARACTER(LEN=5)       ::  count_char
  CHARACTER(LEN=1)       ::  conc_number
  
  
  IF (rank == 0) WRITE(*,'(a)') 'writing fields ...'
  
  !===========================================================================================================
  !=== Ausschrieb-Nr. als String fuer File-Namen =============================================================
  !===========================================================================================================
  CALL num_to_string(5,write_count,count_char)
  !===========================================================================================================
  
  
  CALL exchange(1,1,vel(b1L,b2L,b3L,1))
  CALL exchange(2,2,vel(b1L,b2L,b3L,2))
  CALL exchange(3,3,vel(b1L,b2L,b3L,3))
  
  !===========================================================================================================
  !=== Interpolieren / Schreiben =============================================================================
  !===========================================================================================================
  IF (dimens == 3) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
              DO ii = d1L+1, d1U
                 nl(i,j,k,1) = nl(i,j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1) ! TEST!!! Verallgemeinerte Interpolation verwenden? Mit compute_stats teilen?
              END DO
           END DO
        END DO
     END DO
     IF (write_large) CALL write_hdf('velX_large_'//count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,1)) ! TEST!!! velX --> vel1, etc. umbenennen!!
     IF (write_med  ) CALL write_hdf('velX_med_'  //count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,1))
     IF (write_small) CALL write_hdf('velX_small_'//count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2L,b3L,1))
  ELSE
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              bc33(i,j,1) = bc33(i,j,1) + cIup(ii,i)*vel(i+ii,j,k,1)
           END DO
        END DO
     END DO
     CALL write_2D_hdf('velX_'//count_char,'velX',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,1))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
              DO jj = d2L+1, d2U
                 nl(i,j,k,2) = nl(i,j,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
              END DO
           END DO
        END DO
     END DO
     IF (write_large) CALL write_hdf('velY_large_'//count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,2))
     IF (write_med  ) CALL write_hdf('velY_med_'  //count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,2))
     IF (write_small) CALL write_hdf('velY_small_'//count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2L,b3L,2))
  ELSE
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              bc33(i,j,2) = bc33(i,j,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
           END DO
        END DO
     END DO
     CALL write_2D_hdf('velY_'//count_char,'velY',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,2))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,3) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
              DO kk = d3L+1, d3U
                 nl(i,j,k,3) = nl(i,j,k,3) + cIwp(kk,k)*vel(i,j,k+kk,3)
              END DO
           END DO
        END DO
     END DO
     IF (write_large) CALL write_hdf('velZ_large_'//count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,3))
     IF (write_med  ) CALL write_hdf('velZ_med_'  //count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,3))
     IF (write_small) CALL write_hdf('velZ_small_'//count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2L,b3L,3))
  END IF
  !===========================================================================================================
  IF (1 == 2 .AND. dimens == 3) CALL write_hdf_velall('velA_'//count_char,S1p,S2p,S3p,N1p,N2p,N3p,0,nl)
  !===========================================================================================================
  IF (write_lambda2_yes .AND. dimens == 3) THEN
     CALL lambda
     IF (write_large) CALL write_hdf('lamb_large_'//count_char,'lamb',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,res)
     IF (write_med  ) CALL write_hdf('lamb_med_'  //count_char,'lamb',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,res)
     IF (write_small) CALL write_hdf('lamb_small_'//count_char,'lamb',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,res)
  END IF
  !===========================================================================================================
  IF (dimens == 3) THEN
     IF (write_large) CALL write_hdf('pre_large_'//count_char,'pre',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,pre)
     IF (write_med  ) CALL write_hdf('pre_med_'  //count_char,'pre',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,pre)
     IF (write_small) CALL write_hdf('pre_small_'//count_char,'pre',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,pre)
  ELSE
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,1) = pre(i,j,k)
        END DO
     END DO
     CALL write_2D_hdf('pre_'//count_char,'pre',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,1))
  END IF
  !===========================================================================================================
  IF (concentration_yes) THEN
     IF (dimens == 3) THEN
        DO m = 1, n_conc
           WRITE(conc_number,'(i1.1)') m
           IF (write_large) CALL write_hdf('conc'//conc_number//'_large_'//count_char,'conc'//conc_number,S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,conc(b1L,b2L,b3L,m))
           IF (write_med  ) CALL write_hdf('conc'//conc_number//'_med_'  //count_char,'conc'//conc_number,S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,conc(b1L,b2L,b3L,m))
           IF (write_small) CALL write_hdf('conc'//conc_number//'_small_'//count_char,'conc'//conc_number,S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,conc(b1L,b2L,b3L,m))
           
           CALL write_2D_hdf('depo'//conc_number//'_'//count_char,'depo'//conc_number,N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_conc(1,1,m))
        END DO
     ELSE
        k = 1
        DO m = 1, n_conc
           WRITE(conc_number,'(i1.1)') m
           DO j = S2p, N2p
              DO i = S1p, N1p
                 bc33(i,j,1) = conc(i,j,k,m)
              END DO
           END DO
           CALL write_2D_hdf('conc'//conc_number//'_'//count_char,'conc'//conc_number,N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,1))
        END DO
     END IF
  END IF
  !===========================================================================================================
  IF (particles_yes) THEN
     CALL write_part_hdf('part_'//count_char,.FALSE.)
     IF (dimens == 3) THEN
        ! TEST!!! Schmutziger Umweg:
        res(S1p,1:N2,1:N3) = dep1L_part(1:N2,1:N3)
        CALL exchange_part(2,0,1,res)
        CALL exchange_part(3,0,1,res)
        dep1L_part(1:N2,1:N3) = res(S1p,1:N2,1:N3)
        
        IF (BC_2L == 0 .AND. ls2 ==  0) dep1L_part(1 ,1:N3) = 0.
        IF (BC_2U == 0 .AND. ls2 == -1) dep1L_part(N2,1:N3) = 0.
        IF (BC_3L == 0 .AND. ls3 ==  0) dep1L_part(1:N2,1 ) = 0.
        IF (BC_3U == 0 .AND. ls3 == -1) dep1L_part(1:N2,N3) = 0.
        
        CALL write_2D_hdf('depo_part_'//count_char,'depo_part',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_part(1,1)) ! TEST!!! ok?
     END IF
  END IF
  !===========================================================================================================
  
  write_count    = write_count   + 1
  time_out_vect  = time_out_vect + dtime_out_vect
  write_out_vect = .FALSE.
  
  
  END SUBROUTINE write_fields
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_restart
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  CHARACTER(LEN=3)       ::  next_restart_char
  CHARACTER(LEN=1)       ::  conc_number
  
  
  IF (write_restart_yes) THEN
     
     IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing data for restart',restart,' ...'
     
     !========================================================================================================
     !=== neue Restart-Nr. als String fuer File-Namen ========================================================
     !========================================================================================================
     CALL num_to_string(3,restart,next_restart_char)
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Schreiben ==========================================================================================
     !========================================================================================================
                      CALL write_hdf('velX_restart'//next_restart_char,'velX_restart',S11B,S21B,S31B,N11B,N21B,N31B,1,(/1,1,1/),vel(b1L,b2L,b3L,1))
                      CALL write_hdf('velY_restart'//next_restart_char,'velY_restart',S12B,S22B,S32B,N12B,N22B,N32B,2,(/1,1,1/),vel(b1L,b2L,b3L,2))
     IF (dimens == 3) CALL write_hdf('velZ_restart'//next_restart_char,'velZ_restart',S13B,S23B,S33B,N13B,N23B,N33B,3,(/1,1,1/),vel(b1L,b2L,b3L,3))
     !--------------------------------------------------------------------------------------------------------
     CALL write_hdf('pre_restart'//next_restart_char,'pre_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),pre)
     !--------------------------------------------------------------------------------------------------------
     IF (concentration_yes) THEN
        DO m = 1, n_conc
           WRITE(conc_number,'(i1.1)') m
           CALL write_hdf('conc'//conc_number//'_restart'//next_restart_char,'conc'//conc_number//'_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),conc(b1L,b2L,b3L,m))
           
           ! TEST!!! noch nicht vollstaendig! (betrifft auch andere Routinen)
           IF (dimens == 3) THEN ! TEST!!! "sed" und "depo" koennen in 2D nicht geschrieben werden!!
              CALL write_2D_hdf('depo'//conc_number//'_restart'//next_restart_char,'depo'//conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_conc(1,1,m))
              CALL write_2D_hdf('sed' //conc_number//'_restart'//next_restart_char,'sed' //conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,conc1L    (1,1,m))
           END IF
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (particles_yes) THEN
        !CALL write_part_hdf('part_restart'//next_restart_char,.TRUE.) ! TEST!!!
        CALL write_part_serial('part_restart'//next_restart_char)
        
        IF (dimens == 3) THEN
           ! TEST!!! Schmutziger Umweg:
           res(S1p,1:N2,1:N3) = dep1L_part(1:N2,1:N3)
           CALL exchange_part(2,0,1,res)
           CALL exchange_part(3,0,1,res)
           dep1L_part(1:N2,1:N3) = res(S1p,1:N2,1:N3)
           
           IF (BC_2L == 0 .AND. ls2 ==  0) dep1L_part(1 ,1:N3) = 0.
           IF (BC_2U == 0 .AND. ls2 == -1) dep1L_part(N2,1:N3) = 0.
           IF (BC_3L == 0 .AND. ls3 ==  0) dep1L_part(1:N2,1 ) = 0.
           IF (BC_3U == 0 .AND. ls3 == -1) dep1L_part(1:N2,N3) = 0.
           
           CALL write_2D_hdf('depo_part_restart'//next_restart_char,'depo_part_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_part(1,1)) ! TEST!!! ok?
        END IF
     END IF
     !========================================================================================================
     
  END IF
  
  
  END SUBROUTINE write_restart
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_restart
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk

  REAL                   ::  old_dtime_out_scal
  REAL                   ::  old_dtime_out_vect

  CHARACTER(LEN=1)       ::  conc_number
  
  
  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading data for restart',restart,' ...'
  IF (rank == 0) WRITE(*,*)
  
  !-----------------------------------------------------------------------------------------------------------
  ! Anmerkung: Nur fuer den Notfall ...
  IF (1 == 2) THEN
  
  conc = 0.
  vel  = 0.
  
  CALL read2_hdf('conc1_large_0174','conc1',S1p,S2p,S3p,N1p,N2p,N3p,0,conc(b1L,b2L,b3L,1))
  CALL read2_hdf('velX_large_0174','velX',S1p,S2p,S3p,N1p,N2p,N3p,0,nl(b1L,b2L,b3L,1))
  CALL read2_hdf('velY_large_0174','velY',S1p,S2p,S3p,N1p,N2p,N3p,0,nl(b1L,b2L,b3L,2))
  CALL read2_hdf('velZ_large_0174','velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,nl(b1L,b2L,b3L,3))
  
  CALL exchange(1,1,nl(b1L,b2L,b3L,1))
  CALL exchange(2,2,nl(b1L,b2L,b3L,2))
  CALL exchange(3,3,nl(b1L,b2L,b3L,3))
  
  IF (1 == 1) THEN
  DO k = S31, N31
     DO j = S21, N21
        DO i = S11, N11
           vel(i,j,k,1) = cIpu(g1L,i)*nl(i+g1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              vel(i,j,k,1) = vel(i,j,k,1) + cIpu(ii,i)*nl(i+ii,j,k,1)
           END DO
        END DO
     END DO
  END DO
  DO k = S32, N32
     DO j = S22, N22
        DO i = S12, N12
           vel(i,j,k,2) = cIpv(g2L,j)*nl(i,j+g2L,k,2)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              vel(i,j,k,2) = vel(i,j,k,2) + cIpv(jj,j)*nl(i,j+jj,k,2)
           END DO
        END DO
     END DO
  END DO
  DO k = S33, N33
     DO j = S23, N23
        DO i = S13, N13
           vel(i,j,k,3) = cIpw(g3L,k)*nl(i,j,k+g3L,3)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              vel(i,j,k,3) = vel(i,j,k,3) + cIpw(kk,k)*nl(i,j,k+kk,3)
           END DO
        END DO
     END DO
  END DO
  END IF
  
  ELSE
  
                   CALL read2_hdf('velX_restart'//restart_char,'velX_restart',S11B,S21B,S31B,N11B,N21B,N31B,1,vel(b1L,b2L,b3L,1))
                   CALL read2_hdf('velY_restart'//restart_char,'velY_restart',S12B,S22B,S32B,N12B,N22B,N32B,2,vel(b1L,b2L,b3L,2))
  IF (dimens == 3) CALL read2_hdf('velZ_restart'//restart_char,'velZ_restart',S13B,S23B,S33B,N13B,N23B,N33B,3,vel(b1L,b2L,b3L,3))
  !-----------------------------------------------------------------------------------------------------------
  CALL read2_hdf('pre_restart'//restart_char,'pre_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,pre)
  !-----------------------------------------------------------------------------------------------------------
  IF (concentration_yes) THEN
     DO m = 1, n_conc_old
        WRITE(conc_number,'(i1.1)') m
        CALL read2_hdf('conc'//conc_number//'_restart'//restart_char,'conc'//conc_number//'_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,conc(b1L,b2L,b3L,m))
        
        IF (dimens == 3) THEN ! TEST!!! "sed" und "depo" koennen in 2D nicht geschrieben werden!!
           CALL read2_2D_hdf('depo'//conc_number//'_restart'//restart_char,'depo'//conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_conc(1,1,m))
           CALL read2_2D_hdf('sed' //conc_number//'_restart'//restart_char,'sed' //conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,conc1L    (1,1,m))
        END IF
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (particles_yes) THEN
     !CALL read_part_hdf('part_restart'//restart_char,.TRUE.) ! TEST!!!
     CALL read_part_serial('part_restart'//restart_char)
     IF (dimens == 3) THEN
        dep1L_part = 0.
        CALL read2_2D_hdf('depo_part_restart'//restart_char,'depo_part_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_part(1,1)) ! TEST!!! ok?
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------

         ! --- mjohn 070213 - allow for possibility of adaptation of dtime_out_vect or dtime_out_scal by config.txt after restart
         ! compute next time_out_vect and next time_out_scal by subtracting last dtime (stored in restart file) and adding new dtime (read from config.txt)

         ! Setup file access property list with parallel I/O access
         CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
         CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

         ! Open file collectively
         CALL h5fopen_f('velX_restart'//restart_char//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
         CALL h5pclose_f(plist_id,herror)

         CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
         CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
         !-----------------------------------------------------------------------------------------------------------
         ! Open the dataset:
         CALL h5dopen_f(file_id,'velX_restart',dset_id,herror)
         CALL read_hdf_infoREAL(1,.TRUE. ,.TRUE.,'dtime_out_scal' ,scalar=old_dtime_out_scal)
         CALL read_hdf_infoREAL(1,.TRUE. ,.TRUE.,'dtime_out_vect'  ,scalar=old_dtime_out_vect)
         CALL h5dclose_f(dset_id,herror)
         !===========================================================================================================

         time_out_vect = (time_out_vect-old_dtime_out_vect) + dtime_out_vect
         time_out_scal = (time_out_scal-old_dtime_out_scal) + dtime_out_scal

         ! catch values less than or equal to 'time'
         IF(time_out_vect <= time) THEN
            time_out_vect = time + dtime
         END IF
         IF(time_out_scal <= time) THEN
            time_out_scal = time + dtime
         END IF


  END IF
  
  
  END SUBROUTINE read_restart
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,stride,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname
  
  INTEGER, INTENT(IN)    ::  SS1
  INTEGER, INTENT(IN)    ::  SS2
  INTEGER, INTENT(IN)    ::  SS3
  
  INTEGER, INTENT(IN)    ::  NN1
  INTEGER, INTENT(IN)    ::  NN2
  INTEGER, INTENT(IN)    ::  NN3
  
  INTEGER, INTENT(IN)    ::  vel_dir
  INTEGER, INTENT(IN)    ::  stride(1:3)
  
  REAL   , INTENT(IN)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER(HSIZE_T )      ::  stride_mem(1:3)
  INTEGER(HSIZE_T )      ::  dims_mem  (1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)      ::  offset_mem(1:3)
  
  INTEGER                ::  S1w, M1w, dim1
  INTEGER                ::  S2w, M2w, dim2
  INTEGER                ::  S3w, M3w, dim3
  
  LOGICAL                ::  attr_yes
  
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! Nur fuer Druckgitter! (sonst wird es richtig kompliziert, siehe Multigrid fuer Helmholtz-Problem!).
  IF (vel_dir /= 0) THEN
     IF (stride(vel_dir) /= 1) THEN
        IF (rank == 0) WRITE(*,*) 'ERROR! Cannot write velocity field with stride /= 1 in the corresponding direction!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
  END IF
  IF ((MOD((N1-1),stride(1)) /= 0) .OR. (MOD((N2-1),stride(2)) /= 0) .OR. (MOD((N3-1),stride(3)) /= 0)) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot write field with this stride!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  CALL filespace_props(vel_dir,stride,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  dims_data  = (/((NN1-SS1)/stride(1)+1),((NN2-SS2)/stride(2)+1),((NN3-SS3)/stride(3)+1)/)
  offset_mem = (/(SS1-b1L),(SS2-b2L),(SS3-b3L)/)
  stride_mem = (/stride(1),stride(2),stride(3)/)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1  M2  M3 '      ,array =(/M1 ,M2 ,M3 /)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'S1w S2w S3w'      ,array =(/S1w,S2w,S3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1w M2w M3w'      ,array =(/M1w,M2w,M3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'NB1 NB2 NB3'      ,array =(/NB1,NB2,NB3/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'ls1 ls2 ls3'      ,array =(/ls1,ls2,ls3/)  )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'concentration_yes',scalar=concentration_yes)
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'n_conc'           ,scalar=n_conc           )
  CALL write_hdf_infoREAL(3     ,.FALSE.,attr_yes,'gravity'          ,array =gravity          )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'Re'               ,scalar=Re               )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'Sc'               ,array =Sc               )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'Ric'              ,array =Ric              )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'usc'              ,array =usc              )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'CFL   '           ,scalar=CFL              )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'thetaL'           ,scalar=thetaL           )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'epsU  '           ,scalar=epsU             )
  
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'Euler_yes'        ,scalar=Euler_yes        )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'twostep_yes'      ,scalar=twostep_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'mapping_yes'      ,scalar=mapping_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'upwind_yes'       ,scalar=upwind_yes       )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'upwind_conc_yes'  ,scalar=upwind_conc_yes  )
  
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iL',array=(/BC_1L_global,BC_2L_global,BC_3L_global/))
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iU',array=(/BC_1U_global,BC_2U_global,BC_3U_global/))
  
  CALL write_hdf_infoINT (dim_ncb1c,.FALSE.,attr_yes,'ncb1c',array=ncb1c    )
  CALL write_hdf_infoINT (dim_ncb1g,.FALSE.,attr_yes,'ncb1g',array=ncb1g    )
  CALL write_hdf_infoINT (dim_ncb1d,.FALSE.,attr_yes,'ncb1d',array=ncb1d    )
  
  CALL write_hdf_infoINT (dim_ncb2c,.FALSE.,attr_yes,'ncb2c',array=ncb2c    )
  CALL write_hdf_infoINT (dim_ncb2g,.FALSE.,attr_yes,'ncb2g',array=ncb2g    )
  CALL write_hdf_infoINT (dim_ncb2d,.FALSE.,attr_yes,'ncb2d',array=ncb2d    )
  
  CALL write_hdf_infoINT (dim_ncb3c,.FALSE.,attr_yes,'ncb3c',array=ncb3c    )
  CALL write_hdf_infoINT (dim_ncb3g,.FALSE.,attr_yes,'ncb3g',array=ncb3g    )
  CALL write_hdf_infoINT (dim_ncb3d,.FALSE.,attr_yes,'ncb3d',array=ncb3d    )
  
  CALL write_hdf_infoREAL( M1      ,.FALSE.,attr_yes,'y1p'  ,array=y1p(1:M1))
  CALL write_hdf_infoREAL( M2      ,.FALSE.,attr_yes,'y2p'  ,array=y2p(1:M2))
  CALL write_hdf_infoREAL( M3      ,.FALSE.,attr_yes,'y3p'  ,array=y3p(1:M3))
  
  CALL write_hdf_infoREAL((M1+1)   ,.FALSE.,attr_yes,'y1u'  ,array=y1u(0:M1))
  CALL write_hdf_infoREAL((M2+1)   ,.FALSE.,attr_yes,'y2v'  ,array=y2v(0:M2))
  CALL write_hdf_infoREAL((M3+1)   ,.FALSE.,attr_yes,'y3w'  ,array=y3w(0:M3))
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror,stride_mem)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  IF (vel_dir == 1) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1u(S1w::stride(1)))
  ELSE
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w::stride(1)))
  END IF
  IF (vel_dir == 2) THEN
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2v(S2w::stride(2)))
  ELSE
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2p(S2w::stride(2)))
  END IF
  IF (vel_dir == 3) THEN
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3w(S3w::stride(3)))
  ELSE
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3p(S3w::stride(3)))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_hdf_velall(filename,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  
  INTEGER, INTENT(IN)      ::  SS1
  INTEGER, INTENT(IN)      ::  SS2
  INTEGER, INTENT(IN)      ::  SS3
  
  INTEGER, INTENT(IN)      ::  NN1
  INTEGER, INTENT(IN)      ::  NN2
  INTEGER, INTENT(IN)      ::  NN3
  
  INTEGER, INTENT(IN)      ::  vel_dir
  
  REAL   , INTENT(IN)      ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  INTEGER(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)        ::  offset_mem(1:3)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  INTEGER                  ::  S3w, M3w, dim3
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  CALL filespace_props(vel_dir,(/1,1,1/),S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  dims_data  = (/(NN1-SS1+1   ),(NN2-SS2+1   ),(NN3-SS3+1   )/)
  offset_mem = (/(SS1-b1L     ),(SS2-b2L     ),(SS3-b3L     )/)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'velX',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1  M2  M3 '      ,array =(/M1 ,M2 ,M3 /)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'S1w S2w S3w'      ,array =(/S1w,S2w,S3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1w M2w M3w'      ,array =(/M1w,M2w,M3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'NB1 NB2 NB3'      ,array =(/NB1,NB2,NB3/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'ls1 ls2 ls3'      ,array =(/ls1,ls2,ls3/)  )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'concentration_yes',scalar=concentration_yes)
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'n_conc'           ,scalar=n_conc           )
  CALL write_hdf_infoREAL(3     ,.FALSE.,attr_yes,'gravity'          ,array =gravity          )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'Re'               ,scalar=Re               )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'Sc'               ,array =Sc               )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'Ric'              ,array =Ric              ) ! TEST!!! Rip und usp auch beruecksichtigen!
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'usc'              ,array =usc              )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'CFL   '           ,scalar=CFL              )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'thetaL'           ,scalar=thetaL           )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'epsU  '           ,scalar=epsU             )
  
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'Euler_yes'        ,scalar=Euler_yes        )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'twostep_yes'      ,scalar=twostep_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'mapping_yes'      ,scalar=mapping_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'upwind_yes'       ,scalar=upwind_yes       )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'upwind_conc_yes'  ,scalar=upwind_conc_yes  )
  
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iL',array=(/BC_1L_global,BC_2L_global,BC_3L_global/))
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iU',array=(/BC_1U_global,BC_2U_global,BC_3U_global/))
  
  CALL write_hdf_infoINT (dim_ncb1c,.FALSE.,attr_yes,'ncb1c',array=ncb1c    )
  CALL write_hdf_infoINT (dim_ncb1g,.FALSE.,attr_yes,'ncb1g',array=ncb1g    )
  CALL write_hdf_infoINT (dim_ncb1d,.FALSE.,attr_yes,'ncb1d',array=ncb1d    )
  
  CALL write_hdf_infoINT (dim_ncb2c,.FALSE.,attr_yes,'ncb2c',array=ncb2c    )
  CALL write_hdf_infoINT (dim_ncb2g,.FALSE.,attr_yes,'ncb2g',array=ncb2g    )
  CALL write_hdf_infoINT (dim_ncb2d,.FALSE.,attr_yes,'ncb2d',array=ncb2d    )
  
  CALL write_hdf_infoINT (dim_ncb3c,.FALSE.,attr_yes,'ncb3c',array=ncb3c    )
  CALL write_hdf_infoINT (dim_ncb3g,.FALSE.,attr_yes,'ncb3g',array=ncb3g    )
  CALL write_hdf_infoINT (dim_ncb3d,.FALSE.,attr_yes,'ncb3d',array=ncb3d    )
  
  CALL write_hdf_infoREAL( M1      ,.FALSE.,attr_yes,'y1p'  ,array=y1p(1:M1))
  CALL write_hdf_infoREAL( M2      ,.FALSE.,attr_yes,'y2p'  ,array=y2p(1:M2))
  CALL write_hdf_infoREAL( M3      ,.FALSE.,attr_yes,'y3p'  ,array=y3p(1:M3))
  
  CALL write_hdf_infoREAL((M1+1)   ,.FALSE.,attr_yes,'y1u'  ,array=y1u(0:M1))
  CALL write_hdf_infoREAL((M2+1)   ,.FALSE.,attr_yes,'y2v'  ,array=y2v(0:M2))
  CALL write_hdf_infoREAL((M3+1)   ,.FALSE.,attr_yes,'y3w'  ,array=y3w(0:M3))
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,1),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'velY',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,2),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'velZ',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,3),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (1 == 1) THEN ! TEST!!! wird u.a. gebraucht, um 'velabs' zu speichern.
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'pre',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,pre,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (vel_dir == 1) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1u(S1w:M1w))
  ELSE
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w:M1w))
  END IF
  IF (vel_dir == 2) THEN
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2v(S2w:M2w))
  ELSE
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2p(S2w:M2w))
  END IF
  IF (vel_dir == 3) THEN
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3w(S3w:M3w))
  ELSE
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3p(S3w:M3w))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_hdf_velall
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! evtl. schoener: write_hdf_2D??
  ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
  SUBROUTINE write_2D_hdf(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname
  
  INTEGER, INTENT(IN)      ::  N1
  INTEGER, INTENT(IN)      ::  N2
  
  INTEGER, INTENT(IN)      ::  SS1
  INTEGER, INTENT(IN)      ::  SS2
  
  INTEGER, INTENT(IN)      ::  NN1
  INTEGER, INTENT(IN)      ::  NN2
  
  INTEGER, INTENT(IN)      ::  iShift
  INTEGER, INTENT(IN)      ::  jShift
  
  INTEGER, INTENT(IN)      ::  dir
  
  REAL   , INTENT(IN)      ::  phi(1:N1,1:N2)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  
  LOGICAL                  ::  attr_yes
  
  INTEGER(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
  INTEGER(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  IF (ABS(dir) == 1) THEN
     S1w = 2  + ls2
     S2w = 2  + ls3
     
     M1w = M2 + ls2
     M2w = M3 + ls3
     
     IF (BC_2L_global /= -1) THEN
        S1w = 1
        M1w = M2
     END IF
     IF (BC_3L_global /= -1) THEN
        S2w = 1
        M2w = M3
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (ABS(dir) == 2) THEN
     S1w = 2  + ls1
     S2w = 2  + ls3
     
     M1w = M1 + ls1
     M2w = M3 + ls3
     
     IF (BC_1L_global /= -1) THEN
        S1w = 1
        M1w = M1
     END IF
     IF (BC_3L_global /= -1) THEN
        S2w = 1
        M2w = M3
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (ABS(dir) == 3) THEN
     S1w = 2  + ls1
     S2w = 2  + ls2
     
     M1w = M1 + ls1
     M2w = M2 + ls2
     
     IF (BC_1L_global /= -1) THEN
        S1w = 1
        M1w = M1
     END IF
     IF (BC_2L_global /= -1) THEN
        S2w = 1
        M2w = M2
     END IF
  END IF
  !===========================================================================================================
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  
  dims_file   = (/ dim1 , dim2 /)
  dims_mem    = (/ N1   , N2   /)
  dims_data   = (/(NN1-SS1+1),(NN2-SS2+1)/)
  offset_mem  = (/ SS1-1, SS2-1/)
  offset_file = (/iShift,jShift/)
  
  IF (.NOT. ((dir == -1 .AND. iB(1,1) == 1  ) .OR. (dir == -2 .AND. iB(2,1) == 1  ) .OR. (dir == -3 .AND. iB(3,1) == 1   ) .OR.    &
        &    (dir ==  1 .AND. iB(1,1) == NB1) .OR. (dir ==  2 .AND. iB(2,1) == NB2) .OR. (dir ==  3 .AND. iB(3,1) == NB3))) THEN
     dims_data   = (/0,0/)
     offset_mem  = (/0,0/)
     offset_file = (/0,0/)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(2,dims_file,filespace,herror)
  CALL h5screate_simple_f(2,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
  CALL write_hdf_infoINT (2     ,.FALSE.,attr_yes,'S1w S2w'          ,array =(/S1w,S2w/)      )
  CALL write_hdf_infoINT (2     ,.FALSE.,attr_yes,'M1w M2w'          ,array =(/M1w,M2w/)      )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'concentration_yes',scalar=concentration_yes)
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'n_conc'           ,scalar=n_conc           )
  CALL write_hdf_infoREAL(3     ,.FALSE.,attr_yes,'gravity'          ,array =gravity          )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'Re'               ,scalar=Re               )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'Sc'               ,array =Sc               )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'Ric'              ,array =Ric              )
  CALL write_hdf_infoREAL(n_conc,.FALSE.,attr_yes,'usc'              ,array =usc              )
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  IF (ABS(dir) == 1) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorY',array=y2p(S1w:M1w))
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorZ',array=y3p(S2w:M2w))
  END IF
  IF (ABS(dir) == 2) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w:M1w))
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorZ',array=y3p(S2w:M2w))
  END IF
  IF (ABS(dir) == 3) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w:M1w))
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2p(S2w:M2w))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_2D_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_1D_hdf(filename,dsetname,NN,SSp,NNp,iShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname
  
  INTEGER, INTENT(IN)      ::  NN
  INTEGER, INTENT(IN)      ::  SSp
  INTEGER, INTENT(IN)      ::  NNp
  INTEGER, INTENT(IN)      ::  iShift
  INTEGER, INTENT(IN)      ::  dir
  REAL   , INTENT(IN)      ::  phi(1:NN)
  
  INTEGER                  ::  SSw, MMw, dims
  
  LOGICAL                  ::  attr_yes
  
  INTEGER(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
  INTEGER(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  IF (dir == 1) THEN
     SSw = 2  + ls1
     MMw = M1 + ls1
     IF (BC_1L_global /= -1) THEN
        SSw = 1
        MMw = M1
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dir == 2) THEN
     SSw = 2  + ls2
     MMw = M2 + ls2
     IF (BC_2L_global /= -1) THEN
        SSw = 1
        MMw = M2
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dir == 3) THEN
     SSw = 2  + ls3
     MMw = M3 + ls3
     IF (BC_3L_global /= -1) THEN
        SSw = 1
        MMw = M3
     END IF
  END IF
  !===========================================================================================================
  dims = MMw-SSw+1
  
  dims_file   = (/dims     /)
  dims_mem    = (/NN       /)
  dims_data   = (/NNp-SSp+1/)
  offset_mem  = (/SSp-1    /)
  offset_file = (/iShift   /)
  
  IF (.NOT. ((dir == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) .OR.     &
        &    (dir == 2 .AND. iB(1,1) == 1 .AND. iB(3,1) == 1) .OR.     &
        &    (dir == 3 .AND. iB(1,1) == 1 .AND. iB(2,1) == 1))) THEN
     dims_data   = (/0/)
     offset_mem  = (/0/)
     offset_file = (/0/)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(1,dims_file,filespace,herror)
  CALL h5screate_simple_f(1,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'SSw'              ,scalar=SSw              ) ! TEST!!! SSw, MMw unschoen ...
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'MMw'              ,scalar=MMw              )
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  IF (dir == 1) CALL write_hdf_infoREAL(dims,.FALSE.,.FALSE.,'VectorX',array=y1p(SSw:MMw))
  IF (dir == 2) CALL write_hdf_infoREAL(dims,.FALSE.,.FALSE.,'VectorY',array=y2p(SSw:MMw))
  IF (dir == 3) CALL write_hdf_infoREAL(dims,.FALSE.,.FALSE.,'VectorZ',array=y3p(SSw:MMw))
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_1D_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname
  
  INTEGER, INTENT(IN)      ::  SS1
  INTEGER, INTENT(IN)      ::  SS2
  INTEGER, INTENT(IN)      ::  SS3
  
  INTEGER, INTENT(IN)      ::  NN1
  INTEGER, INTENT(IN)      ::  NN2
  INTEGER, INTENT(IN)      ::  NN3
  
  INTEGER, INTENT(IN)      ::  vel_dir
  
  REAL   , INTENT(OUT)     ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)        ::  offset_mem(1:3)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  INTEGER                  ::  S3w, M3w, dim3
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Hinweis: phi(1:10,15,12) waere somit                                                                     !
  !       h5dump -d /dsetname -s "11,14,0" -c "1,1,10" filename.h5                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Alt, read2_hdf kann deutlich mehr!                                                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  CALL filespace_props(vel_dir,(/1,1,1/),S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  dims_data  = (/(NN1-SS1+1   ),(NN2-SS2+1   ),(NN3-SS3+1   )/)
  offset_mem = (/(SS1-b1L     ),(SS2-b2L     ),(SS3-b3L     )/)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 lesen ============================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'time'          ,scalar=time          )
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'dtime'         ,scalar=dtime         )
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'time_out_vect' ,scalar=time_out_vect )
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'time_out_scal' ,scalar=time_out_scal )
  CALL read_hdf_infoINT (1,.TRUE.,attr_yes,'timestep'      ,scalar=timestep      )
  CALL read_hdf_infoINT (1,.TRUE.,attr_yes,'write_count'   ,scalar=write_count   )
  CALL read_hdf_infoLOG (1,.TRUE.,attr_yes,'write_out_vect',scalar=write_out_vect)
  CALL read_hdf_infoLOG (1,.TRUE.,attr_yes,'write_out_scal',scalar=write_out_scal)
  CALL read_hdf_infoLOG (1,.TRUE.,attr_yes,'new_dtime'     ,scalar=new_dtime     )
  CALL read_hdf_infoINT (1,.TRUE.,attr_yes,'n_conc'        ,scalar=n_conc_old    )
  !-----------------------------------------------------------------------------------------------------------
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  CALL h5dclose_f(dset_id  ,herror)
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  ! TEST!!! neu und ungetestet:
  !===========================================================================================================
  !=== Ghost-cell update =====================================================================================
  !===========================================================================================================
  CALL exchange(1,vel_dir,phi)
  CALL exchange(2,vel_dir,phi)
  CALL exchange(3,vel_dir,phi)
  !===========================================================================================================
  
  
  END SUBROUTINE read_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read2_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname
  
  INTEGER, INTENT(IN)      ::  SS1
  INTEGER, INTENT(IN)      ::  SS2
  INTEGER, INTENT(IN)      ::  SS3
  
  INTEGER, INTENT(IN)      ::  NN1
  INTEGER, INTENT(IN)      ::  NN2
  INTEGER, INTENT(IN)      ::  NN3
  
  INTEGER, INTENT(IN)      ::  vel_dir
  
  REAL   , INTENT(OUT)     ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)        ::  offset_mem(1:3)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  INTEGER                  ::  S3w, M3w, dim3
  
  INTEGER                  ::  S1r, N1r, i0, iGrid
  INTEGER                  ::  S2r, N2r, j0, jGrid
  INTEGER                  ::  S3r, N3r, k0, kGrid
  
  INTEGER                  ::  Siw(1:3), Miw(1:3)
  
  LOGICAL                ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Hinweis: phi(1:10,15,12) waere somit                                                                     !
  !       h5dump -d /dsetname -s "11,14,0" -c "1,1,10" filename.h5                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Kann im Gegensatz zu "read_hdf" das eingelesene Feld beliebg im Raum anordnen.            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Attribute lesen =======================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time'          ,scalar=time          )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'dtime'         ,scalar=dtime         )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_vect' ,scalar=time_out_vect )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_scal' ,scalar=time_out_scal )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'timestep'      ,scalar=timestep      )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_count'   ,scalar=write_count   )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_vect',scalar=write_out_vect)
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_scal',scalar=write_out_scal)
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'new_dtime'     ,scalar=new_dtime     )
  CALL read_hdf_infoINT (3,.FALSE.,attr_yes,'S1w S2w S3w'   ,array =Siw           )
  CALL read_hdf_infoINT (3,.FALSE.,attr_yes,'M1w M2w M3w'   ,array =Miw           )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'n_conc'        ,scalar=n_conc_old    )
  !-----------------------------------------------------------------------------------------------------------
  CALL h5dclose_f(dset_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  
  !************************************************
  ! Raeumliche Verschiebung des einzulesenden Feldes:
  i0 = 0
  j0 = 0
  k0 = 0
  !************************************************
  
  S1w = Siw(1)
  S2w = Siw(2)
  S3w = Siw(3)
  
  M1w = Miw(1)
  M2w = Miw(2)
  M3w = Miw(3)
  
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  dim3 = M3w-S3w+1
  
  dims_file = (/dim1,dim2,dim3/)
  dims_mem  = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  
  !-----------------------------------------------------------------------------------------------------------
  iGrid = 0
  jGrid = 0
  kGrid = 0
  IF (vel_dir == 1 .AND. iB(1,1) > 1) THEN
     IF (BC_1L_global  == -2 .AND. ls1 ==  0) iGrid = 1
     IF (BC_1L_global > 0 .AND. ls1 ==  0) iGrid = 2
     IF (BC_1L_global > 0 .AND. ls1 == -1) iGrid = 1
  END IF
  IF (vel_dir == 2 .AND. iB(2,1) > 1) THEN
     IF (BC_2L_global  == -2 .AND. ls2 ==  0) jGrid = 1
     IF (BC_2L_global > 0 .AND. ls2 ==  0) jGrid = 2
     IF (BC_2L_global > 0 .AND. ls2 == -1) jGrid = 1
  END IF
  IF (vel_dir == 3 .AND. iB(3,1) > 1) THEN
     IF (BC_3L_global  == -2 .AND. ls3 ==  0) kGrid = 1
     IF (BC_3L_global > 0 .AND. ls3 ==  0) kGrid = 2
     IF (BC_3L_global > 0 .AND. ls3 == -1) kGrid = 1
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S1w+i0) <= (NN1+iShift) .AND. (M1w+i0) >= (SS1+iShift)) THEN
     IF ((S1w+i0) >= (SS1+iShift)) THEN
        S1r = S1w+i0-iShift
     ELSE
        S1r = SS1
     END IF
     IF ((M1w+i0) <= (NN1+iShift)) THEN
        N1r = M1w+i0-iShift
     ELSE
        N1r = NN1
     END IF
     
     dims_data  (1) = N1r-S1r+1
     offset_mem (1) = S1r-b1L
     offset_file(1) = iShift+iGrid-i0
     
     IF (offset_file(1) < 0   ) offset_file(1) = 0
     IF (offset_file(1) > dim1) offset_file(1) = dim1
  ELSE
     dims_data  (1) = 0
     offset_mem (1) = 0
     offset_file(1) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S2w+j0) <= (NN2+jShift) .AND. (M2w+j0) >= (SS2+jShift)) THEN
     IF ((S2w+j0) >= (SS2+jShift)) THEN
        S2r = S2w+j0-jShift
     ELSE
        S2r = SS2
     END IF
     IF ((M2w+j0) <= (NN2+jShift)) THEN
        N2r = M2w+j0-jShift
     ELSE
        N2r = NN2
     END IF
     
     dims_data  (2) = N2r-S2r+1
     offset_mem (2) = S2r-b2L
     offset_file(2) = jShift+jGrid-j0
     
     IF (offset_file(2) < 0   ) offset_file(2) = 0
     IF (offset_file(2) > dim2) offset_file(2) = dim2
  ELSE
     dims_data  (2) = 0
     offset_mem (2) = 0
     offset_file(2) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S3w+k0) <= (NN3+kShift) .AND. (M3w+k0) >= (SS3+kShift)) THEN
     IF ((S3w+k0) >= (SS3+kShift)) THEN
        S3r = S3w+k0-kShift
     ELSE
        S3r = SS3
     END IF
     IF ((M3w+k0) <= (NN3+kShift)) THEN
        N3r = M3w+k0-kShift
     ELSE
        N3r = NN3
     END IF
     
     dims_data  (3) = N3r-S3r+1
     offset_mem (3) = S3r-b3L
     offset_file(3) = kShift+kGrid-k0
     
     IF (offset_file(3) < 0   ) offset_file(3) = 0
     IF (offset_file(3) > dim3) offset_file(3) = dim3
  ELSE
     dims_data  (3) = 0
     offset_mem (3) = 0
     offset_file(3) = 0
  END IF
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  ! TEST!!! neu und ungetestet:
  !===========================================================================================================
  !=== Ghost-cell update =====================================================================================
  !===========================================================================================================
  CALL exchange(1,vel_dir,phi)
  CALL exchange(2,vel_dir,phi)
  CALL exchange(3,vel_dir,phi)
  !===========================================================================================================
  
  
  END SUBROUTINE read2_hdf
  
  
  
  
  
  
  
  
  
  
  ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
  SUBROUTINE read2_2D_hdf(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname
  
  INTEGER, INTENT(IN)      ::  N1
  INTEGER, INTENT(IN)      ::  N2
  
  INTEGER, INTENT(IN)      ::  SS1
  INTEGER, INTENT(IN)      ::  SS2
  
  INTEGER, INTENT(IN)      ::  NN1
  INTEGER, INTENT(IN)      ::  NN2
  
  INTEGER, INTENT(IN)      ::  iShift
  INTEGER, INTENT(IN)      ::  jShift
  
  INTEGER, INTENT(IN)      ::  dir
  
  REAL   , INTENT(OUT)     ::  phi(1:N1,1:N2)
  
  INTEGER(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
  INTEGER(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  
  INTEGER                  ::  S1r, N1r, i0
  INTEGER                  ::  S2r, N2r, j0
  
  INTEGER                  ::  Siw(1:2), Miw(1:2)
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen:                                                                                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Attribute lesen =======================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'          ,scalar=time          )
  CALL read_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime'         ,scalar=dtime         )
  CALL read_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect' ,scalar=time_out_vect )
  CALL read_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal' ,scalar=time_out_scal )
  CALL read_hdf_infoINT (1     ,.TRUE. ,attr_yes,'timestep'      ,scalar=timestep      )
  CALL read_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'   ,scalar=write_count   )
  CALL read_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect',scalar=write_out_vect)
  CALL read_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal',scalar=write_out_scal)
  CALL read_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'new_dtime'     ,scalar=new_dtime     )
  CALL read_hdf_infoINT (2     ,.FALSE.,attr_yes,'S1w S2w'       ,array =Siw           )
  CALL read_hdf_infoINT (2     ,.FALSE.,attr_yes,'M1w M2w'       ,array =Miw           )
  CALL read_hdf_infoINT (1     ,.TRUE. ,attr_yes,'n_conc'        ,scalar=n_conc_old    )
  !-----------------------------------------------------------------------------------------------------------
  CALL h5dclose_f(dset_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  
  !************************************************
  ! Raeumliche Verschiebung des einzulesenden Feldes:
  i0 = 0
  j0 = 0
  !************************************************
  
  S1w = Siw(1)
  S2w = Siw(2)
  
  M1w = Miw(1)
  M2w = Miw(2)
  
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  
  dims_file = (/dim1,dim2/)
  dims_mem  = (/N1  ,N2  /)
  
  !-----------------------------------------------------------------------------------------------------------
  IF ((S1w+i0) <= (NN1+iShift) .AND. (M1w+i0) >= (SS1+iShift)) THEN
     IF ((S1w+i0) >= (SS1+iShift)) THEN
        S1r = S1w+i0-iShift
     ELSE
        S1r = SS1
     END IF
     IF ((M1w+i0) <= (NN1+iShift)) THEN
        N1r = M1w+i0-iShift
     ELSE
        N1r = NN1
     END IF
     
     dims_data  (1) = N1r-S1r+1
     offset_mem (1) = S1r-1
     offset_file(1) = iShift-i0
     
     IF (offset_file(1) < 0   ) offset_file(1) = 0
     IF (offset_file(1) > dim1) offset_file(1) = dim1
  ELSE
     dims_data  (1) = 0
     offset_mem (1) = 0
     offset_file(1) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S2w+j0) <= (NN2+jShift) .AND. (M2w+j0) >= (SS2+jShift)) THEN
     IF ((S2w+j0) >= (SS2+jShift)) THEN
        S2r = S2w+j0-jShift
     ELSE
        S2r = SS2
     END IF
     IF ((M2w+j0) <= (NN2+jShift)) THEN
        N2r = M2w+j0-jShift
     ELSE
        N2r = NN2
     END IF
     
     dims_data  (2) = N2r-S2r+1
     offset_mem (2) = S2r-1
     offset_file(2) = jShift-j0
     
     IF (offset_file(2) < 0   ) offset_file(2) = 0
     IF (offset_file(2) > dim2) offset_file(2) = dim2
  ELSE
     dims_data  (2) = 0
     offset_mem (2) = 0
     offset_file(2) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! ACTHUNG!!! Nicht ganz klar, ob dieses Vorgehen im Sinne des Erfinders ist:
  IF (.NOT. ((dir == -1 .AND. iB(1,1) == 1  ) .OR. (dir == -2 .AND. iB(2,1) == 1  ) .OR. (dir == -3 .AND. iB(3,1) == 1   ) .OR.    &
        &    (dir ==  1 .AND. iB(1,1) == NB1) .OR. (dir ==  2 .AND. iB(2,1) == NB2) .OR. (dir ==  3 .AND. iB(3,1) == NB3))) THEN
     dims_data   = (/0,0/)
     offset_mem  = (/0,0/)
     offset_file = (/0,0/)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(2,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  ! TEST!!! fehlt noch:
  !===========================================================================================================
  !=== Ghost-cell update =====================================================================================
  !===========================================================================================================
  !CALL exchange(1,vel_dir,phi)
  !CALL exchange(2,vel_dir,phi)
  !CALL exchange(3,vel_dir,phi)
  !===========================================================================================================
  
  
  END SUBROUTINE read2_2D_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read2_1D_hdf(filename,dsetname,NN,SSp,NNp,iShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname
  
  INTEGER, INTENT(IN)      ::  NN
  INTEGER, INTENT(IN)      ::  SSp
  INTEGER, INTENT(IN)      ::  NNp
  INTEGER, INTENT(IN)      ::  iShift
  INTEGER, INTENT(IN)      ::  dir
  
  REAL   , INTENT(OUT)     ::  phi(1:NN)
  
  INTEGER(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
  INTEGER(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
  INTEGER                  ::  SSw, MMw, dims
  INTEGER                  ::  SSr, NNr, i0
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen:                                                                                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Attribute lesen =======================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoINT (1     ,.TRUE. ,attr_yes,'SSw'           ,scalar=SSw           ) ! TEST!!! SSw, MMw unschoen ...
  CALL read_hdf_infoINT (1     ,.TRUE. ,attr_yes,'MMw'           ,scalar=MMw           )
  !-----------------------------------------------------------------------------------------------------------
  CALL h5dclose_f(dset_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  
  !************************************************
  ! Raeumliche Verschiebung des einzulesenden Feldes:
  i0 = 0
  !************************************************
  
  dims = MMw-SSw+1
  
  dims_file = (/dims/)
  dims_mem  = (/NN /)
  
  !-----------------------------------------------------------------------------------------------------------
  IF ((SSw+i0) <= (NNp+iShift) .AND. (MMw+i0) >= (SSp+iShift)) THEN
     IF ((SSw+i0) >= (SSp+iShift)) THEN
        SSr = SSw+i0-iShift
     ELSE
        SSr = SSp
     END IF
     IF ((MMw+i0) <= (NNp+iShift)) THEN
        NNr = MMw+i0-iShift
     ELSE
        NNr = NNp
     END IF
     
     dims_data  (1) = NNr-SSr+1
     offset_mem (1) = SSr-1
     offset_file(1) = iShift-i0
     
     IF (offset_file(1) < 0   ) offset_file(1) = 0
     IF (offset_file(1) > dims) offset_file(1) = dims
  ELSE
     dims_data  (1) = 0
     offset_mem (1) = 0
     offset_file(1) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! TEST!!! alle Prozesse lesen!
  !IF (.NOT. ((dir == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) .OR.     &
  !      &    (dir == 2 .AND. iB(1,1) == 1 .AND. iB(3,1) == 1) .OR.     &
  !      &    (dir == 3 .AND. iB(1,1) == 1 .AND. iB(2,1) == 1))) THEN
  !   dims_data   = (/0/)
  !   offset_mem  = (/0/)
  !   offset_file = (/0/)
  !END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(1,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE read2_1D_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_part_hdf(filename,all_args_yes)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  LOGICAL     , INTENT(IN) ::  all_args_yes
  
  LOGICAL                  ::  attr_yes
  
  INTEGER(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
  INTEGER(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
  INTEGER(HID_T )          ::  group_id
  INTEGER(SIZE_T)          ::  size_hint
  
  INTEGER                  ::  proc_table       (1:NB1*NB2*NB3)
  INTEGER                  ::  proc_table_global(1:NB1*NB2*NB3)
  
  INTEGER                  ::  m, g
  INTEGER                  ::  dims    (1:n_groups)
  INTEGER                  ::  dims_all(1:n_groups)
  INTEGER                  ::  dims_prev, dims_global, dims_buf, offset
  
  REAL                     ::  bufferA(1:n_args)
  REAL   , ALLOCATABLE     ::  bufferB(:)
  
  CHARACTER(LEN=1)         ::  groupNr
  
  INTEGER                  ::  S1w, M1w
  INTEGER                  ::  S2w, M2w
  INTEGER                  ::  S3w, M3w
  
  
  !--- Tuning-Parameter (noch nicht optimiert!) ---
  size_hint = 20
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------------------------------------
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)  ! TEST!!! Umsortiert! Generell besser so? Koennte auch noch weiter nach oben ...
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  !-----------------------------------------------------------------------------------------------------------
  
  dims_prev = 0
  dims      = 0
  
  DO g = 1, n_groups
     
     !========================================================================================================
     ! Create a group:
     WRITE(groupNr,'(i1.1)') g
     CALL h5gcreate_f(file_id,'group'//groupNr,group_id,herror,size_hint)
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Dimensionen und Offsets ============================================================================
     !========================================================================================================
     !--- Nach Gruppe sortieren ---
     DO m = dims_prev+1, n_part
        IF (NINT(particles(14,m)) == g) THEN
           dims(g) = dims(g) + 1
           
           bufferA  (1:n_args                  ) = particles(1:n_args,dims_prev+dims(g))
           particles(1:n_args,dims_prev+dims(g)) = particles(1:n_args,m                )
           particles(1:n_args,m                ) = bufferA  (1:n_args                  )
        END IF
     END DO
     !--------------------------------------------------------------------------------------------------------
     proc_table = 0
     proc_table(rank+1) = dims(g)
     CALL MPI_ALLREDUCE(proc_table,proc_table_global,NB1*NB2*NB3,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
     
     dims_global = 0
     DO m = 1, NB1*NB2*NB3
        dims_global = dims_global + proc_table_global(m)
     END DO
     
     offset = 0
     DO m = 1, rank
        offset = offset + proc_table_global(m)
     END DO
     
     dims_buf = dims(g)
     IF (dims_buf == 0) dims_buf = 1
     !========================================================================================================
     
     
     IF (dims_global > 0) THEN
        
        !=====================================================================================================
        dims_file   = (/dims_global/)
        dims_mem    = (/dims_buf   /)
        dims_data   = (/dims(g)    /)
        offset_mem  = (/0          /)
        offset_file = (/offset     /)
        !=====================================================================================================
        
        
        !=====================================================================================================
        ! Create file space / memory space:
        CALL h5screate_simple_f(1,dims_file,filespace,herror)
        CALL h5screate_simple_f(1,dims_mem ,memspace ,herror)
        
        ! Select hyperslab in the file space / memory space:
        CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror) ! TEST!!! Umsortiert! Generell besser so?
        CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
        !=====================================================================================================
        
        ALLOCATE(bufferB(1:dims_buf))
        !=====================================================================================================
        !--- create the dataset ---
        CALL h5dcreate_f(group_id,'X',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
        !--- umspeichern ---
        bufferB(1:dims_buf) = particles(1,(dims_prev+1):(dims_prev+dims_buf))
        
        ! Write the dataset collectively:
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        !-----------------------------------------------------------------------------------------------------
        !--- create the dataset ---
        CALL h5dcreate_f(group_id,'Y',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
        !--- umspeichern ---
        bufferB(1:dims_buf) = particles(2,(dims_prev+1):(dims_prev+dims_buf))
        
        ! Write the dataset collectively:
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        !-----------------------------------------------------------------------------------------------------
        IF (dimens == 3) THEN
        !--- create the dataset ---
        CALL h5dcreate_f(group_id,'Z',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
        !--- umspeichern ---
        bufferB(1:dims_buf) = particles(3,(dims_prev+1):(dims_prev+dims_buf))
        
        ! Write the dataset collectively:
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- create the dataset ---
        CALL h5dcreate_f(group_id,'velX',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
        !--- umspeichern ---
        bufferB(1:dims_buf) = particles(4,(dims_prev+1):(dims_prev+dims_buf))
        
        ! Write the dataset collectively:
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        !-----------------------------------------------------------------------------------------------------
        !--- create the dataset ---
        CALL h5dcreate_f(group_id,'velY',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
        !--- umspeichern ---
        bufferB(1:dims_buf) = particles(5,(dims_prev+1):(dims_prev+dims_buf))
        
        ! Write the dataset collectively:
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        !-----------------------------------------------------------------------------------------------------
        IF (dimens == 3) THEN
        !--- create the dataset ---
        CALL h5dcreate_f(group_id,'velZ',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
        !--- umspeichern ---
        bufferB(1:dims_buf) = particles(6,(dims_prev+1):(dims_prev+dims_buf))
        
        ! Write the dataset collectively:
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- create the dataset ---
        CALL h5dcreate_f(group_id,'Number',H5T_NATIVE_DOUBLE,filespace,dset_id,herror) ! TEST!!! Evtl. gleich als Integer schreiben?
        
        !--- umspeichern ---
        bufferB(1:dims_buf) = particles(13,(dims_prev+1):(dims_prev+dims_buf))
        
        ! Write the dataset collectively:
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        !-----------------------------------------------------------------------------------------------------
        IF (all_args_yes) THEN
           !--------------------------------------------------------------------------------------------------
           !--- create the dataset ---
           CALL h5dcreate_f(group_id,'fb_force',H5T_NATIVE_DOUBLE,filespace,dset_id,herror) ! TEST!!! brauchts das wirklich??
           
           !--- umspeichern ---
           bufferB(1:dims_buf) = particles(10,(dims_prev+1):(dims_prev+dims_buf))
           
           ! Write the dataset collectively:
           CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
           CALL h5dclose_f(dset_id,herror)
           !--------------------------------------------------------------------------------------------------
           !--- create the dataset ---
           CALL h5dcreate_f(group_id,'usp',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
           
           !--- umspeichern ---
           bufferB(1:dims_buf) = particles(11,(dims_prev+1):(dims_prev+dims_buf))
           
           ! Write the dataset collectively:
           CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
           CALL h5dclose_f(dset_id,herror)
           !--------------------------------------------------------------------------------------------------
           !--- create the dataset ---
           CALL h5dcreate_f(group_id,'Stp',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
           
           !--- umspeichern ---
           bufferB(1:dims_buf) = particles(12,(dims_prev+1):(dims_prev+dims_buf))
           
           ! Write the dataset collectively:
           CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
           CALL h5dclose_f(dset_id,herror)
           !--------------------------------------------------------------------------------------------------
        END IF
        !=====================================================================================================
        DEALLOCATE(bufferB)
        
        CALL h5sclose_f(memspace ,herror)
        CALL h5sclose_f(filespace,herror)
        
     END IF
     
     CALL h5gclose_f(group_id ,herror)
     
     dims_prev = dims_prev + dims(g)
     
  END DO
  !-----------------------------------------------------------------------------------------------------------
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  CALL h5gcreate_f(file_id,'attributes',dset_id,herror,size_hint)
  
  CALL MPI_ALLREDUCE(dims,dims_all,n_groups,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1       ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1       ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoREAL(1       ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1       ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1       ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1       ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoLOG (1       ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1       ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1       ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
  CALL write_hdf_infoREAL(1       ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL write_hdf_infoREAL(1       ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL write_hdf_infoINT (1       ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL write_hdf_infoLOG (1       ,.TRUE. ,attr_yes,'concentration_yes',scalar=concentration_yes)
  CALL write_hdf_infoINT (1       ,.TRUE. ,attr_yes,'n_conc'           ,scalar=n_conc           )
  CALL write_hdf_infoREAL(3       ,.FALSE.,attr_yes,'gravity'          ,array =gravity          )
  CALL write_hdf_infoREAL(1       ,.TRUE. ,attr_yes,'Re'               ,scalar=Re               )
  CALL write_hdf_infoREAL(n_spec  ,.FALSE.,attr_yes,'Rip'              ,array =Rip              )
  CALL write_hdf_infoREAL(n_spec  ,.FALSE.,attr_yes,'usp'              ,array =usp              )
  CALL write_hdf_infoREAL(n_spec  ,.FALSE.,attr_yes,'Stp'              ,array =Stp              )
  
  CALL write_hdf_infoINT (1       ,.TRUE. ,attr_yes,'n_groups'         ,scalar=n_groups         )
  CALL write_hdf_infoINT (n_groups,.FALSE.,attr_yes,'group dims'       ,array =dims_all         )
  
  CALL h5gclose_f(dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !=== Basisvektoren =========================================================================================
  !===========================================================================================================
  S1w = 2  + ls1
  S2w = 2  + ls2
  S3w = 2  + ls3
  
  M1w = M1 + ls1
  M2w = M2 + ls2
  M3w = M3 + ls3
  
  IF (BC_1L_global /= -1) THEN
     S1w = 1
     M1w = M1
  END IF
  IF (BC_2L_global /= -1) THEN
     S2w = 1
     M2w = M2
  END IF
  IF (BC_3L_global /= -1) THEN
     S3w = 1
     M3w = M3
  END IF
  !-----------------------------------------------------------------------------------------------------------
                   CALL write_hdf_infoREAL(M1w-S1w+1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w:M1w))
                   CALL write_hdf_infoREAL(M2w-S2w+1,.FALSE.,.FALSE.,'VectorY',array=y2p(S2w:M2w))
  IF (dimens == 3) CALL write_hdf_infoREAL(M3w-S3w+1,.FALSE.,.FALSE.,'VectorZ',array=y3p(S3w:M3w))
  !===========================================================================================================
  
  
  CALL h5fclose_f(file_id ,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_part_hdf
    
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_part_serial(dirname)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  dirname
  
  CHARACTER(LEN=3)         ::  bID(1:3)
  
  INTEGER                  ::  i
  
  INTEGER                  ::  S1w, M1w
  INTEGER                  ::  S2w, M2w
  INTEGER                  ::  S3w, M3w
  
  
  IF (rank == 0) CALL SYSTEM('mkdir -p '//dirname) ! TEST!!! funktioniert das?
  
  CALL MPI_BARRIER(COMM_CART,merror)
  
  IF (n_part >= 1) THEN
     
     CALL num_to_string(3,iB(1,1),bID(1))
     CALL num_to_string(3,iB(2,1),bID(2))
     CALL num_to_string(3,iB(3,1),bID(3))
     
     ! TEST!!! ASCII vs. binaer?
     !OPEN(98,FILE=dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt',STATUS='UNKNOWN',FORM='UNFORMATTED')
     OPEN(98,FILE=dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt',STATUS='UNKNOWN')
     
     
     !========================================================================================================
     !=== Schreiben ==========================================================================================
     !========================================================================================================
     WRITE(98,'(1i12)') n_part
     
     DO i = 1, n_part
        WRITE(98,'(10E25.17)') particles(1:6,i), particles(10:14,i) ! TEST!!! ok?
     END DO
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Argumente ==========================================================================================
     !========================================================================================================
     ! TEST!!! ok? mehr Argumente?
     WRITE(98,'(10i12    )') n_part, n_groups, restart, write_count, timestep, n_conc
     WRITE(98,'(10E25.17 )') time, dtime, time_out_vect, time_out_scal, dtime_out_vect, dtime_out_scal
     WRITE(98,'(10i12    )') NB1, NB2, NB3
     WRITE(98,'(10E25.17 )') L1, L2, L3
     WRITE(98,'(100E25.17)') Re, gravity(1:3), Rip(1:n_spec), usp(1:n_spec), Stp(1:n_spec)
     WRITE(98,'(10L1     )') write_out_vect, write_out_scal, new_dtime, concentration_yes
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Basisvektoren ======================================================================================
     !========================================================================================================
     S1w = 2  + ls1
     S2w = 2  + ls2
     S3w = 2  + ls3
     
     M1w = M1 + ls1
     M2w = M2 + ls2
     M3w = M3 + ls3
     
     IF (BC_1L_global /= -1) THEN
        S1w = 1
        M1w = M1
     END IF
     IF (BC_2L_global /= -1) THEN
        S2w = 1
        M2w = M2
     END IF
     IF (BC_3L_global /= -1) THEN
        S3w = 1
        M3w = M3
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO i = S1w, M1w
        WRITE(98,'(1E25.17)') y1p(i)
     END DO
     DO i = S2w, M2w
        WRITE(98,'(1E25.17)') y2p(i)
     END DO
     IF (dimens == 3) THEN
     DO i = S3w, M3w
        WRITE(98,'(1E25.17)') y3p(i)
     END DO
     END IF
     !========================================================================================================
     
     CLOSE(98)
     
  END IF
  
  
  END SUBROUTINE write_part_serial
    
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_part_hdf(filename,all_args_yes)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  filename
  LOGICAL     , INTENT(IN) ::  all_args_yes
  
  LOGICAL                  ::  attr_yes
  
  INTEGER(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  INTEGER(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
  INTEGER(HSIZE_T )        ::  maxdims    (1)
  
  INTEGER(HID_T   )        ::  group_id
  
  INTEGER                  ::  m, g
  INTEGER                  ::  n_part_all, dims, dims_prev, n_groups_old, offset
  
  INTEGER, ALLOCATABLE     ::  dims_old(:)
  REAL   , ALLOCATABLE     ::  bufferB(:)
  
  CHARACTER(LEN=1)         ::  groupNr
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 lesen ============================================================================================
  !===========================================================================================================
  !-- setup file access property list with parallel I/O access ---
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  !--- open the file collectively ---
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  CALL h5gopen_f(file_id,'attributes',dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  !--- read the attributes ---
  attr_yes = .TRUE.
  CALL read_hdf_infoINT (1           ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL read_hdf_infoINT (1           ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL read_hdf_infoREAL(1           ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL read_hdf_infoREAL(1           ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL read_hdf_infoREAL(1           ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL read_hdf_infoREAL(1           ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL read_hdf_infoLOG (1           ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL read_hdf_infoLOG (1           ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL read_hdf_infoLOG (1           ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  CALL read_hdf_infoREAL(1           ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL read_hdf_infoREAL(1           ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL read_hdf_infoINT (1           ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL read_hdf_infoINT (1           ,.TRUE. ,attr_yes,'n_groups'         ,scalar=n_groups_old     )
  !-----------------------------------------------------------------------------------------------------------
  IF (n_groups_old < n_groups) THEN
     IF (rank == 0) WRITE(*,'(a,i3,a,i3,a)') 'Reading only n_groups_old =', n_groups_old, ' particle groups instead of n_groups ='    , n_groups     , ' ...'
  END IF
  IF (n_groups_old > n_groups) THEN
     IF (rank == 0) WRITE(*,'(a,i3,a,i3,a)') 'Reading only n_groups ='    , n_groups    , ' particle groups instead of n_groups_old =', n_groups_old , ' ...'
     n_groups_old = n_groups
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ALLOCATE(dims_old(1:n_groups_old))
  CALL read_hdf_infoINT (n_groups_old,.FALSE.,attr_yes,'group dims'       ,array =dims_old         )
  !-----------------------------------------------------------------------------------------------------------
  n_part_all = 0
  DO g = 1, n_groups_old
     n_part_all = n_part_all + dims_old(g)
  END DO
  
  IF (n_part_all < n_part_max) THEN
     IF (rank == 0) WRITE(*,'(a,i10,a,i10,a)') 'Reading only n_part_all =', n_part_all, ' particles instead of n_part_max =', n_part_max, ' ...'
  END IF
  IF (n_part_all > n_part_max) THEN
     IF (rank == 0) WRITE(*,'(a,i10,a,i10,a)') 'Reading only n_part_max =', n_part_max, ' particles instead of n_part_all =', n_part_all, ' ...'
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL h5gclose_f(dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)  ! TEST!!! Umsortiert! Generell besser so? Koennte auch noch weiter nach oben ...
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_INDEPENDENT_F,herror)
  !-----------------------------------------------------------------------------------------------------------
  
  dims_prev = 0
  particles = 0.
  
  DO g = 1, n_groups_old
     
     dims = dims_old(g)
     IF (dims_prev+dims > n_part_max) dims = n_part_max - dims_prev
     
     
     IF (dims > 0) THEN
        
        !=====================================================================================================
        !--- create a group ---
        WRITE(groupNr,'(i1.1)') g
        CALL h5gopen_f(file_id,'group'//groupNr,group_id,herror)
        !=====================================================================================================
        
        
        ALLOCATE(bufferB(1:dims_old(g)))
        !=====================================================================================================
        !--- open the dataset ---
        CALL h5dopen_f(group_id,'X',dset_id,herror)
        
        CALL h5dget_space_f(dset_id,filespace,herror)
        CALL h5sget_simple_extent_dims_f(filespace,dims_file,maxdims,herror)
        
        !--- read the dataset independently ---
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        
        !--- umspeichern ---
        particles(1,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
        !-----------------------------------------------------------------------------------------------------
        !--- open the dataset ---
        CALL h5dopen_f(group_id,'Y',dset_id,herror)
        
        !--- read the dataset independently ---
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        
        !--- umspeichern ---
        particles(2,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
        !-----------------------------------------------------------------------------------------------------
        IF (dimens == 3) THEN
        !--- open the dataset ---
        CALL h5dopen_f(group_id,'Z',dset_id,herror)
        
        !--- read the dataset independently ---
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        
        !--- umspeichern ---
        particles(3,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- open the dataset ---
        CALL h5dopen_f(group_id,'velX',dset_id,herror)
        !--- read the dataset independently ---
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        
        !--- umspeichern ---
        particles(4,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
        !-----------------------------------------------------------------------------------------------------
        !--- open the dataset ---
        CALL h5dopen_f(group_id,'velY',dset_id,herror)
        !--- read the dataset independently ---
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        
        !--- umspeichern ---
        particles(5,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
        !-----------------------------------------------------------------------------------------------------
        IF (dimens == 3) THEN
        !--- open the dataset ---
        CALL h5dopen_f(group_id,'velZ',dset_id,herror)
        !--- read the dataset independently ---
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        
        !--- umspeichern ---
        particles(6,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
        END IF
        !-----------------------------------------------------------------------------------------------------
        !--- open the dataset ---
        CALL h5dopen_f(group_id,'Number',dset_id,herror)
        
        !--- read the dataset independently ---
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
        CALL h5dclose_f(dset_id,herror)
        
        !--- umspeichern ---
        particles(13,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
        !-----------------------------------------------------------------------------------------------------
        IF (all_args_yes) THEN
           !--------------------------------------------------------------------------------------------------
           ! Open the dataset: ! TEST!!! brauchts das wirklich??
           CALL h5dopen_f(group_id,'fb_force',dset_id,herror)
           !--- read the dataset independently ---
           CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
           CALL h5dclose_f(dset_id,herror)
           
           !--- umspeichern ---
           particles(10,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
           !--------------------------------------------------------------------------------------------------
           ! Open the dataset:
           CALL h5dopen_f(group_id,'usp'    ,dset_id,herror)
           !--- read the dataset independently ---
           CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
           CALL h5dclose_f(dset_id,herror)
           
           !--- umspeichern ---
           particles(11,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
           !--------------------------------------------------------------------------------------------------
           ! Open the dataset:
           CALL h5dopen_f(group_id,'Stp'    ,dset_id,herror)
           !--- read the dataset independently ---
           CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
           CALL h5dclose_f(dset_id,herror)
           
           !--- umspeichern ---
           particles(12,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
           !--------------------------------------------------------------------------------------------------
        END IF
        !=====================================================================================================
        DEALLOCATE(bufferB)
        
        particles(14,(dims_prev+1):(dims_prev+dims)) = REAL(g)
        
        CALL h5gclose_f(group_id,herror)
        
        dims_prev = dims_prev + dims
        
     END IF
     
  END DO
  !-----------------------------------------------------------------------------------------------------------
  CALL h5pclose_f(plist_id,herror)
  CALL h5fclose_f(file_id ,herror)
  !===========================================================================================================
  DEALLOCATE(dims_old)
  
  
  END SUBROUTINE read_part_hdf

  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_part_serial(dirname)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) ::  dirname
  
  CHARACTER(LEN=3)         ::  bID(1:3)
  
  INTEGER                  ::  i
  INTEGER                  ::  ios
  
  
  CALL num_to_string(3,iB(1,1),bID(1))
  CALL num_to_string(3,iB(2,1),bID(2))
  CALL num_to_string(3,iB(3,1),bID(3))
  
  
  ! TEST!!! ASCII vs. binaer?
  !OPEN(98,FILE=TRIM(dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt'),STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',IOSTAT=ios)
  OPEN(98,FILE=TRIM(dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt'),STATUS='OLD',FORM='FORMATTED',ACTION='READ',IOSTAT=ios)
          
  IF (ios == 0) THEN
     
     !===========================================================================================================
     !=== Lesen =================================================================================================
     !===========================================================================================================
     READ(98,'(1i12)') n_part
     
     DO i = 1, n_part
        READ(98,'(10E25.17)') particles(1:6,i), particles(10:14,i) ! TEST!!! ok?
     END DO
     !===========================================================================================================
     
     
     !!===========================================================================================================
     !!=== weitere Argumente =====================================================================================
     !!===========================================================================================================
     ! TEST!!! ok? mehr Argumente?
     !WRITE(98,'(10i12    )') n_part, n_groups, restart, write_count, timestep, n_conc
     !WRITE(98,'(10E25.17 )') time, dtime, time_out_vect, time_out_scal, dtime_out_vect, dtime_out_scal
     !WRITE(98,'(10i12    )') NB1, NB2, NB3
     !WRITE(98,'(10E25.17 )') L1, L2, L3
     !WRITE(98,'(100E25.17)') Re, gravity(1:3), Rip(1:n_spec), usp(1:n_spec), Stp(1:n_spec)
     !WRITE(98,'(10L1     )') write_out_vect, write_out_scal, new_dtime, concentration_yes
     !!===========================================================================================================
     
     CLOSE(98)
     
  ELSE
     
     n_part = 0
     
  END IF
  
  
  END SUBROUTINE read_part_serial

  
  
  
  
  

  SUBROUTINE write_stats_hdf_4D(filename,dsetname, SS1,SS2,SS3,SS4, NN1,NN2,NN3,NN4, phi)

  ! ---mjohn 230412

  ! function is suited to write 4D data, specifically designed for spatial mode analysis
  ! dims 1 and 2 are spatial directions, dim1 is x1, dim2 is x2, data is COLLOCATED (count from SSi to NNi).
  ! dims 3 and 4 are integer values, typically Fourier and Hermite modes (count from SSi to NNi).

  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_lib
  USE HDF5

  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname

  INTEGER, INTENT(IN)    ::  SS1
  INTEGER, INTENT(IN)    ::  SS2
  INTEGER, INTENT(IN)    ::  SS3       ! supposed to be 1 (1:NN3)
  INTEGER, INTENT(IN)    ::  SS4       ! supposed to be 0 (0:NN4)

  INTEGER, INTENT(IN)    ::  NN1       ! supposed to be size in 1 direction (wall-normal)
  INTEGER, INTENT(IN)    ::  NN2       ! supposed to be size in 2 direction (spanwise)
  INTEGER, INTENT(IN)    ::  NN3       ! supposed to be n_frequencymodes (1:NN3)
  INTEGER, INTENT(IN)    ::  NN4       ! supposed to be n_Hermitemodes-1 (0:NN4)

  REAL   , INTENT(IN)    ::  phi(SS1:NN1, SS2:NN2, SS3:NN3, SS4:NN4)

  INTEGER                ::  Fourier(SS3:NN3), Hermite(SS4:NN4)
  INTEGER                ::  stride(4)
  INTEGER(HSIZE_T )      ::  dims_file  (4), dims_mem  (4), dims_data(4), stride_mem(4)
  INTEGER(HSSIZE_T)      ::  offset_file(4), offset_mem(4)

  LOGICAL                ::  attr_yes
  INTEGER                ::  i

  INCLUDE 'mpif.h'


  ! this function is for stride == 1 only
  stride(1:4) = 1

  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================


  !===========================================================================================================
  !=== Dimensionen und Offsets (hard-coded, da nur eine Struktur von Analyse ausgeschrieben wird) ============
  !===========================================================================================================
  dims_file  = (/ M1,        M2,        NN3-SS3+1, NN4-SS4+1 /)
  dims_mem   = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)
  dims_data  = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)

  offset_mem = (/ SS1-1,     SS2-1,     0,         0         /)
  offset_file= (/ iShift,    jShift,    0,         0         /)
  stride_mem = (/ 1,         1,         1,         1         /)


  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)

  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(4, dims_file, filespace, herror)
  CALL h5screate_simple_f(4, dims_mem , memspace , herror)

  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )

  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror,stride_mem)

  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)

  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)

  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)

  DO i = SS3, NN3
     Fourier(i) = i
  END DO
  DO i = SS4, NN4
     Hermite(i) = i
  END DO

  CALL write_hdf_infoREAL(M1,       .FALSE.,.FALSE.,'VectorX',  array=y1p(SS1:NN1)     )
  CALL write_hdf_infoREAL(M2,       .FALSE.,.FALSE.,'VectorY',  array=y2p(SS2:NN2)     )
  CALL write_hdf_infoINT (NN3-SS3+1,.FALSE.,.FALSE.,'VectorDFT',array=Fourier(SS3:NN3) )
  CALL write_hdf_infoINT (NN4-SS4+1,.FALSE.,.FALSE.,'VectorHe', array=Hermite(SS4:NN4) )
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================


  END SUBROUTINE write_stats_hdf_4D




  SUBROUTINE read_stats_hdf_4D(filename,dsetname, SS1,SS2,SS3,SS4, NN1,NN2,NN3,NN4, phi)


  ! ---mjohn 230412

  ! function is suited to read 4D data, specifically designed for spatial mode analysis
  ! dims 1 and 2 are spatial directions, dim1 is x1, dim2 is x2, data is COLLOCATED (count from SSi to NNi).
  ! dims 3 and 4 are integer values, typically Fourier and Hermite modes (count from SSi to NNi).


  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_lib
  USE HDF5

  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) ::  filename
  CHARACTER(*), INTENT(IN) ::  dsetname

  INTEGER, INTENT(IN)    ::  SS1
  INTEGER, INTENT(IN)    ::  SS2
  INTEGER, INTENT(IN)    ::  SS3       ! supposed to be 1 (1:NN3)
  INTEGER, INTENT(IN)    ::  SS4       ! supposed to be 0 (0:NN4)

  INTEGER, INTENT(IN)    ::  NN1       ! supposed to be size in 1 direction (wall-normal)
  INTEGER, INTENT(IN)    ::  NN2       ! supposed to be size in 2 direction (spanwise)
  INTEGER, INTENT(IN)    ::  NN3       ! supposed to be n_frequencymodes (1:NN3)
  INTEGER, INTENT(IN)    ::  NN4       ! supposed to be n_Hermitemodes-1 (0:NN4)

  REAL   , INTENT(OUT)    ::  phi(SS1:NN1, SS2:NN2, SS3:NN3, SS4:NN4)

  INTEGER(HSIZE_T )      ::  dims_file  (4), dims_mem  (4), dims_data(4)
  INTEGER(HSSIZE_T)      ::  offset_file(4), offset_mem(4)

  INCLUDE 'mpif.h'



  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)

  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF

  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================

  !===========================================================================================================
  !=== Dimensionen und Offsets (hard-coded, da nur eine Struktur von Analyse ausgeschrieben wird) ============
  !===========================================================================================================
  dims_file  = (/ M1,          M2,          NN3-SS3+1, NN4-SS4+1 /)
  dims_mem   = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)
  dims_data  = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)

  offset_mem = (/ SS1-1,       SS2-1,       0,         0         /)
  offset_file= (/ iShift,      jShift,      0,         0         /)


  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)

  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(4,dims_mem,memspace ,herror)

  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)

  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)

  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)

  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================

  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================


  END SUBROUTINE read_stats_hdf_4D



  
  
  
  
  SUBROUTINE write_hdf_infoREAL(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(IN)  ::  NN
  LOGICAL          , INTENT(IN)  ::  scalar_yes
  LOGICAL          , INTENT(IN)  ::  attr_yes
  CHARACTER(*)     , INTENT(IN)  ::  name
  REAL   , OPTIONAL, INTENT(IN)  ::  array(1:NN)
  REAL   , OPTIONAL, INTENT(IN)  ::  scalar
  
  INTEGER(HID_T)                 ::  memspace
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  REAL                           ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (scalar_yes) THEN
     value = scalar
  ELSE
     value = array
  END IF
  
  CALL h5screate_simple_f(1,dim_mem,memspace,herror)
  
  IF (attr_yes) THEN
     CALL h5acreate_f(dset_id,name,memtypeREAL,memspace,attr_id,herror)
     CALL h5awrite_f(attr_id,memtypeREAL,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,memspace,attr_id,herror)
     CALL H5dwrite_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror,memspace,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id ,herror)
  END IF
  
  CALL h5sclose_f(memspace,herror)
  
  
  END SUBROUTINE write_hdf_infoREAL
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_hdf_infoINT(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(IN)  ::  NN
  LOGICAL          , INTENT(IN)  ::  scalar_yes
  LOGICAL          , INTENT(IN)  ::  attr_yes
  CHARACTER(*)     , INTENT(IN)  ::  name
  INTEGER, OPTIONAL, INTENT(IN)  ::  array(1:NN)
  INTEGER, OPTIONAL, INTENT(IN)  ::  scalar
  
  INTEGER(HID_T)                 ::  memspace
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (scalar_yes) THEN
     value = scalar
  ELSE
     value = array
  END IF
  
  CALL h5screate_simple_f(1,dim_mem,memspace,herror)
  
  IF (attr_yes) THEN
     CALL h5acreate_f(dset_id,name,memtypeINT ,memspace,attr_id,herror)
     CALL h5awrite_f(attr_id,memtypeINT ,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,memspace,attr_id,herror)
     CALL H5dwrite_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,memspace,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id ,herror)
  END IF
  
  CALL h5sclose_f(memspace,herror)
  
  
  END SUBROUTINE write_hdf_infoINT
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_hdf_infoLOG(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(IN)  ::  NN
  LOGICAL          , INTENT(IN)  ::  scalar_yes
  LOGICAL          , INTENT(IN)  ::  attr_yes
  CHARACTER(*)     , INTENT(IN)  ::  name
  LOGICAL, OPTIONAL, INTENT(IN)  ::  array(1:NN)
  LOGICAL, OPTIONAL, INTENT(IN)  ::  scalar
  
  INTEGER(HID_T)                 ::  memspace
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  m
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (scalar_yes) THEN
     IF (scalar) THEN
        value = 1
     ELSE
        value = 0
     END IF
  ELSE
     DO m = 1, NN
        IF (array(m)) THEN
           value(m) = 1
        ELSE
           value(m) = 0
        END IF
     END DO
  END IF
  
  CALL h5screate_simple_f(1,dim_mem,memspace,herror)
  
  IF (attr_yes) THEN
     CALL h5acreate_f(dset_id,name,memtypeINT ,memspace,attr_id,herror)
     CALL h5awrite_f(attr_id,memtypeINT ,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,memspace,attr_id,herror)
     CALL H5dwrite_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,memspace,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id ,herror)
  END IF
  
  CALL h5sclose_f(memspace,herror)
  
  
  END SUBROUTINE write_hdf_infoLOG
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_hdf_infoREAL(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(IN)  ::  NN
  LOGICAL          , INTENT(IN)  ::  scalar_yes
  LOGICAL          , INTENT(IN)  ::  attr_yes
  CHARACTER(*)     , INTENT(IN)  ::  name
  REAL   , OPTIONAL, INTENT(OUT) ::  array(1:NN)
  REAL   , OPTIONAL, INTENT(OUT) ::  scalar
  
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  REAL                           ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (attr_yes) THEN
     CALL h5aopen_name_f(dset_id,name,attr_id,herror)
     IF (rank == 0) CALL h5aread_f(attr_id,memtypeREAL,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dopen_f(file_id,name,attr_id,herror)
     IF (rank == 0) CALL H5dread_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id,herror)
  END IF
  
  CALL MPI_BCAST(value,NN,MPI_REAL8,0,COMM_CART,merror)
  
  IF (scalar_yes) THEN
     scalar = value(1)
  ELSE
     array  = value
  END IF
  
  
  END SUBROUTINE read_hdf_infoREAL
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_hdf_infoINT(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(IN)  ::  NN
  LOGICAL          , INTENT(IN)  ::  scalar_yes
  LOGICAL          , INTENT(IN)  ::  attr_yes
  CHARACTER(*)     , INTENT(IN)  ::  name
  INTEGER, OPTIONAL, INTENT(OUT) ::  array(1:NN)
  INTEGER, OPTIONAL, INTENT(OUT) ::  scalar
  
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (attr_yes) THEN
     CALL h5aopen_name_f(dset_id,name,attr_id,herror)
     IF (rank == 0) CALL h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dopen_f(file_id,name,attr_id,herror)
     IF (rank == 0) CALL H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id,herror)
  END IF
  
  CALL MPI_BCAST(value,NN,MPI_INTEGER,0,COMM_CART,merror)
  
  IF (scalar_yes) THEN
     scalar = value(1)
  ELSE
     array  = value
  END IF
  
  
  END SUBROUTINE read_hdf_infoINT
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_hdf_infoLOG(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(IN)  ::  NN
  LOGICAL          , INTENT(IN)  ::  scalar_yes
  LOGICAL          , INTENT(IN)  ::  attr_yes
  CHARACTER(*)     , INTENT(IN)  ::  name
  LOGICAL, OPTIONAL, INTENT(OUT) ::  array(1:NN)
  LOGICAL, OPTIONAL, INTENT(OUT) ::  scalar
  
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  m
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (attr_yes) THEN
     CALL h5aopen_name_f(dset_id,name,attr_id,herror)
     IF (rank == 0) CALL h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dopen_f(file_id,name,attr_id,herror)
     IF (rank == 0) CALL H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id,herror)
  END IF
  
  CALL MPI_BCAST(value,NN,MPI_LOGICAL,0,COMM_CART,merror)
  
  IF (scalar_yes) THEN
     IF (value(1) == 1) THEN
        scalar = .TRUE.
     ELSE
        scalar = .FALSE.
     END IF
  ELSE
     DO m = 1, NN
        IF (value(m) == 1) THEN
           array(m) = .TRUE.
        ELSE
           array(m) = .FALSE.
        END IF
     END DO
  END IF
  
  
  END SUBROUTINE read_hdf_infoLOG
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE filespace_props(vel_dir,stride,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(IN)  ::  vel_dir
  INTEGER          , INTENT(IN)  ::  stride(1:3)
  INTEGER          , INTENT(OUT) ::  S1w, M1w, dim1
  INTEGER          , INTENT(OUT) ::  S2w, M2w, dim2
  INTEGER          , INTENT(OUT) ::  S3w, M3w, dim3
  
  
  !-----------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen "vel_dir" nicht global eingefuehrt, so   !
  !                dass diese Routine immer wieder von Neuem aufgerufen werden muss.                                      !
  !              - Routine dient in erster Linie der kuerzeren Schreibweise / besseren Uebersicht                         !
  !              - Default bei Intervallgrenzen sind periodische / Nachbarblock-RB                                        !
  !-----------------------------------------------------------------------------------------------------------------------!
  
  
  S1w = 2 + ls1
  S2w = 2 + ls2
  S3w = 2 + ls3
  
  M1w = (M1-1)/stride(1) + 1 + ls1
  M2w = (M2-1)/stride(2) + 1 + ls2
  M3w = (M3-1)/stride(3) + 1 + ls3
  
  offset_file = (/iShift/stride(1),jShift/stride(2),kShift/stride(3)/)
  
  
  !===========================================================================================================
  IF (BC_1L_global /= -1) THEN
     IF (vel_dir == 1) THEN
        IF (BC_1L_global == -2) THEN
           S1w = 1
           IF (iB(1,1) > 1 .AND. ls1 ==  0) offset_file(1) = offset_file(1) + 1
        ELSE
           S1w = 0
           IF (iB(1,1) > 1 .AND. ls1 ==  0) offset_file(1) = offset_file(1) + 2
           IF (iB(1,1) > 1 .AND. ls1 == -1) offset_file(1) = offset_file(1) + 1
        END IF
        IF (BC_1U_global == -2) THEN
           M1w = M1-1
        ELSE
           M1w = M1
        END IF
     ELSE
        S1w = 1
        M1w = (M1-1)/stride(1) + 1
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2L_global /= -1) THEN
     IF (vel_dir == 2) THEN
        IF (BC_2L_global == -2) THEN
           S2w = 1
           IF (iB(2,1) > 1 .AND. ls2 ==  0) offset_file(2) = offset_file(2) + 1
        ELSE
           S2w = 0
           IF (iB(2,1) > 1 .AND. ls2 ==  0) offset_file(2) = offset_file(2) + 2
           IF (iB(2,1) > 1 .AND. ls2 == -1) offset_file(2) = offset_file(2) + 1
        END IF
        IF (BC_2U_global == -2) THEN
           M2w = M2-1
        ELSE
           M2w = M2
        END IF
     ELSE
        S2w = 1
        M2w = (M2-1)/stride(2) + 1
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_3L_global /= -1) THEN
     IF (vel_dir == 3) THEN
        IF (BC_3L_global == -2) THEN
           S3w = 1
           IF (iB(3,1) > 1 .AND. ls3 ==  0) offset_file(3) = offset_file(3) + 1
        ELSE
           S3w = 0
           IF (iB(3,1) > 1 .AND. ls3 ==  0) offset_file(3) = offset_file(3) + 2
           IF (iB(3,1) > 1 .AND. ls3 == -1) offset_file(3) = offset_file(3) + 1
        END IF
        IF (BC_3U_global == -2) THEN
           M3w = M3-1
        ELSE
           M3w = M3
        END IF
     ELSE
        S3w = 1
        M3w = (M3-1)/stride(3) + 1
     END IF
  END IF
  !===========================================================================================================
  
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  dim3 = M3w-S3w+1
  
  dims_file = (/dim1,dim2,dim3/)
  
  
  END SUBROUTINE filespace_props
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE lambda
  
  IMPLICIT NONE
  
  REAL                   ::  Dvel(1:3,1:3)
  REAL                   ::  TT(1:6)
  
  REAL                   ::  PP, QQ, RR
  REAL                   ::  rho, theta
  REAL                   ::  eps
  
  REAL                   ::  temp
  REAL                   ::  lam(1:3)
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  INTEGER                ::  m, n
  REAL                   ::  pi
  
  
  !-----------------------------------------------------------------!
  ! Anmerkung: [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]    !
  !-----------------------------------------------------------------!
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  eps = 10.**(-5)
  
  !===========================================================================================================
  !=== Geschwindigkeiten interpolieren =======================================================================
  !===========================================================================================================
  
  ! TEST!!! kann man auch weglassen (wird schon vorher ausgefuehrt), bzw. exchange-Parameter einfuehren!
  CALL exchange(1,1,vel(b1L,b2L,b3L,1))
  CALL exchange(2,2,vel(b1L,b2L,b3L,2))
  CALL exchange(3,3,vel(b1L,b2L,b3L,3))
  
  !===========================================================================================================
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              nl(i,j,k,1) = nl(i,j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1)
           END DO
        END DO
     END DO
  END DO
  
  !===========================================================================================================
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              nl(i,j,k,2) = nl(i,j,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
           END DO
        END DO
     END DO
  END DO
  
  !===========================================================================================================
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,3) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              nl(i,j,k,3) = nl(i,j,k,3) + cIwp(kk,k)*vel(i,j,k+kk,3)
           END DO
        END DO
     END DO
  END DO
  
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== lambda ================================================================================================
  !===========================================================================================================
  
  CALL exchange_all_all(.FALSE.,nl)
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           
           Dvel = 0.
           
           !--- d/dx -----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO ii = b1L, b1U
              Dvel(:,1) = Dvel(:,1) + cp1(ii,i)*nl(i+ii,j,k,:)
           END DO
           
           !--- d/dy -----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO jj = b2L, b2U
              Dvel(:,2) = Dvel(:,2) + cp2(jj,j)*nl(i,j+jj,k,:)
           END DO
           
           !--- d/dz -----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO kk = b3L, b3U
              Dvel(:,3) = Dvel(:,3) + cp3(kk,k)*nl(i,j,k+kk,:)
           END DO
           
           
           !--- Tensor (symmetrisch) -------------------------------------------------------------------------
           TT(1) = 2.*(Dvel(1,1)*Dvel(1,1) + Dvel(1,2)*Dvel(2,1) + Dvel(1,3)*Dvel(3,1))
           TT(2) = 2.*(Dvel(2,1)*Dvel(1,2) + Dvel(2,2)*Dvel(2,2) + Dvel(2,3)*Dvel(3,2))
           TT(3) = 2.*(Dvel(3,1)*Dvel(1,3) + Dvel(3,2)*Dvel(2,3) + Dvel(3,3)*Dvel(3,3))
           
           TT(4) = (Dvel(1,1) + Dvel(2,2))*(Dvel(1,2) + Dvel(2,1)) + Dvel(1,3)*Dvel(3,2) + Dvel(2,3)*Dvel(3,1)
           TT(5) = (Dvel(1,1) + Dvel(3,3))*(Dvel(1,3) + Dvel(3,1)) + Dvel(1,2)*Dvel(2,3) + Dvel(3,2)*Dvel(2,1)
           TT(6) = (Dvel(2,2) + Dvel(3,3))*(Dvel(2,3) + Dvel(3,2)) + Dvel(2,1)*Dvel(1,3) + Dvel(3,1)*Dvel(1,2)
           
           
           !--- Invarianten ----------------------------------------------------------------------------------
           PP = TT(1) + TT(2) + TT(3)
           QQ = TT(1)*TT(2) + TT(1)*TT(3) + TT(2)*TT(3) - TT(4)**2 - TT(5)**2 - TT(6)**2
           RR = TT(1)*TT(2)*TT(3) + 2.*TT(4)*TT(5)*TT(6) - TT(1)*TT(6)**2 - TT(2)*TT(5)**2 - TT(3)*TT(4)**2
           
           
           !--- Eigenwerte -----------------------------------------------------------------------------------
           rho = (PP/3.)**2 - QQ/3.
           IF (rho <= eps) THEN
              !----------------------------------------------------------------------------------------------!
              ! y:=lam-PP/3.                                                                                 !
              ! y**3+p*y+q=0.                                                                                !
              ! p:=QQ-PP**2/3.                                                                               !
              ! q:=-2*(PP/3)**3+PP*QQ/3.-RR                                                                  !
              ! TT ist symmetrisch ==> lam(1:3) reell <==> D=(p/3.)**3+(q/2.)**2 .LE. 0.                     !
              !                    ==> falls p=QQ-PP**2/3.=0.                                                !
              !                    ==> q:=-2*(PP/3)**3+PP**3/9.-RR = (PP/3)**3-RR                            !
              !                    ==> D=(q/2.)**2 .LE. 0. <==> q=0. <==> RR=(PP/3)**3 <==> lam(1:3)=lam(3)  !
              !                    ==> y=0.                                                                  !
              !                    ==> lam=PP/3.=RR**(1./3.)                                                 !
              !----------------------------------------------------------------------------------------------!
              res(i,j,k) = PP/3.
           ELSE
              rho = SQRT(rho)
              QQ  = ((PP/3.)**3 - QQ*(PP/6.) + RR/2.)/rho**3
              
              IF (ABS(QQ) < 1.) THEN
                 theta = ACOS(QQ)/3.
              ELSE
                 IF (QQ > 0.) THEN
                    theta = 0.
                 ELSE
                    theta = pi/3.
                 END IF
              END IF
              
              lam(1) = COS(theta           )
              lam(2) = COS(theta + 2.*pi/3.)
              lam(3) = COS(theta + 4.*pi/3.)
              
              
              !--- sortieren ---------------------------------------------------------------------------------
              DO m = 1, 3
!pgi$ unroll = n:8
                 DO n = m+1, 3
                    IF (lam(m) > lam(n)) THEN
                       temp   = lam(n)
                       lam(n) = lam(m)
                       lam(m) = temp
                    END IF
                 END DO
              END DO
              
              ! Faktor 1/2, da bei TT(1:6) mit 2 durchmultipliziert wurde ...
              res(i,j,k) = rho*lam(2) + PP/6.
           END IF
           
        END DO
     END DO
  END DO
  
  
  END SUBROUTINE lambda
  
  
  
END MODULE mod_inout
