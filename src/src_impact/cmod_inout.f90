!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* Apr 2012                                                                                                  *
!*************************************************************************************************************

module cmod_inout
  
    use mod_dims
    use mod_vars
    use mod_exchange
    use mod_lib
    use HDF5
  
    use mpi

    use iso_c_binding

    !    private
  
    public write_fields
    public write_restart, read_restart
    public write_hdf, read_hdf, read2_hdf, write_hdf_velall
    public write_2D_hdf, read2_2D_hdf
    public write_1D_hdf, read2_1D_hdf
    public write_part_hdf, read_part_hdf
    public write_part_serial, read_part_serial
    public write_stats_hdf_4D, read_stats_hdf_4D
  
    public write_hdf_infoREAL, read_hdf_infoREAL
    public write_hdf_infoINT , read_hdf_infoINT
    public write_hdf_infoLOG , read_hdf_infoLOG
  
    public filespace_props
    public lambda
  
  
  
    integer(HID_T)                 ::  file_id, plist_id, dset_id
    integer(HID_T)                 ::  filespace, memspace
    integer(HID_T)                 ::  memtypeREAL, memtypeINT
  
    integer(HSIZE_T )              ::  dims_file  (1:3) ! wird global eingefuehrt, um Probleme bei Uebergabe als Argument nach "CALL filespace_props" zu vermeiden.
    integer(HSSIZE_T)              ::  offset_file(1:3)
  
    integer                        ::  i, j, k ! TEST!!! wo brauchts das global?
  
contains

    subroutine SF_write3D(  &
        N,                  &
        bL, bU,             &
        SS, NN,             &
        phi,                &
        filecount ) bind (c,name='SF_write3D')

        implicit none

        integer(c_int), intent(in)      ::  N(3)

        integer(c_int), intent(in)      ::  bL(3)
        integer(c_int), intent(in)      ::  bU(3)

        integer(c_int), intent(in)      ::  SS(3)

        integer(c_int), intent(in)      ::  NN(3)

        !        integer                ::  m

        integer                ::  i!, ii
        integer                ::  j!, jj
        integer                ::  k!, kk


        real(c_double), intent(out)   ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        integer(c_int)  , intent(in)    :: filecount


        character(len=5)       ::  count_char
        !  character(len=1)       ::  conc_number



        !===========================================================================================================
        !=== Ausschrieb-Nr. als String fuer File-Namen =============================================================
        !===========================================================================================================
        !  call num_to_string(5,write_count,count_char)
        call num_to_string(5,filecount,count_char)
        !===========================================================================================================



        !===========================================================================================================
        if (write_large) call write_hdf('pre_large_'//count_char,'pre',SS(1),SS(2),SS(3),NN(1),NN(2),NN(3),0,stride_large,phi)
        if (write_med  ) call write_hdf('pre_med_'  //count_char,'pre',SS(1),SS(2),SS(3),NN(1),NN(2),NN(3),0,stride_med  ,phi)
        if (write_small) call write_hdf('pre_small_'//count_char,'pre',SS(1),SS(2),SS(3),NN(1),NN(2),NN(3),0,stride_small,phi)
        !===========================================================================================================

    end subroutine SF_write3D


!      !> \todo include ls1, BC_global, IB, y1p ...
!    ! TEST!!! evtl. schoener: write_hdf_2D??
!    ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
!    subroutine SF_write_HDF5_2D(    &
!        filename,                   &
!        dsetname,                   &
!        ls,                         &
!        M,                          &
!        N,                          &
!        SS,                         &
!        NN,                         &
!        BC_L,                       &
!        BC_U,                       &
!        NB,                         &
!        iB,                         &
!        shift,                      &
!        dir,phi)
!
!        implicit none
!
!        character(*), intent(in) ::  filename
!        character(*), intent(in) ::  dsetname
!
!        integer, intent(in)      ::  ls(1:3)
!
!        integer, intent(in)      ::  M(1:3)
!        integer, intent(in)      ::  N(1:3)
!
!        integer, intent(in)      ::  SS(1:3)
!        integer, intent(in)      ::  NN(1:3)
!
!        integer, intent(in)      ::  BC_L(1:3)
!        integer, intent(in)      ::  BC_U(1:3)
!
!        integer, intent(in)      ::  NB(1:3)
!
!        integer, intent(in)      ::  iB(1:3)
!        integer, intent(in)      ::  shift(1:3)
!
!        integer, intent(in)      ::  dir
!
!        real   , intent(in)      ::  phi(1:N1,1:N2)
!
!        integer                  ::  S1w, M1w, dim1
!        integer                  ::  S2w, M2w, dim2
!
!        logical                  ::  attr_yes
!
!        integer(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
!        integer(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
!        !                eingefuehrt.                                                                              !
!        !----------------------------------------------------------------------------------------------------------!
!
!
!        !===========================================================================================================
!        !=== Attribut-Datentyp =====================================================================================
!        !===========================================================================================================
!        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
!        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Dimensionen und Offsets ===============================================================================
!        !===========================================================================================================
!        if (ABS(dir) == 1) then
!            S1w = 2  + ls(2)
!            S2w = 2  + ls(3)
!
!            M1w = M(2) + ls(2)
!            M2w = M(3) + ls(3)
!
!            if (BC_L(2) /= -1) then
!                S1w = 1
!                M1w = M(2)
!            end if
!            if (BC_L(3) /= -1) then
!                S2w = 1
!                M2w = M(3)
!            end if
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        if (ABS(dir) == 2) then
!            S1w = 2  + ls(1)
!            S2w = 2  + ls(3)
!
!            M1w = M(1) + ls(1)
!            M2w = M(3) + ls(3)
!
!            if (BC_L(1) /= -1) then
!                S1w = 1
!                M1w = M(1)
!            end if
!            if (BC_L(3) /= -1) then
!                S2w = 1
!                M2w = M(3)
!            end if
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        if (ABS(dir) == 3) then
!            S1w = 2  + ls(1)
!            S2w = 2  + ls(2)
!
!            M1w = M(1) + ls(1)
!            M2w = M(2) + ls(2)
!
!            if (BC_L(1) /= -1) then
!                S1w = 1
!                M1w = M(1)
!            end if
!            if (BC_L(2) /= -1) then
!                S2w = 1
!                M2w = M(2)
!            end if
!        end if
!        !===========================================================================================================
!        dim1 = M1w-S1w+1
!        dim2 = M2w-S2w+1
!
!        dims_file   = (/ dim1 , dim2 /)
!        dims_mem    = (/ N(1) , N(2)   /)
!        dims_data   = (/(NN(1)-SS(1)+1),(NN(2)-SS(2)+1)/)
!        offset_mem  = (/ SS(1)-1, SS(2)-1/)
!        offset_file = (/shift(i),shift(i)/)
!
!        if (.not. ((dir == -1 .and. iB(1) == 1    )   &
!            .or.   (dir == -2 .and. iB(2) == 1    )   &
!            .or.   (dir == -3 .and. iB(3) == 1    )   &
!            .or.   (dir ==  1 .and. iB(1) == NB(1))   &
!            .or.   (dir ==  2 .and. iB(2) == NB(2))   &
!            .or.   (dir ==  3 .and. iB(3) == NB(3)))) then
!            dims_data   = (/0,0/)
!            offset_mem  = (/0,0/)
!            offset_file = (/0,0/)
!        end if
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== HDF5 schreiben ========================================================================================
!        !===========================================================================================================
!        ! Setup file access property list with parallel I/O access:
!        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
!        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
!
!        ! Create the file collectively:
!        call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
!        call h5pclose_f(plist_id,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        ! Create file space / memory space:
!        call h5screate_simple_f(2,dims_file,filespace,herror)
!        call h5screate_simple_f(2,dims_mem ,memspace ,herror)
!
!        ! Create the dataset:
!        call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        ! Write the attributes:
!        attr_yes = .true.
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'restart'          ,scalar=restart          )
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
!
!        call write_hdf_infoINT (2     ,.false.,attr_yes,'S1w S2w'          ,array =(/S1w,S2w/)      )
!        call write_hdf_infoINT (2     ,.false.,attr_yes,'M1w M2w'          ,array =(/M1w,M2w/)      )
!
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time'             ,scalar=time             )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )
!
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'concentration_yes',scalar=concentration_yes)
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'n_conc'           ,scalar=n_conc           )
!        call write_hdf_infoREAL(3     ,.false.,attr_yes,'gravity'          ,array =gravity          )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'Re'               ,scalar=Re               )
!        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Sc'               ,array =Sc               )
!        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Ric'              ,array =Ric              )
!        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'usc'              ,array =usc              )
!        !-----------------------------------------------------------------------------------------------------------
!        ! Select hyperslab in the file space / memory space:
!        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
!        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
!
!        ! Create property list for collective dataset write:
!        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
!        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
!
!        ! Write the dataset collectively:
!        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
!
!        call h5pclose_f(plist_id ,herror)
!        call h5dclose_f(dset_id  ,herror)
!        call h5sclose_f(memspace ,herror)
!        call h5sclose_f(filespace,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        if (ABS(dir) == 1) then
!            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorY',array=y2p(S1w:M1w))
!            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorZ',array=y3p(S2w:M2w))
!        end if
!        if (ABS(dir) == 2) then
!            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
!            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorZ',array=y3p(S2w:M2w))
!        end if
!        if (ABS(dir) == 3) then
!            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
!            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorY',array=y2p(S2w:M2w))
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        call h5fclose_f(file_id,herror)
!    !===========================================================================================================
!
!
!    end subroutine SF_write_HDF5_2D


    !pgi$g unroll = n:8
    !!pgi$r unroll = n:8
    !!pgi$l unroll = n:8
  
    subroutine SF_write2D(  &
        N,                  &
        bL, bU,             &
        SS, NN,             &
        shift,              &
        phi,                &
        filecount ) bind (c,name='SF_write2D')

        implicit none

        integer(c_int), intent(in)      ::  N(3)

        integer(c_int), intent(in)      ::  bL(3)
        integer(c_int), intent(in)      ::  bU(3)

        integer(c_int), intent(in)      ::  SS(3)

        integer(c_int), intent(in)      ::  NN(3)

        integer(c_int), intent(in)      ::  shift(3)

        real(c_double), intent(inout)   ::  phi( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ) )

        integer(c_int)  , intent(in)    :: filecount

        real(c_double)                  ::  temp( 1:N(1), 1:N(2) )


        integer                ::  i!, ii
        integer                ::  j!, jj
        integer                ::  k!, kk

        character(len=5)       ::  count_char
        !  character(len=1)       ::  conc_number


        if (rank == 0) write(*,'(a)') 'writing pressure field ...'

        !===========================================================================================================
        !=== Ausschrieb-Nr. als String fuer File-Namen =============================================================
        !===========================================================================================================
        !  call num_to_string(5,write_count,count_char)
        call num_to_string(5,filecount,count_char)
        !===========================================================================================================

        !===========================================================================================================
        k = 1
        do j = SS(2), NN(2)
            do i = SS(1), NN(1)
                temp(i,j) = phi(i,j,k)
            end do
        end do
        call write_2D_hdf('pre_'//count_char,'pre',N(1),N(2),SS(1),SS(2),NN(1),NN(2),Shift(1),Shift(2),-3,temp(1,1))
    !===========================================================================================================

    end subroutine SF_write2D


  
!    subroutine SF_write(   &
!        N,                 &
!        bL, bU,            &
!        SS, NN,            &
!        phi,               &
!        filecount ) bind (c,name='SF_write')
!
!        implicit none
!
!        integer(c_int), intent(in)      ::  N(3)
!
!        integer(c_int), intent(in)      ::  bL(3)
!        integer(c_int), intent(in)      ::  bU(3)
!
!        integer(c_int), intent(in)      ::  SS(3)
!
!        integer(c_int), intent(in)      ::  NN(3)
!
!        !        integer                ::  m
!
!        integer                ::  i!, ii
!        integer                ::  j!, jj
!        integer                ::  k!, kk
!
!
!        real(c_double), intent(out)   ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))
!        real(c_double)                ::  temp( 1:N(1), 1:N(2) )
!
!        integer(c_int)  , intent(in)    :: filecount
!
!
!        character(len=5)       ::  count_char
!        !  character(len=1)       ::  conc_number
!
!
!        if (rank == 0) write(*,'(a)') 'writing pressure field ...'
!
!        !===========================================================================================================
!        !=== Ausschrieb-Nr. als String fuer File-Namen =============================================================
!        !===========================================================================================================
!        !  call num_to_string(5,write_count,count_char)
!        call num_to_string(5,filecount,count_char)
!        !===========================================================================================================
!
!
!
!        !===========================================================================================================
!        if( dimens == 3 )then
!            if (write_large) call write_hdf('pre_large_'//count_char,'pre',SS(1),SS(2),SS(3),NN(1),NN(2),NN(3),0,stride_large,phi)
!            if (write_med  ) call write_hdf('pre_med_'  //count_char,'pre',SS(1),SS(2),SS(3),NN(1),NN(2),NN(3),0,stride_med  ,phi)
!            if (write_small) call write_hdf('pre_small_'//count_char,'pre',SS(1),SS(2),SS(3),NN(1),NN(2),NN(3),0,stride_small,phi)
!        else
!            k = 1
!            do j = SS(2), NN(2)
!                do i = SS(1), NN(1)
!                    temp(i,j) = phi(i,j,k)
!                end do
!            end do
!            call write_2D_hdf('pre_'//count_char,'pre',N(1),N(2),SS(1),SS(2),NN(1),NN(2),iShift,jShift,-3,temp(1,1))
!        end if
!    !===========================================================================================================
!
!    end subroutine SF_write
!
  
  
    subroutine VF_write( phiU, phiV, phiW, filecount ) bind (c,name='VF_write')

        implicit none

        !        integer                ::  m

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk

        real(c_double)  , intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double)  , intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double)  , intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        integer(c_int)  , intent(in)   :: filecount

        character(len=5)       ::  count_char
        !  character(len=1)       ::  conc_number

        real                   ::  tempV (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
        real                   ::  temp(1:N1,1:N2,1:2)



        if (rank == 0) write(*,'(a)') 'writing velocity field ...'


        !===========================================================================================================
        !=== Ausschrieb-Nr. als String fuer File-Namen =============================================================
        !===========================================================================================================
        !  call num_to_string(5,write_count,count_char)
        call num_to_string(5,filecount,count_char)
        !===========================================================================================================

        !===========================================================================================================
        !=== Interpolieren / Schreiben =============================================================================
        !===========================================================================================================
        if( dimens == 3 ) then
            do k = S3p, N3p
                do j = S2p, N2p
                    do i = S1p, N1p
                        tempV(i,j,k,1) = cIup(d1L,i)*phiU(i+d1L,j,k)
                        !pgi$ unroll = n:8
                        do ii = d1L+1, d1U
                            tempV(i,j,k,1) = tempV(i,j,k,1) + cIup(ii,i)*phiU(i+ii,j,k) ! TEST!!! Verallgemeinerte Interpolation verwenden? Mit compute_stats teilen?
                        end do
                    end do
                end do
            end do
            if (write_large) call write_hdf('velX_large_'//count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,tempV(b1L,b2L,b3L,1)) ! TEST!!! velX --> vel1, etc. umbenennen!!
            if (write_med  ) call write_hdf('velX_med_'  //count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,tempV(b1L,b2L,b3L,1))
            if (write_small) call write_hdf('velX_small_'//count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,tempV(b1L,b2L,b3L,1))
        else
            k = 1
            do j = S2p, N2p
                do i = S1p, N1p
                    temp(i,j,1) = cIup(d1L,i)*phiU(i+d1L,j,k)
                    !pgi$ unroll = n:8
                    do ii = d1L+1, d1U
                        temp(i,j,1) = temp(i,j,1) + cIup(ii,i)*phiU(i+ii,j,k)
                    end do
                end do
            end do
            call write_2D_hdf('velX_'//count_char,'velX',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,temp(1,1,1))
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (dimens == 3 ) then
            do k = S3p, N3p
                do j = S2p, N2p
                    do i = S1p, N1p
                        tempV(i,j,k,2) = cIvp(d2L,j)*phiV(i,j+d2L,k)
                        !pgi$ unroll = n:8
                        do jj = d2L+1, d2U
                            tempV(i,j,k,2) = tempV(i,j,k,2) + cIvp(jj,j)*phiV(i,j+jj,k)
                        end do
                    end do
                end do
            end do
            if (write_large) call write_hdf('velY_large_'//count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,tempV(b1L,b2L,b3L,2))
            if (write_med  ) call write_hdf('velY_med_'  //count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,tempV(b1L,b2L,b3L,2))
            if (write_small) call write_hdf('velY_small_'//count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,tempV(b1L,b2L,b3L,2))
        else
            k = 1
            do j = S2p, N2p
                do i = S1p, N1p
                    temp(i,j,2) = cIvp(d2L,j)*phiV(i,j+d2L,k)
                    !pgi$ unroll = n:8
                    do jj = d2L+1, d2U
                        temp(i,j,2) = temp(i,j,2) + cIvp(jj,j)*phiV(i,j+jj,k)
                    end do
                end do
            end do
            call write_2D_hdf('velY_'//count_char,'velY',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,temp(1,1,2))
        end if
        !-----------------------------------------------------------------------------------------------------------
        if( dimens == 3 ) then
            do k = S3p, N3p
                do j = S2p, N2p
                    do i = S1p, N1p
                        tempV(i,j,k,3) = cIwp(d3L,k)*phiW(i,j,k+d3L)
                        !pgi$ unroll = n:8
                        do kk = d3L+1, d3U
                            tempV(i,j,k,3) = tempV(i,j,k,3) + cIwp(kk,k)*phiW(i,j,k+kk)
                        end do
                    end do
                end do
            end do
            if (write_large) call write_hdf('velZ_large_'//count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,tempV(b1L,b2L,b3L,3))
            if (write_med  ) call write_hdf('velZ_med_'  //count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,tempV(b1L,b2L,b3L,3))
            if (write_small) call write_hdf('velZ_small_'//count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,tempV(b1L,b2L,b3L,3))
        end if
        !===========================================================================================================
        if( 1 == 2 .and. dimens == 3) call write_hdf_velall('velA_'//count_char,S1p,S2p,S3p,N1p,N2p,N3p,0,tempV)
    !===========================================================================================================
    !===========================================================================================================

    end subroutine VF_write
  
  
  
  
  
  
  
    !> \todo implement
    subroutine write_restart
  
        implicit none
  
        integer                ::  m
        character(len=3)       ::  next_restart_char
        character(len=1)       ::  conc_number
  
  
        if (write_restart_yes) then
     
            if (rank == 0) write(*,'(a,i3,a)') 'writing data for restart',restart,' ...'
     
            !========================================================================================================
            !=== neue Restart-Nr. als String fuer File-Namen ========================================================
            !========================================================================================================
            call num_to_string(3,restart,next_restart_char)
            !========================================================================================================
     
     
            !========================================================================================================
            !=== Schreiben ==========================================================================================
            !========================================================================================================
            call write_hdf('velX_restart'//next_restart_char,'velX_restart',S11B,S21B,S31B,N11B,N21B,N31B,1,(/1,1,1/),vel(b1L,b2L,b3L,1))
            call write_hdf('velY_restart'//next_restart_char,'velY_restart',S12B,S22B,S32B,N12B,N22B,N32B,2,(/1,1,1/),vel(b1L,b2L,b3L,2))
            if (dimens == 3) call write_hdf('velZ_restart'//next_restart_char,'velZ_restart',S13B,S23B,S33B,N13B,N23B,N33B,3,(/1,1,1/),vel(b1L,b2L,b3L,3))
            !--------------------------------------------------------------------------------------------------------
            call write_hdf('pre_restart'//next_restart_char,'pre_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),pre)
            !--------------------------------------------------------------------------------------------------------
            if (concentration_yes) then
                do m = 1, n_conc
                    write(conc_number,'(i1.1)') m
                    call write_hdf('conc'//conc_number//'_restart'//next_restart_char,'conc'//conc_number//'_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),conc(b1L,b2L,b3L,m))
           
                    ! TEST!!! noch nicht vollstaendig! (betrifft auch andere Routinen)
                    if (dimens == 3) then ! TEST!!! "sed" und "depo" koennen in 2D nicht geschrieben werden!!
                        call write_2D_hdf('depo'//conc_number//'_restart'//next_restart_char,'depo'//conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_conc(1,1,m))
                        call write_2D_hdf('sed' //conc_number//'_restart'//next_restart_char,'sed' //conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,conc1L    (1,1,m))
                    end if
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            if (particles_yes) then
                !CALL write_part_hdf('part_restart'//next_restart_char,.TRUE.) ! TEST!!!
                call write_part_serial('part_restart'//next_restart_char)
        
                if (dimens == 3) then
                    ! TEST!!! Schmutziger Umweg:
                    res(S1p,1:N2,1:N3) = dep1L_part(1:N2,1:N3)
                    call exchange_part(2,0,1,res)
                    call exchange_part(3,0,1,res)
                    dep1L_part(1:N2,1:N3) = res(S1p,1:N2,1:N3)
           
                    if (BC_2L == 0 .and. ls2 ==  0) dep1L_part(1 ,1:N3) = 0.
                    if (BC_2U == 0 .and. ls2 == -1) dep1L_part(N2,1:N3) = 0.
                    if (BC_3L == 0 .and. ls3 ==  0) dep1L_part(1:N2,1 ) = 0.
                    if (BC_3U == 0 .and. ls3 == -1) dep1L_part(1:N2,N3) = 0.
           
                    call write_2D_hdf('depo_part_restart'//next_restart_char,'depo_part_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_part(1,1)) ! TEST!!! ok?
                end if
            end if
           !========================================================================================================

        end if
  
  
    end subroutine write_restart
  
  
  
  
  
  
  
  
  
  
    !> \todo implement
    subroutine read_restart
  
        implicit none
  
        integer                ::  m
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk

        real                   ::  old_dtime_out_scal
        real                   ::  old_dtime_out_vect

        character(len=1)       ::  conc_number
  
  
        if (rank == 0) write(*,'(a,i3,a)') 'reading data for restart',restart,' ...'
        if (rank == 0) write(*,*)
  
        !-----------------------------------------------------------------------------------------------------------
        ! Anmerkung: Nur fuer den Notfall ...
        if (1 == 2) then
  
            conc = 0.
            vel  = 0.
  
            call read2_hdf('conc1_large_0174','conc1',S1p,S2p,S3p,N1p,N2p,N3p,0,conc(b1L,b2L,b3L,1))
            call read2_hdf('velX_large_0174','velX',S1p,S2p,S3p,N1p,N2p,N3p,0,nl(b1L,b2L,b3L,1))
            call read2_hdf('velY_large_0174','velY',S1p,S2p,S3p,N1p,N2p,N3p,0,nl(b1L,b2L,b3L,2))
            call read2_hdf('velZ_large_0174','velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,nl(b1L,b2L,b3L,3))
  
            call exchange(1,1,nl(b1L,b2L,b3L,1))
            call exchange(2,2,nl(b1L,b2L,b3L,2))
            call exchange(3,3,nl(b1L,b2L,b3L,3))
  
            if (1 == 1) then
                do k = S31, N31
                    do j = S21, N21
                        do i = S11, N11
                            vel(i,j,k,1) = cIpu(g1L,i)*nl(i+g1L,j,k,1)
                            !pgi$ unroll = n:8
                            do ii = g1L+1, g1U
                                vel(i,j,k,1) = vel(i,j,k,1) + cIpu(ii,i)*nl(i+ii,j,k,1)
                            end do
                        end do
                    end do
                end do
                do k = S32, N32
                    do j = S22, N22
                        do i = S12, N12
                            vel(i,j,k,2) = cIpv(g2L,j)*nl(i,j+g2L,k,2)
                            !pgi$ unroll = n:8
                            do jj = g2L+1, g2U
                                vel(i,j,k,2) = vel(i,j,k,2) + cIpv(jj,j)*nl(i,j+jj,k,2)
                            end do
                        end do
                    end do
                end do
                do k = S33, N33
                    do j = S23, N23
                        do i = S13, N13
                            vel(i,j,k,3) = cIpw(g3L,k)*nl(i,j,k+g3L,3)
                            !pgi$ unroll = n:8
                            do kk = g3L+1, g3U
                                vel(i,j,k,3) = vel(i,j,k,3) + cIpw(kk,k)*nl(i,j,k+kk,3)
                            end do
                        end do
                    end do
                end do
            end if
  
        else
  
            call read2_hdf('velX_restart'//restart_char,'velX_restart',S11B,S21B,S31B,N11B,N21B,N31B,1,vel(b1L,b2L,b3L,1))
            call read2_hdf('velY_restart'//restart_char,'velY_restart',S12B,S22B,S32B,N12B,N22B,N32B,2,vel(b1L,b2L,b3L,2))
            if (dimens == 3) call read2_hdf('velZ_restart'//restart_char,'velZ_restart',S13B,S23B,S33B,N13B,N23B,N33B,3,vel(b1L,b2L,b3L,3))
            !-----------------------------------------------------------------------------------------------------------
            call read2_hdf('pre_restart'//restart_char,'pre_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,pre)
            !-----------------------------------------------------------------------------------------------------------
            if (concentration_yes) then
                do m = 1, n_conc_old
                    write(conc_number,'(i1.1)') m
                    call read2_hdf('conc'//conc_number//'_restart'//restart_char,'conc'//conc_number//'_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,conc(b1L,b2L,b3L,m))
        
                    if (dimens == 3) then ! TEST!!! "sed" und "depo" koennen in 2D nicht geschrieben werden!!
                        call read2_2D_hdf('depo'//conc_number//'_restart'//restart_char,'depo'//conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_conc(1,1,m))
                        call read2_2D_hdf('sed' //conc_number//'_restart'//restart_char,'sed' //conc_number//'_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,conc1L    (1,1,m))
                    end if
                end do
            end if
            !-----------------------------------------------------------------------------------------------------------
            if (particles_yes) then
                !CALL read_part_hdf('part_restart'//restart_char,.TRUE.) ! TEST!!!
                call read_part_serial('part_restart'//restart_char)
                if (dimens == 3) then
                    dep1L_part = 0.
                    call read2_2D_hdf('depo_part_restart'//restart_char,'depo_part_restart',N2,N3,S2p,S3p,N2p,N3p,jShift,kShift,-1,dep1L_part(1,1)) ! TEST!!! ok?
                end if
            end if
            !-----------------------------------------------------------------------------------------------------------

            ! --- mjohn 070213 - allow for possibility of adaptation of dtime_out_vect or dtime_out_scal by config.txt after restart
            ! compute next time_out_vect and next time_out_scal by subtracting last dtime (stored in restart file) and adding new dtime (read from config.txt)

            ! Setup file access property list with parallel I/O access
            call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
            call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

            ! Open file collectively
            call h5fopen_f('velX_restart'//restart_char//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
            call h5pclose_f(plist_id,herror)

            call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
            call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
            !-----------------------------------------------------------------------------------------------------------
            ! Open the dataset:
            call h5dopen_f(file_id,'velX_restart',dset_id,herror)
            call read_hdf_infoREAL(1,.true. ,.true.,'dtime_out_scal' ,scalar=old_dtime_out_scal)
            call read_hdf_infoREAL(1,.true. ,.true.,'dtime_out_vect'  ,scalar=old_dtime_out_vect)
            call h5dclose_f(dset_id,herror)
            !===========================================================================================================

            time_out_vect = (time_out_vect-old_dtime_out_vect) + dtime_out_vect
            time_out_scal = (time_out_scal-old_dtime_out_scal) + dtime_out_scal

            ! catch values less than or equal to 'time'
            if(time_out_vect <= time) then
                time_out_vect = time + dtime
            end if
            if(time_out_scal <= time) then
                time_out_scal = time + dtime
            end if


        end if
  
  
    end subroutine read_restart
  
  
  
  
  
  
  
  
  
  
  
    !> \todo
    subroutine write_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,stride,phi)
  
        implicit none
  
        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        integer, intent(in)    ::  vel_dir
        integer, intent(in)    ::  stride(1:3)
  
        real   , intent(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer(HSIZE_T )      ::  stride_mem(1:3)
        integer(HSIZE_T )      ::  dims_mem  (1:3), dims_data(1:3)
        integer(HSSIZE_T)      ::  offset_mem(1:3)
  
        integer                ::  S1w, M1w, dim1
        integer                ::  S2w, M2w, dim2
        integer                ::  S3w, M3w, dim3
  
        logical                ::  attr_yes
  
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
        !                eingefuehrt.                                                                              !
        !----------------------------------------------------------------------------------------------------------!
  
  
        ! Nur fuer Druckgitter! (sonst wird es richtig kompliziert, siehe Multigrid fuer Helmholtz-Problem!).
        if (vel_dir /= 0) then
            if (stride(vel_dir) /= 1) then
                if (rank == 0) write(*,*) 'ERROR! Cannot write velocity field with stride /= 1 in the corresponding direction!'
                call MPI_FINALIZE(merror)
                stop
            end if
        end if
        if ((MOD((N1-1),stride(1)) /= 0) .or. (MOD((N2-1),stride(2)) /= 0) .or. (MOD((N3-1),stride(3)) /= 0)) then
            if (rank == 0) write(*,*) 'ERROR! Cannot write field with this stride!'
            call MPI_FINALIZE(merror)
            stop
        end if
  
  
        !===========================================================================================================
        !=== Attribut-Datentyp =====================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Dimensionen und Offsets ===============================================================================
        !===========================================================================================================
        call filespace_props(vel_dir,stride,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
        dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
        dims_data  = (/((NN1-SS1)/stride(1)+1),((NN2-SS2)/stride(2)+1),((NN3-SS3)/stride(3)+1)/)
        offset_mem = (/(SS1-b1L),(SS2-b2L),(SS3-b3L)/)
        stride_mem = (/stride(1),stride(2),stride(3)/)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== HDF5 schreiben ========================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access:
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Create the file collectively:
        call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
        call h5pclose_f(plist_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Create file space / memory space:
        call h5screate_simple_f(3,dims_file,filespace,herror)
        call h5screate_simple_f(3,dims_mem ,memspace ,herror)
  
        ! Create the dataset:
        call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Write the attributes:
        attr_yes = .true.
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'restart'          ,scalar=restart          )
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
        call write_hdf_infoINT (3     ,.false.,attr_yes,'M1  M2  M3 '      ,array =(/M1 ,M2 ,M3 /)  )
        call write_hdf_infoINT (3     ,.false.,attr_yes,'S1w S2w S3w'      ,array =(/S1w,S2w,S3w/)  )
        call write_hdf_infoINT (3     ,.false.,attr_yes,'M1w M2w M3w'      ,array =(/M1w,M2w,M3w/)  )
        call write_hdf_infoINT (3     ,.false.,attr_yes,'NB1 NB2 NB3'      ,array =(/NB1,NB2,NB3/)  )
        call write_hdf_infoINT (3     ,.false.,attr_yes,'ls1 ls2 ls3'      ,array =(/ls1,ls2,ls3/)  )
  
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time'             ,scalar=time             )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )
  
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'concentration_yes',scalar=concentration_yes)
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'n_conc'           ,scalar=n_conc           )
        call write_hdf_infoREAL(3     ,.false.,attr_yes,'gravity'          ,array =gravity          )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'Re'               ,scalar=Re               )
        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Sc'               ,array =Sc               )
        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Ric'              ,array =Ric              )
        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'usc'              ,array =usc              )
  
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'CFL   '           ,scalar=CFL              )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'thetaL'           ,scalar=thetaL           )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'epsU  '           ,scalar=epsU             )
  
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'Euler_yes'        ,scalar=Euler_yes        )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'twostep_yes'      ,scalar=twostep_yes      )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'mapping_yes'      ,scalar=mapping_yes      )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'upwind_yes'       ,scalar=upwind_yes       )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'upwind_conc_yes'  ,scalar=upwind_conc_yes  )
  
        call write_hdf_infoINT (3,.false.,attr_yes,'BC_iL',array=(/BC_1L_global,BC_2L_global,BC_3L_global/))
        call write_hdf_infoINT (3,.false.,attr_yes,'BC_iU',array=(/BC_1U_global,BC_2U_global,BC_3U_global/))
  
        call write_hdf_infoINT (dim_ncb1c,.false.,attr_yes,'ncb1c',array=ncb1c    )
        call write_hdf_infoINT (dim_ncb1g,.false.,attr_yes,'ncb1g',array=ncb1g    )
        call write_hdf_infoINT (dim_ncb1d,.false.,attr_yes,'ncb1d',array=ncb1d    )
  
        call write_hdf_infoINT (dim_ncb2c,.false.,attr_yes,'ncb2c',array=ncb2c    )
        call write_hdf_infoINT (dim_ncb2g,.false.,attr_yes,'ncb2g',array=ncb2g    )
        call write_hdf_infoINT (dim_ncb2d,.false.,attr_yes,'ncb2d',array=ncb2d    )
  
        call write_hdf_infoINT (dim_ncb3c,.false.,attr_yes,'ncb3c',array=ncb3c    )
        call write_hdf_infoINT (dim_ncb3g,.false.,attr_yes,'ncb3g',array=ncb3g    )
        call write_hdf_infoINT (dim_ncb3d,.false.,attr_yes,'ncb3d',array=ncb3d    )
  
        call write_hdf_infoREAL( M1      ,.false.,attr_yes,'y1p'  ,array=y1p(1:M1))
        call write_hdf_infoREAL( M2      ,.false.,attr_yes,'y2p'  ,array=y2p(1:M2))
        call write_hdf_infoREAL( M3      ,.false.,attr_yes,'y3p'  ,array=y3p(1:M3))
  
        call write_hdf_infoREAL((M1+1)   ,.false.,attr_yes,'y1u'  ,array=y1u(0:M1))
        call write_hdf_infoREAL((M2+1)   ,.false.,attr_yes,'y2v'  ,array=y2v(0:M2))
        call write_hdf_infoREAL((M3+1)   ,.false.,attr_yes,'y3w'  ,array=y3w(0:M3))
        !-----------------------------------------------------------------------------------------------------------
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror,stride_mem)
  
        ! Create property list for collective dataset write:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
        ! Write the dataset collectively:
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
        call h5pclose_f(plist_id ,herror)
        call h5dclose_f(dset_id  ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
        !-----------------------------------------------------------------------------------------------------------
        if (vel_dir == 1) then
            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1u(S1w::stride(1)))
        else
            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1p(S1w::stride(1)))
        end if
        if (vel_dir == 2) then
            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorY',array=y2v(S2w::stride(2)))
        else
            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorY',array=y2p(S2w::stride(2)))
        end if
        if (vel_dir == 3) then
            call write_hdf_infoREAL(dim3,.false.,.false.,'VectorZ',array=y3w(S3w::stride(3)))
        else
            call write_hdf_infoREAL(dim3,.false.,.false.,'VectorZ',array=y3p(S3w::stride(3)))
        end if
        !-----------------------------------------------------------------------------------------------------------
        call h5fclose_f(file_id,herror)
    !===========================================================================================================
  
  
    end subroutine write_hdf
  
  
  
  
  
  
  
  
  
  
  
!    subroutine write_hdf_velall(filename,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
!
!        implicit none
!
!        character(*), intent(in) ::  filename
!
!        integer, intent(in)      ::  SS1
!        integer, intent(in)      ::  SS2
!        integer, intent(in)      ::  SS3
!
!        integer, intent(in)      ::  NN1
!        integer, intent(in)      ::  NN2
!        integer, intent(in)      ::  NN3
!
!        integer, intent(in)      ::  vel_dir
!
!        real   , intent(in)      ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
!
!        integer(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
!        integer(HSSIZE_T)        ::  offset_mem(1:3)
!
!        integer                  ::  S1w, M1w, dim1
!        integer                  ::  S2w, M2w, dim2
!        integer                  ::  S3w, M3w, dim3
!
!        logical                  ::  attr_yes
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
!        !                eingefuehrt.                                                                              !
!        !----------------------------------------------------------------------------------------------------------!
!
!
!        !===========================================================================================================
!        !=== Attribut-Datentyp =====================================================================================
!        !===========================================================================================================
!        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
!        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Dimensionen und Offsets ===============================================================================
!        !===========================================================================================================
!        call filespace_props(vel_dir,(/1,1,1/),S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
!
!        dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
!        dims_data  = (/(NN1-SS1+1   ),(NN2-SS2+1   ),(NN3-SS3+1   )/)
!        offset_mem = (/(SS1-b1L     ),(SS2-b2L     ),(SS3-b3L     )/)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== HDF5 schreiben ========================================================================================
!        !===========================================================================================================
!        ! Setup file access property list with parallel I/O access:
!        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
!        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
!
!        ! Create the file collectively:
!        call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
!        call h5pclose_f(plist_id,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        !-----------------------------------------------------------------------------------------------------------
!        ! Create file space / memory space:
!        call h5screate_simple_f(3,dims_file,filespace,herror)
!        call h5screate_simple_f(3,dims_mem ,memspace ,herror)
!
!        ! Create the dataset:
!        call h5dcreate_f(file_id,'velX',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        ! Write the attributes:
!        attr_yes = .true.
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'restart'          ,scalar=restart          )
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
!
!        call write_hdf_infoINT (3     ,.false.,attr_yes,'M1  M2  M3 '      ,array =(/M1 ,M2 ,M3 /)  )
!        call write_hdf_infoINT (3     ,.false.,attr_yes,'S1w S2w S3w'      ,array =(/S1w,S2w,S3w/)  )
!        call write_hdf_infoINT (3     ,.false.,attr_yes,'M1w M2w M3w'      ,array =(/M1w,M2w,M3w/)  )
!        call write_hdf_infoINT (3     ,.false.,attr_yes,'NB1 NB2 NB3'      ,array =(/NB1,NB2,NB3/)  )
!        call write_hdf_infoINT (3     ,.false.,attr_yes,'ls1 ls2 ls3'      ,array =(/ls1,ls2,ls3/)  )
!
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time'             ,scalar=time             )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )
!
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'concentration_yes',scalar=concentration_yes)
!        call write_hdf_infoINT (1     ,.true. ,attr_yes,'n_conc'           ,scalar=n_conc           )
!        call write_hdf_infoREAL(3     ,.false.,attr_yes,'gravity'          ,array =gravity          )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'Re'               ,scalar=Re               )
!        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Sc'               ,array =Sc               )
!        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Ric'              ,array =Ric              ) ! TEST!!! Rip und usp auch beruecksichtigen!
!        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'usc'              ,array =usc              )
!
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'CFL   '           ,scalar=CFL              )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'thetaL'           ,scalar=thetaL           )
!        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'epsU  '           ,scalar=epsU             )
!
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'Euler_yes'        ,scalar=Euler_yes        )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'twostep_yes'      ,scalar=twostep_yes      )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'mapping_yes'      ,scalar=mapping_yes      )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'upwind_yes'       ,scalar=upwind_yes       )
!        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'upwind_conc_yes'  ,scalar=upwind_conc_yes  )
!
!        call write_hdf_infoINT (3,.false.,attr_yes,'BC_iL',array=(/BC_1L_global,BC_2L_global,BC_3L_global/))
!        call write_hdf_infoINT (3,.false.,attr_yes,'BC_iU',array=(/BC_1U_global,BC_2U_global,BC_3U_global/))
!
!        call write_hdf_infoINT (dim_ncb1c,.false.,attr_yes,'ncb1c',array=ncb1c    )
!        call write_hdf_infoINT (dim_ncb1g,.false.,attr_yes,'ncb1g',array=ncb1g    )
!        call write_hdf_infoINT (dim_ncb1d,.false.,attr_yes,'ncb1d',array=ncb1d    )
!
!        call write_hdf_infoINT (dim_ncb2c,.false.,attr_yes,'ncb2c',array=ncb2c    )
!        call write_hdf_infoINT (dim_ncb2g,.false.,attr_yes,'ncb2g',array=ncb2g    )
!        call write_hdf_infoINT (dim_ncb2d,.false.,attr_yes,'ncb2d',array=ncb2d    )
!
!        call write_hdf_infoINT (dim_ncb3c,.false.,attr_yes,'ncb3c',array=ncb3c    )
!        call write_hdf_infoINT (dim_ncb3g,.false.,attr_yes,'ncb3g',array=ncb3g    )
!        call write_hdf_infoINT (dim_ncb3d,.false.,attr_yes,'ncb3d',array=ncb3d    )
!
!        call write_hdf_infoREAL( M1      ,.false.,attr_yes,'y1p'  ,array=y1p(1:M1))
!        call write_hdf_infoREAL( M2      ,.false.,attr_yes,'y2p'  ,array=y2p(1:M2))
!        call write_hdf_infoREAL( M3      ,.false.,attr_yes,'y3p'  ,array=y3p(1:M3))
!
!        call write_hdf_infoREAL((M1+1)   ,.false.,attr_yes,'y1u'  ,array=y1u(0:M1))
!        call write_hdf_infoREAL((M2+1)   ,.false.,attr_yes,'y2v'  ,array=y2v(0:M2))
!        call write_hdf_infoREAL((M3+1)   ,.false.,attr_yes,'y3w'  ,array=y3w(0:M3))
!        !-----------------------------------------------------------------------------------------------------------
!        ! Select hyperslab in the file space / memory space:
!        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
!        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
!
!        ! Create property list for collective dataset write:
!        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
!        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
!
!        ! Write the dataset collectively:
!        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,1),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
!
!        call h5pclose_f(plist_id ,herror)
!        call h5dclose_f(dset_id  ,herror)
!        call h5sclose_f(memspace ,herror)
!        call h5sclose_f(filespace,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        !-----------------------------------------------------------------------------------------------------------
!        ! Create file space / memory space:
!        call h5screate_simple_f(3,dims_file,filespace,herror)
!        call h5screate_simple_f(3,dims_mem ,memspace ,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        ! Create the dataset:
!        call h5dcreate_f(file_id,'velY',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        ! Select hyperslab in the file space / memory space:
!        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
!        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
!
!        ! Create property list for collective dataset write:
!        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
!        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
!
!        ! Write the dataset collectively:
!        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,2),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
!
!        call h5pclose_f(plist_id ,herror)
!        call h5dclose_f(dset_id  ,herror)
!        call h5sclose_f(memspace ,herror)
!        call h5sclose_f(filespace,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        !-----------------------------------------------------------------------------------------------------------
!        ! Create file space / memory space:
!        call h5screate_simple_f(3,dims_file,filespace,herror)
!        call h5screate_simple_f(3,dims_mem ,memspace ,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        ! Create the dataset:
!        call h5dcreate_f(file_id,'velZ',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        ! Select hyperslab in the file space / memory space:
!        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
!        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
!
!        ! Create property list for collective dataset write:
!        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
!        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
!
!        ! Write the dataset collectively:
!        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,3),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
!
!        call h5pclose_f(plist_id ,herror)
!        call h5dclose_f(dset_id  ,herror)
!        call h5sclose_f(memspace ,herror)
!        call h5sclose_f(filespace,herror)
!        !-----------------------------------------------------------------------------------------------------------
!        !-----------------------------------------------------------------------------------------------------------
!        if (1 == 1) then ! TEST!!! wird u.a. gebraucht, um 'velabs' zu speichern.
!            ! Create file space / memory space:
!            call h5screate_simple_f(3,dims_file,filespace,herror)
!            call h5screate_simple_f(3,dims_mem ,memspace ,herror)
!            !-----------------------------------------------------------------------------------------------------------
!            ! Create the dataset:
!            call h5dcreate_f(file_id,'pre',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
!            !-----------------------------------------------------------------------------------------------------------
!            ! Select hyperslab in the file space / memory space:
!            call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
!            call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
!
!            ! Create property list for collective dataset write:
!            call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
!            call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
!
!            ! Write the dataset collectively:
!            call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,pre,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
!
!            call h5pclose_f(plist_id ,herror)
!            call h5dclose_f(dset_id  ,herror)
!            call h5sclose_f(memspace ,herror)
!            call h5sclose_f(filespace,herror)
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        !-----------------------------------------------------------------------------------------------------------
!        if (vel_dir == 1) then
!            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1u(S1w:M1w))
!        else
!            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
!        end if
!        if (vel_dir == 2) then
!            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorY',array=y2v(S2w:M2w))
!        else
!            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorY',array=y2p(S2w:M2w))
!        end if
!        if (vel_dir == 3) then
!            call write_hdf_infoREAL(dim3,.false.,.false.,'VectorZ',array=y3w(S3w:M3w))
!        else
!            call write_hdf_infoREAL(dim3,.false.,.false.,'VectorZ',array=y3p(S3w:M3w))
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        call h5fclose_f(file_id,herror)
!    !===========================================================================================================
!
!
!    end subroutine write_hdf_velall
  
  
  
  
  
  
  
  
  
  
    !> \todo include ls1, BC_global, IB, y1p ...
    ! TEST!!! evtl. schoener: write_hdf_2D??
    ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
    subroutine write_2D_hdf(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,dir,phi)
  
        implicit none
  
        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname
  
        integer, intent(in)      ::  N1
        integer, intent(in)      ::  N2
  
        integer, intent(in)      ::  SS1
        integer, intent(in)      ::  SS2
  
        integer, intent(in)      ::  NN1
        integer, intent(in)      ::  NN2
  
        integer, intent(in)      ::  iShift
        integer, intent(in)      ::  jShift
  
        integer, intent(in)      ::  dir
  
        real   , intent(in)      ::  phi(1:N1,1:N2)
  
        integer                  ::  S1w, M1w, dim1
        integer                  ::  S2w, M2w, dim2
  
        logical                  ::  attr_yes
  
        integer(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
        integer(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
        !                eingefuehrt.                                                                              !
        !----------------------------------------------------------------------------------------------------------!
  
  
        !===========================================================================================================
        !=== Attribut-Datentyp =====================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Dimensionen und Offsets ===============================================================================
        !===========================================================================================================
        if (ABS(dir) == 1) then
            S1w = 2  + ls2
            S2w = 2  + ls3
     
            M1w = M2 + ls2
            M2w = M3 + ls3
     
            if (BC_2L_global /= -1) then
                S1w = 1
                M1w = M2
            end if
            if (BC_3L_global /= -1) then
                S2w = 1
                M2w = M3
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (ABS(dir) == 2) then
            S1w = 2  + ls1
            S2w = 2  + ls3
     
            M1w = M1 + ls1
            M2w = M3 + ls3
     
            if (BC_1L_global /= -1) then
                S1w = 1
                M1w = M1
            end if
            if (BC_3L_global /= -1) then
                S2w = 1
                M2w = M3
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (ABS(dir) == 3) then
            S1w = 2  + ls1
            S2w = 2  + ls2
     
            M1w = M1 + ls1
            M2w = M2 + ls2
     
            if (BC_1L_global /= -1) then
                S1w = 1
                M1w = M1
            end if
            if (BC_2L_global /= -1) then
                S2w = 1
                M2w = M2
            end if
        end if
        !===========================================================================================================
        dim1 = M1w-S1w+1
        dim2 = M2w-S2w+1
  
        dims_file   = (/ dim1 , dim2 /)
        dims_mem    = (/ N1   , N2   /)
        dims_data   = (/(NN1-SS1+1),(NN2-SS2+1)/)
        offset_mem  = (/ SS1-1, SS2-1/)
        offset_file = (/iShift,jShift/)
  
        if (.not. ((dir == -1 .and. iB(1,1) == 1  ) .or. (dir == -2 .and. iB(2,1) == 1  ) .or. (dir == -3 .and. iB(3,1) == 1   ) .or.    &
            &    (dir ==  1 .and. iB(1,1) == NB1) .or. (dir ==  2 .and. iB(2,1) == NB2) .or. (dir ==  3 .and. iB(3,1) == NB3))) then
            dims_data   = (/0,0/)
            offset_mem  = (/0,0/)
            offset_file = (/0,0/)
        end if
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== HDF5 schreiben ========================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access:
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Create the file collectively:
        call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
        call h5pclose_f(plist_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Create file space / memory space:
        call h5screate_simple_f(2,dims_file,filespace,herror)
        call h5screate_simple_f(2,dims_mem ,memspace ,herror)
  
        ! Create the dataset:
        call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Write the attributes:
        attr_yes = .true.
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'restart'          ,scalar=restart          )
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
        call write_hdf_infoINT (2     ,.false.,attr_yes,'S1w S2w'          ,array =(/S1w,S2w/)      )
        call write_hdf_infoINT (2     ,.false.,attr_yes,'M1w M2w'          ,array =(/M1w,M2w/)      )
  
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time'             ,scalar=time             )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )
  
        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'concentration_yes',scalar=concentration_yes)
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'n_conc'           ,scalar=n_conc           )
        call write_hdf_infoREAL(3     ,.false.,attr_yes,'gravity'          ,array =gravity          )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'Re'               ,scalar=Re               )
        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Sc'               ,array =Sc               )
        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Ric'              ,array =Ric              )
        call write_hdf_infoREAL(n_conc,.false.,attr_yes,'usc'              ,array =usc              )
        !-----------------------------------------------------------------------------------------------------------
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
        ! Create property list for collective dataset write:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
        ! Write the dataset collectively:
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
        call h5pclose_f(plist_id ,herror)
        call h5dclose_f(dset_id  ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
        !-----------------------------------------------------------------------------------------------------------
        if (ABS(dir) == 1) then
            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorY',array=y2p(S1w:M1w))
            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorZ',array=y3p(S2w:M2w))
        end if
        if (ABS(dir) == 2) then
            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorZ',array=y3p(S2w:M2w))
        end if
        if (ABS(dir) == 3) then
            call write_hdf_infoREAL(dim1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
            call write_hdf_infoREAL(dim2,.false.,.false.,'VectorY',array=y2p(S2w:M2w))
        end if
        !-----------------------------------------------------------------------------------------------------------
        call h5fclose_f(file_id,herror)
    !===========================================================================================================
  
  
    end subroutine write_2D_hdf
  
  
  
  
  
  
  
  
  
  
  
    subroutine write_1D_hdf(filename,dsetname,NN,SSp,NNp,iShift,dir,phi)
  
        implicit none
  
        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname
  
        integer, intent(in)      ::  NN
        integer, intent(in)      ::  SSp
        integer, intent(in)      ::  NNp
        integer, intent(in)      ::  iShift
        integer, intent(in)      ::  dir
        real   , intent(in)      ::  phi(1:NN)
  
        integer                  ::  SSw, MMw, dims
  
        logical                  ::  attr_yes
  
        integer(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
        integer(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
        !                eingefuehrt.                                                                              !
        !----------------------------------------------------------------------------------------------------------!
  
  
        !===========================================================================================================
        !=== Attribut-Datentyp =====================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Dimensionen und Offsets ===============================================================================
        !===========================================================================================================
        if (dir == 1) then
            SSw = 2  + ls1
            MMw = M1 + ls1
            if (BC_1L_global /= -1) then
                SSw = 1
                MMw = M1
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (dir == 2) then
            SSw = 2  + ls2
            MMw = M2 + ls2
            if (BC_2L_global /= -1) then
                SSw = 1
                MMw = M2
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (dir == 3) then
            SSw = 2  + ls3
            MMw = M3 + ls3
            if (BC_3L_global /= -1) then
                SSw = 1
                MMw = M3
            end if
        end if
        !===========================================================================================================
        dims = MMw-SSw+1
  
        dims_file   = (/dims     /)
        dims_mem    = (/NN       /)
        dims_data   = (/NNp-SSp+1/)
        offset_mem  = (/SSp-1    /)
        offset_file = (/iShift   /)
  
        if (.not. ((dir == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) .or.     &
            &    (dir == 2 .and. iB(1,1) == 1 .and. iB(3,1) == 1) .or.     &
            &    (dir == 3 .and. iB(1,1) == 1 .and. iB(2,1) == 1))) then
            dims_data   = (/0/)
            offset_mem  = (/0/)
            offset_file = (/0/)
        end if
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== HDF5 schreiben ========================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access:
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Create the file collectively:
        call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
        call h5pclose_f(plist_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Create file space / memory space:
        call h5screate_simple_f(1,dims_file,filespace,herror)
        call h5screate_simple_f(1,dims_mem ,memspace ,herror)
  
        ! Create the dataset:
        call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Write the attributes:
        attr_yes = .true.
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'SSw'              ,scalar=SSw              ) ! TEST!!! SSw, MMw unschoen ...
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'MMw'              ,scalar=MMw              )
        !-----------------------------------------------------------------------------------------------------------
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
        ! Create property list for collective dataset write:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
        ! Write the dataset collectively:
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
        call h5pclose_f(plist_id ,herror)
        call h5dclose_f(dset_id  ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
        !-----------------------------------------------------------------------------------------------------------
        if (dir == 1) call write_hdf_infoREAL(dims,.false.,.false.,'VectorX',array=y1p(SSw:MMw))
        if (dir == 2) call write_hdf_infoREAL(dims,.false.,.false.,'VectorY',array=y2p(SSw:MMw))
        if (dir == 3) call write_hdf_infoREAL(dims,.false.,.false.,'VectorZ',array=y3p(SSw:MMw))
        !-----------------------------------------------------------------------------------------------------------
        call h5fclose_f(file_id,herror)
    !===========================================================================================================
  
  
    end subroutine write_1D_hdf
  
  
  
  
  
  
  
  
  
  
  
    subroutine read_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
        implicit none
  
        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname
  
        integer, intent(in)      ::  SS1
        integer, intent(in)      ::  SS2
        integer, intent(in)      ::  SS3
  
        integer, intent(in)      ::  NN1
        integer, intent(in)      ::  NN2
        integer, intent(in)      ::  NN3
  
        integer, intent(in)      ::  vel_dir
  
        real   , intent(out)     ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
        integer(HSSIZE_T)        ::  offset_mem(1:3)
  
        integer                  ::  S1w, M1w, dim1
        integer                  ::  S2w, M2w, dim2
        integer                  ::  S3w, M3w, dim3
  
        logical                  ::  attr_yes
  
  
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
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Dimensionen und Offsets ===============================================================================
        !===========================================================================================================
        call filespace_props(vel_dir,(/1,1,1/),S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
        dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
        dims_data  = (/(NN1-SS1+1   ),(NN2-SS2+1   ),(NN3-SS3+1   )/)
        offset_mem = (/(SS1-b1L     ),(SS2-b2L     ),(SS3-b3L     )/)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== HDF5 lesen ============================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Open the file collectively
        call h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
            call MPI_FINALIZE(merror)
            stop
        end if
  
        call h5pclose_f(plist_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Open the dataset:
        call h5dopen_f(file_id,dsetname,dset_id,herror)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
            call MPI_FINALIZE(merror)
            stop
        end if
        !-----------------------------------------------------------------------------------------------------------
        ! Read the attributes:
        attr_yes = .true.
        call read_hdf_infoREAL(1,.true.,attr_yes,'time'          ,scalar=time          )
        call read_hdf_infoREAL(1,.true.,attr_yes,'dtime'         ,scalar=dtime         )
        call read_hdf_infoREAL(1,.true.,attr_yes,'time_out_vect' ,scalar=time_out_vect )
        call read_hdf_infoREAL(1,.true.,attr_yes,'time_out_scal' ,scalar=time_out_scal )
        call read_hdf_infoINT (1,.true.,attr_yes,'timestep'      ,scalar=timestep      )
        call read_hdf_infoINT (1,.true.,attr_yes,'write_count'   ,scalar=write_count   )
        call read_hdf_infoLOG (1,.true.,attr_yes,'write_out_vect',scalar=write_out_vect)
        call read_hdf_infoLOG (1,.true.,attr_yes,'write_out_scal',scalar=write_out_scal)
        call read_hdf_infoLOG (1,.true.,attr_yes,'new_dtime'     ,scalar=new_dtime     )
        call read_hdf_infoINT (1,.true.,attr_yes,'n_conc'        ,scalar=n_conc_old    )
        !-----------------------------------------------------------------------------------------------------------
        ! Get file space / create memory space:
        call h5dget_space_f    (dset_id   ,filespace,herror)
        call h5screate_simple_f(3,dims_mem,memspace ,herror)
  
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
        ! Create property list for collective dataset read:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
        ! Read the dataset collectively:
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
        call h5pclose_f(plist_id ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
        call h5dclose_f(dset_id  ,herror)
        !-----------------------------------------------------------------------------------------------------------
        call h5fclose_f(file_id,herror)
        !===========================================================================================================
  
  
        ! TEST!!! neu und ungetestet:
        !===========================================================================================================
        !=== Ghost-cell update =====================================================================================
        !===========================================================================================================
        call exchange(1,vel_dir,phi)
        call exchange(2,vel_dir,phi)
        call exchange(3,vel_dir,phi)
    !===========================================================================================================
  
  
    end subroutine read_hdf
  
  
  
  
  
  
  
  
  
  
  
    subroutine read2_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
        implicit none
  
        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname
  
        integer, intent(in)      ::  SS1
        integer, intent(in)      ::  SS2
        integer, intent(in)      ::  SS3
  
        integer, intent(in)      ::  NN1
        integer, intent(in)      ::  NN2
        integer, intent(in)      ::  NN3
  
        integer, intent(in)      ::  vel_dir
  
        real   , intent(out)     ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
        integer(HSSIZE_T)        ::  offset_mem(1:3)
  
        integer                  ::  S1w, M1w, dim1
        integer                  ::  S2w, M2w, dim2
        integer                  ::  S3w, M3w, dim3
  
        integer                  ::  S1r, N1r, i0, iGrid
        integer                  ::  S2r, N2r, j0, jGrid
        integer                  ::  S3r, N3r, k0, kGrid
  
        integer                  ::  Siw(1:3), Miw(1:3)
  
        logical                ::  attr_yes
  
  
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
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Open the file collectively
        call h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
            call MPI_FINALIZE(merror)
            stop
        end if
  
        call h5pclose_f(plist_id,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Attribute lesen =======================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Open the dataset:
        call h5dopen_f(file_id,dsetname,dset_id,herror)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
            call MPI_FINALIZE(merror)
            stop
        end if
        !-----------------------------------------------------------------------------------------------------------
        ! Read the attributes:
        attr_yes = .true.
        call read_hdf_infoREAL(1,.true. ,attr_yes,'time'          ,scalar=time          )
        call read_hdf_infoREAL(1,.true. ,attr_yes,'dtime'         ,scalar=dtime         )
        call read_hdf_infoREAL(1,.true. ,attr_yes,'time_out_vect' ,scalar=time_out_vect )
        call read_hdf_infoREAL(1,.true. ,attr_yes,'time_out_scal' ,scalar=time_out_scal )
        call read_hdf_infoINT (1,.true. ,attr_yes,'timestep'      ,scalar=timestep      )
        call read_hdf_infoINT (1,.true. ,attr_yes,'write_count'   ,scalar=write_count   )
        call read_hdf_infoLOG (1,.true. ,attr_yes,'write_out_vect',scalar=write_out_vect)
        call read_hdf_infoLOG (1,.true. ,attr_yes,'write_out_scal',scalar=write_out_scal)
        call read_hdf_infoLOG (1,.true. ,attr_yes,'new_dtime'     ,scalar=new_dtime     )
        call read_hdf_infoINT (3,.false.,attr_yes,'S1w S2w S3w'   ,array =Siw           )
        call read_hdf_infoINT (3,.false.,attr_yes,'M1w M2w M3w'   ,array =Miw           )
        call read_hdf_infoINT (1,.true. ,attr_yes,'n_conc'        ,scalar=n_conc_old    )
        !-----------------------------------------------------------------------------------------------------------
        call h5dclose_f(dset_id,herror)
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
        if (vel_dir == 1 .and. iB(1,1) > 1) then
            if (BC_1L_global  == -2 .and. ls1 ==  0) iGrid = 1
            if (BC_1L_global > 0 .and. ls1 ==  0) iGrid = 2
            if (BC_1L_global > 0 .and. ls1 == -1) iGrid = 1
        end if
        if (vel_dir == 2 .and. iB(2,1) > 1) then
            if (BC_2L_global  == -2 .and. ls2 ==  0) jGrid = 1
            if (BC_2L_global > 0 .and. ls2 ==  0) jGrid = 2
            if (BC_2L_global > 0 .and. ls2 == -1) jGrid = 1
        end if
        if (vel_dir == 3 .and. iB(3,1) > 1) then
            if (BC_3L_global  == -2 .and. ls3 ==  0) kGrid = 1
            if (BC_3L_global > 0 .and. ls3 ==  0) kGrid = 2
            if (BC_3L_global > 0 .and. ls3 == -1) kGrid = 1
        end if
        !-----------------------------------------------------------------------------------------------------------
        if ((S1w+i0) <= (NN1+iShift) .and. (M1w+i0) >= (SS1+iShift)) then
            if ((S1w+i0) >= (SS1+iShift)) then
                S1r = S1w+i0-iShift
            else
                S1r = SS1
            end if
            if ((M1w+i0) <= (NN1+iShift)) then
                N1r = M1w+i0-iShift
            else
                N1r = NN1
            end if
     
            dims_data  (1) = N1r-S1r+1
            offset_mem (1) = S1r-b1L
            offset_file(1) = iShift+iGrid-i0
     
            if (offset_file(1) < 0   ) offset_file(1) = 0
            if (offset_file(1) > dim1) offset_file(1) = dim1
        else
            dims_data  (1) = 0
            offset_mem (1) = 0
            offset_file(1) = 0
        end if
        !-----------------------------------------------------------------------------------------------------------
        if ((S2w+j0) <= (NN2+jShift) .and. (M2w+j0) >= (SS2+jShift)) then
            if ((S2w+j0) >= (SS2+jShift)) then
                S2r = S2w+j0-jShift
            else
                S2r = SS2
            end if
            if ((M2w+j0) <= (NN2+jShift)) then
                N2r = M2w+j0-jShift
            else
                N2r = NN2
            end if
     
            dims_data  (2) = N2r-S2r+1
            offset_mem (2) = S2r-b2L
            offset_file(2) = jShift+jGrid-j0
     
            if (offset_file(2) < 0   ) offset_file(2) = 0
            if (offset_file(2) > dim2) offset_file(2) = dim2
        else
            dims_data  (2) = 0
            offset_mem (2) = 0
            offset_file(2) = 0
        end if
        !-----------------------------------------------------------------------------------------------------------
        if ((S3w+k0) <= (NN3+kShift) .and. (M3w+k0) >= (SS3+kShift)) then
            if ((S3w+k0) >= (SS3+kShift)) then
                S3r = S3w+k0-kShift
            else
                S3r = SS3
            end if
            if ((M3w+k0) <= (NN3+kShift)) then
                N3r = M3w+k0-kShift
            else
                N3r = NN3
            end if
     
            dims_data  (3) = N3r-S3r+1
            offset_mem (3) = S3r-b3L
            offset_file(3) = kShift+kGrid-k0
     
            if (offset_file(3) < 0   ) offset_file(3) = 0
            if (offset_file(3) > dim3) offset_file(3) = dim3
        else
            dims_data  (3) = 0
            offset_mem (3) = 0
            offset_file(3) = 0
        end if
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        !=== Feld lesen ============================================================================================
        !===========================================================================================================
        phi = 0.
        !-----------------------------------------------------------------------------------------------------------
        ! Open the dataset:
        call h5dopen_f(file_id,dsetname,dset_id,herror)
  
        ! Get file space / create memory space:
        call h5dget_space_f    (dset_id   ,filespace,herror)
        call h5screate_simple_f(3,dims_mem,memspace ,herror)
  
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
        ! Create property list for collective dataset read:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
        ! Read the dataset collectively:
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
        call h5pclose_f(plist_id ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
        call h5dclose_f(dset_id  ,herror)
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        !=== File schliessen =======================================================================================
        !===========================================================================================================
        call h5fclose_f(file_id,herror)
        !===========================================================================================================
  
  
        ! TEST!!! neu und ungetestet:
        !===========================================================================================================
        !=== Ghost-cell update =====================================================================================
        !===========================================================================================================
        call exchange(1,vel_dir,phi)
        call exchange(2,vel_dir,phi)
        call exchange(3,vel_dir,phi)
    !===========================================================================================================
  
  
    end subroutine read2_hdf
  
  
  
  
  
  
  
  
  
  
    ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
    subroutine read2_2D_hdf(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,dir,phi)
  
        implicit none
  
        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname
  
        integer, intent(in)      ::  N1
        integer, intent(in)      ::  N2
  
        integer, intent(in)      ::  SS1
        integer, intent(in)      ::  SS2
  
        integer, intent(in)      ::  NN1
        integer, intent(in)      ::  NN2
  
        integer, intent(in)      ::  iShift
        integer, intent(in)      ::  jShift
  
        integer, intent(in)      ::  dir
  
        real   , intent(out)     ::  phi(1:N1,1:N2)
  
        integer(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
        integer(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
        integer                  ::  S1w, M1w, dim1
        integer                  ::  S2w, M2w, dim2
  
        integer                  ::  S1r, N1r, i0
        integer                  ::  S2r, N2r, j0
  
        integer                  ::  Siw(1:2), Miw(1:2)
  
        logical                  ::  attr_yes
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen:                                                                                             !
        !----------------------------------------------------------------------------------------------------------!
  
  
        !===========================================================================================================
        !=== File oeffnen ==========================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Open the file collectively
        call h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
            call MPI_FINALIZE(merror)
            stop
        end if
  
        call h5pclose_f(plist_id,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Attribute lesen =======================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Open the dataset:
        call h5dopen_f(file_id,dsetname,dset_id,herror)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
            call MPI_FINALIZE(merror)
            stop
        end if
        !-----------------------------------------------------------------------------------------------------------
        ! Read the attributes:
        attr_yes = .true.
        call read_hdf_infoREAL(1     ,.true. ,attr_yes,'time'          ,scalar=time          )
        call read_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime'         ,scalar=dtime         )
        call read_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_vect' ,scalar=time_out_vect )
        call read_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_scal' ,scalar=time_out_scal )
        call read_hdf_infoINT (1     ,.true. ,attr_yes,'timestep'      ,scalar=timestep      )
        call read_hdf_infoINT (1     ,.true. ,attr_yes,'write_count'   ,scalar=write_count   )
        call read_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_vect',scalar=write_out_vect)
        call read_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_scal',scalar=write_out_scal)
        call read_hdf_infoLOG (1     ,.true. ,attr_yes,'new_dtime'     ,scalar=new_dtime     )
        call read_hdf_infoINT (2     ,.false.,attr_yes,'S1w S2w'       ,array =Siw           )
        call read_hdf_infoINT (2     ,.false.,attr_yes,'M1w M2w'       ,array =Miw           )
        call read_hdf_infoINT (1     ,.true. ,attr_yes,'n_conc'        ,scalar=n_conc_old    )
        !-----------------------------------------------------------------------------------------------------------
        call h5dclose_f(dset_id,herror)
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
        if ((S1w+i0) <= (NN1+iShift) .and. (M1w+i0) >= (SS1+iShift)) then
            if ((S1w+i0) >= (SS1+iShift)) then
                S1r = S1w+i0-iShift
            else
                S1r = SS1
            end if
            if ((M1w+i0) <= (NN1+iShift)) then
                N1r = M1w+i0-iShift
            else
                N1r = NN1
            end if
     
            dims_data  (1) = N1r-S1r+1
            offset_mem (1) = S1r-1
            offset_file(1) = iShift-i0
     
            if (offset_file(1) < 0   ) offset_file(1) = 0
            if (offset_file(1) > dim1) offset_file(1) = dim1
        else
            dims_data  (1) = 0
            offset_mem (1) = 0
            offset_file(1) = 0
        end if
        !-----------------------------------------------------------------------------------------------------------
        if ((S2w+j0) <= (NN2+jShift) .and. (M2w+j0) >= (SS2+jShift)) then
            if ((S2w+j0) >= (SS2+jShift)) then
                S2r = S2w+j0-jShift
            else
                S2r = SS2
            end if
            if ((M2w+j0) <= (NN2+jShift)) then
                N2r = M2w+j0-jShift
            else
                N2r = NN2
            end if
     
            dims_data  (2) = N2r-S2r+1
            offset_mem (2) = S2r-1
            offset_file(2) = jShift-j0
     
            if (offset_file(2) < 0   ) offset_file(2) = 0
            if (offset_file(2) > dim2) offset_file(2) = dim2
        else
            dims_data  (2) = 0
            offset_mem (2) = 0
            offset_file(2) = 0
        end if
        !-----------------------------------------------------------------------------------------------------------
        ! ACTHUNG!!! Nicht ganz klar, ob dieses Vorgehen im Sinne des Erfinders ist:
        if (.not. ((dir == -1 .and. iB(1,1) == 1  ) .or. (dir == -2 .and. iB(2,1) == 1  ) .or. (dir == -3 .and. iB(3,1) == 1   ) .or.    &
            &    (dir ==  1 .and. iB(1,1) == NB1) .or. (dir ==  2 .and. iB(2,1) == NB2) .or. (dir ==  3 .and. iB(3,1) == NB3))) then
            dims_data   = (/0,0/)
            offset_mem  = (/0,0/)
            offset_file = (/0,0/)
        end if
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Feld lesen ============================================================================================
        !===========================================================================================================
        phi = 0.
        !-----------------------------------------------------------------------------------------------------------
        ! Open the dataset:
        call h5dopen_f(file_id,dsetname,dset_id,herror)
  
        ! Get file space / create memory space:
        call h5dget_space_f    (dset_id   ,filespace,herror)
        call h5screate_simple_f(2,dims_mem,memspace ,herror)
  
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
        ! Create property list for collective dataset read:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
        ! Read the dataset collectively:
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
        call h5pclose_f(plist_id ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
  
        call h5dclose_f(dset_id  ,herror)
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        !=== File schliessen =======================================================================================
        !===========================================================================================================
        call h5fclose_f(file_id,herror)
    !===========================================================================================================
  
  
    ! TEST!!! fehlt noch:
    !===========================================================================================================
    !=== Ghost-cell update =====================================================================================
    !===========================================================================================================
    !CALL exchange(1,vel_dir,phi)
    !CALL exchange(2,vel_dir,phi)
    !CALL exchange(3,vel_dir,phi)
    !===========================================================================================================
  
  
    end subroutine read2_2D_hdf
  
  
  
  
  
  
  
  
  
  
  
    subroutine read2_1D_hdf(filename,dsetname,NN,SSp,NNp,iShift,dir,phi)
  
        implicit none
  
        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname
  
        integer, intent(in)      ::  NN
        integer, intent(in)      ::  SSp
        integer, intent(in)      ::  NNp
        integer, intent(in)      ::  iShift
        integer, intent(in)      ::  dir
  
        real   , intent(out)     ::  phi(1:NN)
  
        integer(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
        integer(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
        integer                  ::  SSw, MMw, dims
        integer                  ::  SSr, NNr, i0
  
        logical                  ::  attr_yes
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen:                                                                                             !
        !----------------------------------------------------------------------------------------------------------!
  
  
        !===========================================================================================================
        !=== File oeffnen ==========================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Open the file collectively
        call h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
            call MPI_FINALIZE(merror)
            stop
        end if
  
        call h5pclose_f(plist_id,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Attribute lesen =======================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Open the dataset:
        call h5dopen_f(file_id,dsetname,dset_id,herror)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
            call MPI_FINALIZE(merror)
            stop
        end if
        !-----------------------------------------------------------------------------------------------------------
        ! Read the attributes:
        attr_yes = .true.
        call read_hdf_infoINT (1     ,.true. ,attr_yes,'SSw'           ,scalar=SSw           ) ! TEST!!! SSw, MMw unschoen ...
        call read_hdf_infoINT (1     ,.true. ,attr_yes,'MMw'           ,scalar=MMw           )
        !-----------------------------------------------------------------------------------------------------------
        call h5dclose_f(dset_id,herror)
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
        if ((SSw+i0) <= (NNp+iShift) .and. (MMw+i0) >= (SSp+iShift)) then
            if ((SSw+i0) >= (SSp+iShift)) then
                SSr = SSw+i0-iShift
            else
                SSr = SSp
            end if
            if ((MMw+i0) <= (NNp+iShift)) then
                NNr = MMw+i0-iShift
            else
                NNr = NNp
            end if
     
            dims_data  (1) = NNr-SSr+1
            offset_mem (1) = SSr-1
            offset_file(1) = iShift-i0
     
            if (offset_file(1) < 0   ) offset_file(1) = 0
            if (offset_file(1) > dims) offset_file(1) = dims
        else
            dims_data  (1) = 0
            offset_mem (1) = 0
            offset_file(1) = 0
        end if
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
        call h5dopen_f(file_id,dsetname,dset_id,herror)
  
        ! Get file space / create memory space:
        call h5dget_space_f    (dset_id   ,filespace,herror)
        call h5screate_simple_f(1,dims_mem,memspace ,herror)
  
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
        ! Create property list for collective dataset read:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
        ! Read the dataset collectively:
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
        call h5pclose_f(plist_id ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
  
        call h5dclose_f(dset_id  ,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== File schliessen =======================================================================================
        !===========================================================================================================
        call h5fclose_f(file_id,herror)
    !===========================================================================================================
  
  
    end subroutine read2_1D_hdf
  
  
  
  
  
  
  
  
  
  
  
    subroutine write_part_hdf(filename,all_args_yes)
  
        implicit none
  
        character(*), intent(in) ::  filename
        logical     , intent(in) ::  all_args_yes
  
        logical                  ::  attr_yes
  
        integer(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
        integer(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
        integer(HID_T )          ::  group_id
        integer(SIZE_T)          ::  size_hint
  
        integer                  ::  proc_table       (1:NB1*NB2*NB3)
        integer                  ::  proc_table_global(1:NB1*NB2*NB3)
  
        integer                  ::  m, g
        integer                  ::  dims    (1:n_groups)
        integer                  ::  dims_all(1:n_groups)
        integer                  ::  dims_prev, dims_global, dims_buf, offset
  
        real                     ::  bufferA(1:n_args)
        real   , allocatable     ::  bufferB(:)
  
        character(len=1)         ::  groupNr
  
        integer                  ::  S1w, M1w
        integer                  ::  S2w, M2w
        integer                  ::  S3w, M3w
  
  
        !--- Tuning-Parameter (noch nicht optimiert!) ---
        size_hint = 20
  
  
        !===========================================================================================================
        !=== Attribut-Datentyp =====================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== HDF5 schreiben ========================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access:
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        ! Create the file collectively:
        call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
        call h5pclose_f(plist_id,herror)
        !-----------------------------------------------------------------------------------------------------------
  
        !-----------------------------------------------------------------------------------------------------------
        ! Create property list for collective dataset write:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)  ! TEST!!! Umsortiert! Generell besser so? Koennte auch noch weiter nach oben ...
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
        !-----------------------------------------------------------------------------------------------------------
  
        dims_prev = 0
        dims      = 0
  
        do g = 1, n_groups
     
            !========================================================================================================
            ! Create a group:
            write(groupNr,'(i1.1)') g
            call h5gcreate_f(file_id,'group'//groupNr,group_id,herror,size_hint)
            !========================================================================================================
     
     
            !========================================================================================================
            !=== Dimensionen und Offsets ============================================================================
            !========================================================================================================
            !--- Nach Gruppe sortieren ---
            do m = dims_prev+1, n_part
                if (NINT(particles(14,m)) == g) then
                    dims(g) = dims(g) + 1
           
                    bufferA  (1:n_args                  ) = particles(1:n_args,dims_prev+dims(g))
                    particles(1:n_args,dims_prev+dims(g)) = particles(1:n_args,m                )
                    particles(1:n_args,m                ) = bufferA  (1:n_args                  )
                end if
            end do
            !--------------------------------------------------------------------------------------------------------
            proc_table = 0
            proc_table(rank+1) = dims(g)
            call MPI_ALLREDUCE(proc_table,proc_table_global,NB1*NB2*NB3,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
     
            dims_global = 0
            do m = 1, NB1*NB2*NB3
                dims_global = dims_global + proc_table_global(m)
            end do
     
            offset = 0
            do m = 1, rank
                offset = offset + proc_table_global(m)
            end do
     
            dims_buf = dims(g)
            if (dims_buf == 0) dims_buf = 1
            !========================================================================================================
     
     
            if (dims_global > 0) then
        
                !=====================================================================================================
                dims_file   = (/dims_global/)
                dims_mem    = (/dims_buf   /)
                dims_data   = (/dims(g)    /)
                offset_mem  = (/0          /)
                offset_file = (/offset     /)
                !=====================================================================================================
        
        
                !=====================================================================================================
                ! Create file space / memory space:
                call h5screate_simple_f(1,dims_file,filespace,herror)
                call h5screate_simple_f(1,dims_mem ,memspace ,herror)
        
                ! Select hyperslab in the file space / memory space:
                call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror) ! TEST!!! Umsortiert! Generell besser so?
                call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
                !=====================================================================================================
        
                allocate(bufferB(1:dims_buf))
                !=====================================================================================================
                !--- create the dataset ---
                call h5dcreate_f(group_id,'X',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
                !--- umspeichern ---
                bufferB(1:dims_buf) = particles(1,(dims_prev+1):(dims_prev+dims_buf))
        
                ! Write the dataset collectively:
                call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
                !-----------------------------------------------------------------------------------------------------
                !--- create the dataset ---
                call h5dcreate_f(group_id,'Y',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
                !--- umspeichern ---
                bufferB(1:dims_buf) = particles(2,(dims_prev+1):(dims_prev+dims_buf))
        
                ! Write the dataset collectively:
                call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
                !-----------------------------------------------------------------------------------------------------
                if (dimens == 3) then
                    !--- create the dataset ---
                    call h5dcreate_f(group_id,'Z',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
                    !--- umspeichern ---
                    bufferB(1:dims_buf) = particles(3,(dims_prev+1):(dims_prev+dims_buf))
        
                    ! Write the dataset collectively:
                    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
                end if
                !-----------------------------------------------------------------------------------------------------
                !--- create the dataset ---
                call h5dcreate_f(group_id,'velX',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
                !--- umspeichern ---
                bufferB(1:dims_buf) = particles(4,(dims_prev+1):(dims_prev+dims_buf))
        
                ! Write the dataset collectively:
                call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
                !-----------------------------------------------------------------------------------------------------
                !--- create the dataset ---
                call h5dcreate_f(group_id,'velY',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
                !--- umspeichern ---
                bufferB(1:dims_buf) = particles(5,(dims_prev+1):(dims_prev+dims_buf))
        
                ! Write the dataset collectively:
                call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
                !-----------------------------------------------------------------------------------------------------
                if (dimens == 3) then
                    !--- create the dataset ---
                    call h5dcreate_f(group_id,'velZ',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        
                    !--- umspeichern ---
                    bufferB(1:dims_buf) = particles(6,(dims_prev+1):(dims_prev+dims_buf))
        
                    ! Write the dataset collectively:
                    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
                end if
                !-----------------------------------------------------------------------------------------------------
                !--- create the dataset ---
                call h5dcreate_f(group_id,'Number',H5T_NATIVE_DOUBLE,filespace,dset_id,herror) ! TEST!!! Evtl. gleich als Integer schreiben?
        
                !--- umspeichern ---
                bufferB(1:dims_buf) = particles(13,(dims_prev+1):(dims_prev+dims_buf))
        
                ! Write the dataset collectively:
                call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
                !-----------------------------------------------------------------------------------------------------
                if (all_args_yes) then
                    !--------------------------------------------------------------------------------------------------
                    !--- create the dataset ---
                    call h5dcreate_f(group_id,'fb_force',H5T_NATIVE_DOUBLE,filespace,dset_id,herror) ! TEST!!! brauchts das wirklich??
           
                    !--- umspeichern ---
                    bufferB(1:dims_buf) = particles(10,(dims_prev+1):(dims_prev+dims_buf))
           
                    ! Write the dataset collectively:
                    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
                    !--------------------------------------------------------------------------------------------------
                    !--- create the dataset ---
                    call h5dcreate_f(group_id,'usp',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
           
                    !--- umspeichern ---
                    bufferB(1:dims_buf) = particles(11,(dims_prev+1):(dims_prev+dims_buf))
           
                    ! Write the dataset collectively:
                    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
                    !--------------------------------------------------------------------------------------------------
                    !--- create the dataset ---
                    call h5dcreate_f(group_id,'Stp',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
           
                    !--- umspeichern ---
                    bufferB(1:dims_buf) = particles(12,(dims_prev+1):(dims_prev+dims_buf))
           
                    ! Write the dataset collectively:
                    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
                   !--------------------------------------------------------------------------------------------------
                end if
                !=====================================================================================================
                deallocate(bufferB)
        
                call h5sclose_f(memspace ,herror)
                call h5sclose_f(filespace,herror)
        
            end if
     
            call h5gclose_f(group_id ,herror)
     
            dims_prev = dims_prev + dims(g)
     
        end do
        !-----------------------------------------------------------------------------------------------------------
        call h5pclose_f(plist_id,herror)
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        call h5gcreate_f(file_id,'attributes',dset_id,herror,size_hint)
  
        call MPI_ALLREDUCE(dims,dims_all,n_groups,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
        ! Write the attributes:
        attr_yes = .true.
        call write_hdf_infoINT (1       ,.true. ,attr_yes,'restart'          ,scalar=restart          )
        call write_hdf_infoINT (1       ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
        call write_hdf_infoREAL(1       ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
        call write_hdf_infoREAL(1       ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
        call write_hdf_infoREAL(1       ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
        call write_hdf_infoREAL(1       ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
        call write_hdf_infoLOG (1       ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
        call write_hdf_infoLOG (1       ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
        call write_hdf_infoLOG (1       ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
        call write_hdf_infoREAL(1       ,.true. ,attr_yes,'time'             ,scalar=time             )
        call write_hdf_infoREAL(1       ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
        call write_hdf_infoINT (1       ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )
  
        call write_hdf_infoLOG (1       ,.true. ,attr_yes,'concentration_yes',scalar=concentration_yes)
        call write_hdf_infoINT (1       ,.true. ,attr_yes,'n_conc'           ,scalar=n_conc           )
        call write_hdf_infoREAL(3       ,.false.,attr_yes,'gravity'          ,array =gravity          )
        call write_hdf_infoREAL(1       ,.true. ,attr_yes,'Re'               ,scalar=Re               )
        call write_hdf_infoREAL(n_spec  ,.false.,attr_yes,'Rip'              ,array =Rip              )
        call write_hdf_infoREAL(n_spec  ,.false.,attr_yes,'usp'              ,array =usp              )
        call write_hdf_infoREAL(n_spec  ,.false.,attr_yes,'Stp'              ,array =Stp              )
  
        call write_hdf_infoINT (1       ,.true. ,attr_yes,'n_groups'         ,scalar=n_groups         )
        call write_hdf_infoINT (n_groups,.false.,attr_yes,'group dims'       ,array =dims_all         )
  
        call h5gclose_f(dset_id,herror)
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
  
        if (BC_1L_global /= -1) then
            S1w = 1
            M1w = M1
        end if
        if (BC_2L_global /= -1) then
            S2w = 1
            M2w = M2
        end if
        if (BC_3L_global /= -1) then
            S3w = 1
            M3w = M3
        end if
        !-----------------------------------------------------------------------------------------------------------
        call write_hdf_infoREAL(M1w-S1w+1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
        call write_hdf_infoREAL(M2w-S2w+1,.false.,.false.,'VectorY',array=y2p(S2w:M2w))
        if (dimens == 3) call write_hdf_infoREAL(M3w-S3w+1,.false.,.false.,'VectorZ',array=y3p(S3w:M3w))
        !===========================================================================================================
  
  
        call h5fclose_f(file_id ,herror)
    !===========================================================================================================
  
  
    end subroutine write_part_hdf
    
  
  
  
  
  
  
  
  
  
  
    subroutine write_part_serial(dirname)
  
        implicit none
  
        character(*), intent(in) ::  dirname
  
        character(len=3)         ::  bID(1:3)
  
        integer                  ::  i
  
        integer                  ::  S1w, M1w
        integer                  ::  S2w, M2w
        integer                  ::  S3w, M3w
  
  
        if (rank == 0) call SYSTEM('mkdir -p '//dirname) ! TEST!!! funktioniert das?
  
        call MPI_BARRIER(COMM_CART,merror)
  
        if (n_part >= 1) then
     
            call num_to_string(3,iB(1,1),bID(1))
            call num_to_string(3,iB(2,1),bID(2))
            call num_to_string(3,iB(3,1),bID(3))
     
            ! TEST!!! ASCII vs. binaer?
            !OPEN(98,FILE=dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt',STATUS='UNKNOWN',FORM='UNFORMATTED')
            open(98,file=dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt',status='UNKNOWN')
     
     
            !========================================================================================================
            !=== Schreiben ==========================================================================================
            !========================================================================================================
            write(98,'(1i12)') n_part
     
            do i = 1, n_part
                write(98,'(10E25.17)') particles(1:6,i), particles(10:14,i) ! TEST!!! ok?
            end do
            !========================================================================================================
     
     
            !========================================================================================================
            !=== Argumente ==========================================================================================
            !========================================================================================================
            ! TEST!!! ok? mehr Argumente?
            write(98,'(10i12    )') n_part, n_groups, restart, write_count, timestep, n_conc
            write(98,'(10E25.17 )') time, dtime, time_out_vect, time_out_scal, dtime_out_vect, dtime_out_scal
            write(98,'(10i12    )') NB1, NB2, NB3
            write(98,'(10E25.17 )') L1, L2, L3
            write(98,'(100E25.17)') Re, gravity(1:3), Rip(1:n_spec), usp(1:n_spec), Stp(1:n_spec)
            write(98,'(10L1     )') write_out_vect, write_out_scal, new_dtime, concentration_yes
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
     
            if (BC_1L_global /= -1) then
                S1w = 1
                M1w = M1
            end if
            if (BC_2L_global /= -1) then
                S2w = 1
                M2w = M2
            end if
            if (BC_3L_global /= -1) then
                S3w = 1
                M3w = M3
            end if
            !--------------------------------------------------------------------------------------------------------
            do i = S1w, M1w
                write(98,'(1E25.17)') y1p(i)
            end do
            do i = S2w, M2w
                write(98,'(1E25.17)') y2p(i)
            end do
            if (dimens == 3) then
                do i = S3w, M3w
                    write(98,'(1E25.17)') y3p(i)
                end do
            end if
            !========================================================================================================
     
            close(98)
     
        end if
  
  
    end subroutine write_part_serial
    
  
  
  
  
  
  
  
  
  
  
    subroutine read_part_hdf(filename,all_args_yes)
  
        implicit none
  
        character(*), intent(in) ::  filename
        logical     , intent(in) ::  all_args_yes
  
        logical                  ::  attr_yes
  
        !        integer(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
        integer(HSIZE_T )        ::  dims_file  (1)!, dims_mem  (1)!, dims_data(1)
        integer(HSIZE_T )        ::  maxdims    (1)
  
        integer(HID_T   )        ::  group_id
  
        integer                  ::  g !m,
        integer                  ::  n_part_all, dims, dims_prev, n_groups_old!, offset
  
        integer, allocatable     ::  dims_old(:)
        real   , allocatable     ::  bufferB(:)
  
        character(len=1)         ::  groupNr
  
  
        !===========================================================================================================
        !=== Attribut-Datentyp =====================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== HDF5 lesen ============================================================================================
        !===========================================================================================================
        !-- setup file access property list with parallel I/O access ---
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
        !--- open the file collectively ---
        call h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
            call MPI_FINALIZE(merror)
            stop
        end if
  
        call h5pclose_f(plist_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        call h5gopen_f(file_id,'attributes',dset_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        !--- read the attributes ---
        attr_yes = .true.
        call read_hdf_infoINT (1           ,.true. ,attr_yes,'restart'          ,scalar=restart          )
        call read_hdf_infoINT (1           ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
        call read_hdf_infoREAL(1           ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
        call read_hdf_infoREAL(1           ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
        call read_hdf_infoREAL(1           ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
        call read_hdf_infoREAL(1           ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
        call read_hdf_infoLOG (1           ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
        call read_hdf_infoLOG (1           ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
        call read_hdf_infoLOG (1           ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
        call read_hdf_infoREAL(1           ,.true. ,attr_yes,'time'             ,scalar=time             )
        call read_hdf_infoREAL(1           ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
        call read_hdf_infoINT (1           ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )
  
        call read_hdf_infoINT (1           ,.true. ,attr_yes,'n_groups'         ,scalar=n_groups_old     )
        !-----------------------------------------------------------------------------------------------------------
        if (n_groups_old < n_groups) then
            if (rank == 0) write(*,'(a,i3,a,i3,a)') 'Reading only n_groups_old =', n_groups_old, ' particle groups instead of n_groups ='    , n_groups     , ' ...'
        end if
        if (n_groups_old > n_groups) then
            if (rank == 0) write(*,'(a,i3,a,i3,a)') 'Reading only n_groups ='    , n_groups    , ' particle groups instead of n_groups_old =', n_groups_old , ' ...'
            n_groups_old = n_groups
        end if
        !-----------------------------------------------------------------------------------------------------------
        allocate(dims_old(1:n_groups_old))
        call read_hdf_infoINT (n_groups_old,.false.,attr_yes,'group dims'       ,array =dims_old         )
        !-----------------------------------------------------------------------------------------------------------
        n_part_all = 0
        do g = 1, n_groups_old
            n_part_all = n_part_all + dims_old(g)
        end do
  
        if (n_part_all < n_part_max) then
            if (rank == 0) write(*,'(a,i10,a,i10,a)') 'Reading only n_part_all =', n_part_all, ' particles instead of n_part_max =', n_part_max, ' ...'
        end if
        if (n_part_all > n_part_max) then
            if (rank == 0) write(*,'(a,i10,a,i10,a)') 'Reading only n_part_max =', n_part_max, ' particles instead of n_part_all =', n_part_all, ' ...'
        end if
        !-----------------------------------------------------------------------------------------------------------
        call h5gclose_f(dset_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Create property list for collective dataset write:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)  ! TEST!!! Umsortiert! Generell besser so? Koennte auch noch weiter nach oben ...
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_INDEPENDENT_F,herror)
        !-----------------------------------------------------------------------------------------------------------
  
        dims_prev = 0
        particles = 0.
  
        do g = 1, n_groups_old
     
            dims = dims_old(g)
            if (dims_prev+dims > n_part_max) dims = n_part_max - dims_prev
     
     
            if (dims > 0) then
        
                !=====================================================================================================
                !--- create a group ---
                write(groupNr,'(i1.1)') g
                call h5gopen_f(file_id,'group'//groupNr,group_id,herror)
                !=====================================================================================================
        
        
                allocate(bufferB(1:dims_old(g)))
                !=====================================================================================================
                !--- open the dataset ---
                call h5dopen_f(group_id,'X',dset_id,herror)
        
                call h5dget_space_f(dset_id,filespace,herror)
                call h5sget_simple_extent_dims_f(filespace,dims_file,maxdims,herror)
        
                !--- read the dataset independently ---
                call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
        
                !--- umspeichern ---
                particles(1,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                !-----------------------------------------------------------------------------------------------------
                !--- open the dataset ---
                call h5dopen_f(group_id,'Y',dset_id,herror)
        
                !--- read the dataset independently ---
                call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
        
                !--- umspeichern ---
                particles(2,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                !-----------------------------------------------------------------------------------------------------
                if (dimens == 3) then
                    !--- open the dataset ---
                    call h5dopen_f(group_id,'Z',dset_id,herror)
        
                    !--- read the dataset independently ---
                    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
        
                    !--- umspeichern ---
                    particles(3,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                end if
                !-----------------------------------------------------------------------------------------------------
                !--- open the dataset ---
                call h5dopen_f(group_id,'velX',dset_id,herror)
                !--- read the dataset independently ---
                call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
        
                !--- umspeichern ---
                particles(4,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                !-----------------------------------------------------------------------------------------------------
                !--- open the dataset ---
                call h5dopen_f(group_id,'velY',dset_id,herror)
                !--- read the dataset independently ---
                call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
        
                !--- umspeichern ---
                particles(5,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                !-----------------------------------------------------------------------------------------------------
                if (dimens == 3) then
                    !--- open the dataset ---
                    call h5dopen_f(group_id,'velZ',dset_id,herror)
                    !--- read the dataset independently ---
                    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
        
                    !--- umspeichern ---
                    particles(6,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                end if
                !-----------------------------------------------------------------------------------------------------
                !--- open the dataset ---
                call h5dopen_f(group_id,'Number',dset_id,herror)
        
                !--- read the dataset independently ---
                call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                call h5dclose_f(dset_id,herror)
        
                !--- umspeichern ---
                particles(13,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                !-----------------------------------------------------------------------------------------------------
                if (all_args_yes) then
                    !--------------------------------------------------------------------------------------------------
                    ! Open the dataset: ! TEST!!! brauchts das wirklich??
                    call h5dopen_f(group_id,'fb_force',dset_id,herror)
                    !--- read the dataset independently ---
                    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
           
                    !--- umspeichern ---
                    particles(10,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                    !--------------------------------------------------------------------------------------------------
                    ! Open the dataset:
                    call h5dopen_f(group_id,'usp'    ,dset_id,herror)
                    !--- read the dataset independently ---
                    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
           
                    !--- umspeichern ---
                    particles(11,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                    !--------------------------------------------------------------------------------------------------
                    ! Open the dataset:
                    call h5dopen_f(group_id,'Stp'    ,dset_id,herror)
                    !--- read the dataset independently ---
                    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,bufferB,dims_file,herror,xfer_prp=plist_id)
                    call h5dclose_f(dset_id,herror)
           
                    !--- umspeichern ---
                    particles(12,(dims_prev+1):(dims_prev+dims)) = bufferB(1:dims)
                   !--------------------------------------------------------------------------------------------------
                end if
                !=====================================================================================================
                deallocate(bufferB)
        
                particles(14,(dims_prev+1):(dims_prev+dims)) = REAL(g)
        
                call h5gclose_f(group_id,herror)
        
                dims_prev = dims_prev + dims
        
            end if
     
        end do
        !-----------------------------------------------------------------------------------------------------------
        call h5pclose_f(plist_id,herror)
        call h5fclose_f(file_id ,herror)
        !===========================================================================================================
        deallocate(dims_old)
  
  
    end subroutine read_part_hdf

  
  
  
  
  
  
  
  
  
  
    subroutine read_part_serial(dirname)
  
        implicit none
  
        character(*), intent(in) ::  dirname
  
        character(len=3)         ::  bID(1:3)
  
        integer                  ::  i
        integer                  ::  ios
  
  
        call num_to_string(3,iB(1,1),bID(1))
        call num_to_string(3,iB(2,1),bID(2))
        call num_to_string(3,iB(3,1),bID(3))
  
  
        ! TEST!!! ASCII vs. binaer?
        !OPEN(98,FILE=TRIM(dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt'),STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',IOSTAT=ios)
        open(98,file=TRIM(dirname//'/'//bID(1)//bID(2)//bID(3)//'.txt'),status='OLD',form='FORMATTED',action='READ',iostat=ios)
          
        if (ios == 0) then
     
            !===========================================================================================================
            !=== Lesen =================================================================================================
            !===========================================================================================================
            read(98,'(1i12)') n_part
     
            do i = 1, n_part
                read(98,'(10E25.17)') particles(1:6,i), particles(10:14,i) ! TEST!!! ok?
            end do
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
     
            close(98)
     
        else
     
            n_part = 0
     
        end if
  
  
    end subroutine read_part_serial

  
  
  
  
  

    subroutine write_stats_hdf_4D(filename,dsetname, SS1,SS2,SS3,SS4, NN1,NN2,NN3,NN4, phi)

        ! ---mjohn 230412

        ! function is suited to write 4D data, specifically designed for spatial mode analysis
        ! dims 1 and 2 are spatial directions, dim1 is x1, dim2 is x2, data is COLLOCATED (count from SSi to NNi).
        ! dims 3 and 4 are integer values, typically Fourier and Hermite modes (count from SSi to NNi).

        use mod_dims
        use mod_vars
        use mod_exchange
        use mod_lib
        use HDF5

        implicit none

        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname

        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3       ! supposed to be 1 (1:NN3)
        integer, intent(in)    ::  SS4       ! supposed to be 0 (0:NN4)

        integer, intent(in)    ::  NN1       ! supposed to be size in 1 direction (wall-normal)
        integer, intent(in)    ::  NN2       ! supposed to be size in 2 direction (spanwise)
        integer, intent(in)    ::  NN3       ! supposed to be n_frequencymodes (1:NN3)
        integer, intent(in)    ::  NN4       ! supposed to be n_Hermitemodes-1 (0:NN4)

        real   , intent(in)    ::  phi(SS1:NN1, SS2:NN2, SS3:NN3, SS4:NN4)

        integer                ::  Fourier(SS3:NN3), Hermite(SS4:NN4)
        integer                ::  stride(4)
        integer(HSIZE_T )      ::  dims_file  (4), dims_mem  (4), dims_data(4), stride_mem(4)
        integer(HSSIZE_T)      ::  offset_file(4), offset_mem(4)

        logical                ::  attr_yes
        integer                ::  i

  INCLUDE 'mpif.h'


        ! this function is for stride == 1 only
        stride(1:4) = 1

        !===========================================================================================================
        !=== Attribut-Datentyp =====================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
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
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

        ! Create the file collectively:
        call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
        call h5pclose_f(plist_id,herror)

        !-----------------------------------------------------------------------------------------------------------
        ! Create file space / memory space:
        call h5screate_simple_f(4, dims_file, filespace, herror)
        call h5screate_simple_f(4, dims_mem , memspace , herror)

        ! Create the dataset:
        call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
        !-----------------------------------------------------------------------------------------------------------
        ! Write the attributes:
        attr_yes = .true.
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'restart'          ,scalar=restart          )
        call write_hdf_infoINT (1     ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time'             ,scalar=time             )

        !-----------------------------------------------------------------------------------------------------------
        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror,stride_mem)

        ! Create property list for collective dataset write:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)

        ! Write the dataset collectively:
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)

        call h5pclose_f(plist_id ,herror)
        call h5dclose_f(dset_id  ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)

        do i = SS3, NN3
            Fourier(i) = i
        end do
        do i = SS4, NN4
            Hermite(i) = i
        end do

        call write_hdf_infoREAL(M1,       .false.,.false.,'VectorX',  array=y1p(SS1:NN1)     )
        call write_hdf_infoREAL(M2,       .false.,.false.,'VectorY',  array=y2p(SS2:NN2)     )
        call write_hdf_infoINT (NN3-SS3+1,.false.,.false.,'VectorDFT',array=Fourier(SS3:NN3) )
        call write_hdf_infoINT (NN4-SS4+1,.false.,.false.,'VectorHe', array=Hermite(SS4:NN4) )
        !-----------------------------------------------------------------------------------------------------------
        call h5fclose_f(file_id,herror)
    !===========================================================================================================


    end subroutine write_stats_hdf_4D




    subroutine read_stats_hdf_4D(filename,dsetname, SS1,SS2,SS3,SS4, NN1,NN2,NN3,NN4, phi)


        ! ---mjohn 230412

        ! function is suited to read 4D data, specifically designed for spatial mode analysis
        ! dims 1 and 2 are spatial directions, dim1 is x1, dim2 is x2, data is COLLOCATED (count from SSi to NNi).
        ! dims 3 and 4 are integer values, typically Fourier and Hermite modes (count from SSi to NNi).


        use mod_dims
        use mod_vars
        use mod_exchange
        use mod_lib
        use HDF5

        implicit none

        character(*), intent(in) ::  filename
        character(*), intent(in) ::  dsetname

        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3       ! supposed to be 1 (1:NN3)
        integer, intent(in)    ::  SS4       ! supposed to be 0 (0:NN4)

        integer, intent(in)    ::  NN1       ! supposed to be size in 1 direction (wall-normal)
        integer, intent(in)    ::  NN2       ! supposed to be size in 2 direction (spanwise)
        integer, intent(in)    ::  NN3       ! supposed to be n_frequencymodes (1:NN3)
        integer, intent(in)    ::  NN4       ! supposed to be n_Hermitemodes-1 (0:NN4)

        real   , intent(out)    ::  phi(SS1:NN1, SS2:NN2, SS3:NN3, SS4:NN4)

        integer(HSIZE_T )      ::  dims_file  (4), dims_mem  (4), dims_data(4)
        integer(HSSIZE_T)      ::  offset_file(4), offset_mem(4)

  INCLUDE 'mpif.h'



        !===========================================================================================================
        !=== File oeffnen ==========================================================================================
        !===========================================================================================================
        ! Setup file access property list with parallel I/O access
        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
        call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

        ! Open the file collectively
        call h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)

        if (herror == -1) then
            if (rank == 0) write(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
            call MPI_FINALIZE(merror)
            stop
        end if

        call h5pclose_f(plist_id,herror)
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
        call h5dopen_f(file_id,dsetname,dset_id,herror)

        ! Get file space / create memory space:
        call h5dget_space_f    (dset_id   ,filespace,herror)
        call h5screate_simple_f(4,dims_mem,memspace ,herror)

        ! Select hyperslab in the file space / memory space:
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
        call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)

        ! Create property list for collective dataset read:
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)

        ! Read the dataset collectively:
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)

        call h5pclose_f(plist_id ,herror)
        call h5sclose_f(memspace ,herror)
        call h5sclose_f(filespace,herror)
        call h5dclose_f(dset_id  ,herror)
        !===========================================================================================================

        !===========================================================================================================
        !=== File schliessen =======================================================================================
        !===========================================================================================================
        call h5fclose_f(file_id,herror)
    !===========================================================================================================


    end subroutine read_stats_hdf_4D



  
  
  
  
    subroutine write_hdf_infoREAL(NN,scalar_yes,attr_yes,name,array,scalar)
  
        implicit none
  
        integer          , intent(in)  ::  NN
        logical          , intent(in)  ::  scalar_yes
        logical          , intent(in)  ::  attr_yes
        character(*)     , intent(in)  ::  name
        real   , optional, intent(in)  ::  array(1:NN)
        real   , optional, intent(in)  ::  scalar
  
        integer(HID_T)                 ::  memspace
        integer(HID_T)                 ::  attr_id
        integer(HSIZE_T)               ::  dim_mem(1)
        real                           ::  value(1:NN)
  
  
        dim_mem = (/NN/)
  
        if (scalar_yes) then
            value = scalar
        else
            value = array
        end if
  
        call h5screate_simple_f(1,dim_mem,memspace,herror)
  
        if (attr_yes) then
            call h5acreate_f(dset_id,name,memtypeREAL,memspace,attr_id,herror)
            call h5awrite_f(attr_id,memtypeREAL,value,dim_mem,herror)
            call h5aclose_f(attr_id,herror)
        else
            call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,memspace,attr_id,herror)
            call H5dwrite_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror,memspace,H5S_ALL_F,H5P_DEFAULT_F)
            call h5dclose_f(attr_id ,herror)
        end if
  
        call h5sclose_f(memspace,herror)
  
  
    end subroutine write_hdf_infoREAL
  
  
  
  
  
  
  
  
  
  
  
    subroutine write_hdf_infoINT(NN,scalar_yes,attr_yes,name,array,scalar)
  
        implicit none
  
        integer          , intent(in)  ::  NN
        logical          , intent(in)  ::  scalar_yes
        logical          , intent(in)  ::  attr_yes
        character(*)     , intent(in)  ::  name
        integer, optional, intent(in)  ::  array(1:NN)
        integer, optional, intent(in)  ::  scalar
  
        integer(HID_T)                 ::  memspace
        integer(HID_T)                 ::  attr_id
        integer(HSIZE_T)               ::  dim_mem(1)
        integer                        ::  value(1:NN)
  
  
        dim_mem = (/NN/)
  
        if (scalar_yes) then
            value = scalar
        else
            value = array
        end if
  
        call h5screate_simple_f(1,dim_mem,memspace,herror)
  
        if (attr_yes) then
            call h5acreate_f(dset_id,name,memtypeINT ,memspace,attr_id,herror)
            call h5awrite_f(attr_id,memtypeINT ,value,dim_mem,herror)
            call h5aclose_f(attr_id,herror)
        else
            call h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,memspace,attr_id,herror)
            call H5dwrite_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,memspace,H5S_ALL_F,H5P_DEFAULT_F)
            call h5dclose_f(attr_id ,herror)
        end if
  
        call h5sclose_f(memspace,herror)
  
  
    end subroutine write_hdf_infoINT
  
  
  
  
  
  
  
  
  
  
  
    subroutine write_hdf_infoLOG(NN,scalar_yes,attr_yes,name,array,scalar)
  
        implicit none
  
        integer          , intent(in)  ::  NN
        logical          , intent(in)  ::  scalar_yes
        logical          , intent(in)  ::  attr_yes
        character(*)     , intent(in)  ::  name
        logical, optional, intent(in)  ::  array(1:NN)
        logical, optional, intent(in)  ::  scalar
  
        integer(HID_T)                 ::  memspace
        integer(HID_T)                 ::  attr_id
        integer(HSIZE_T)               ::  dim_mem(1)
        integer                        ::  m
        integer                        ::  value(1:NN)
  
  
        dim_mem = (/NN/)
  
        if (scalar_yes) then
            if (scalar) then
                value = 1
            else
                value = 0
            end if
        else
            do m = 1, NN
                if (array(m)) then
                    value(m) = 1
                else
                    value(m) = 0
                end if
            end do
        end if
  
        call h5screate_simple_f(1,dim_mem,memspace,herror)
  
        if (attr_yes) then
            call h5acreate_f(dset_id,name,memtypeINT ,memspace,attr_id,herror)
            call h5awrite_f(attr_id,memtypeINT ,value,dim_mem,herror)
            call h5aclose_f(attr_id,herror)
        else
            call h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,memspace,attr_id,herror)
            call H5dwrite_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,memspace,H5S_ALL_F,H5P_DEFAULT_F)
            call h5dclose_f(attr_id ,herror)
        end if
  
        call h5sclose_f(memspace,herror)
  
  
    end subroutine write_hdf_infoLOG
  
  
  
  
  
  
  
  
  
  
  
    subroutine read_hdf_infoREAL(NN,scalar_yes,attr_yes,name,array,scalar)
  
        implicit none
  
        integer          , intent(in)  ::  NN
        logical          , intent(in)  ::  scalar_yes
        logical          , intent(in)  ::  attr_yes
        character(*)     , intent(in)  ::  name
        real   , optional, intent(out) ::  array(1:NN)
        real   , optional, intent(out) ::  scalar
  
        integer(HID_T)                 ::  attr_id
        integer(HSIZE_T)               ::  dim_mem(1)
        real                           ::  value(1:NN)
  
  
        dim_mem = (/NN/)
  
        if (attr_yes) then
            call h5aopen_name_f(dset_id,name,attr_id,herror)
            if (rank == 0) call h5aread_f(attr_id,memtypeREAL,value,dim_mem,herror)
            call h5aclose_f(attr_id,herror)
        else
            call h5dopen_f(file_id,name,attr_id,herror)
            if (rank == 0) call H5dread_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
            call h5dclose_f(attr_id,herror)
        end if
  
        call MPI_BCAST(value,NN,MPI_REAL8,0,COMM_CART,merror)
  
        if (scalar_yes) then
            scalar = value(1)
        else
            array  = value
        end if
  
  
    end subroutine read_hdf_infoREAL
  
  
  
  
  
  
  
  
  
  
  
    subroutine read_hdf_infoINT(NN,scalar_yes,attr_yes,name,array,scalar)
  
        implicit none
  
        integer          , intent(in)  ::  NN
        logical          , intent(in)  ::  scalar_yes
        logical          , intent(in)  ::  attr_yes
        character(*)     , intent(in)  ::  name
        integer, optional, intent(out) ::  array(1:NN)
        integer, optional, intent(out) ::  scalar
  
        integer(HID_T)                 ::  attr_id
        integer(HSIZE_T)               ::  dim_mem(1)
        integer                        ::  value(1:NN)
  
  
        dim_mem = (/NN/)
  
        if (attr_yes) then
            call h5aopen_name_f(dset_id,name,attr_id,herror)
            if (rank == 0) call h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
            call h5aclose_f(attr_id,herror)
        else
            call h5dopen_f(file_id,name,attr_id,herror)
            if (rank == 0) call H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
            call h5dclose_f(attr_id,herror)
        end if
  
        call MPI_BCAST(value,NN,MPI_INTEGER,0,COMM_CART,merror)
  
        if (scalar_yes) then
            scalar = value(1)
        else
            array  = value
        end if
  
  
    end subroutine read_hdf_infoINT
  
  
  
  
  
  
  
  
  
  
  
    subroutine read_hdf_infoLOG(NN,scalar_yes,attr_yes,name,array,scalar)
  
        implicit none
  
        integer          , intent(in)  ::  NN
        logical          , intent(in)  ::  scalar_yes
        logical          , intent(in)  ::  attr_yes
        character(*)     , intent(in)  ::  name
        logical, optional, intent(out) ::  array(1:NN)
        logical, optional, intent(out) ::  scalar
  
        integer(HID_T)                 ::  attr_id
        integer(HSIZE_T)               ::  dim_mem(1)
        integer                        ::  m
        integer                        ::  value(1:NN)
  
  
        dim_mem = (/NN/)
  
        if (attr_yes) then
            call h5aopen_name_f(dset_id,name,attr_id,herror)
            if (rank == 0) call h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
            call h5aclose_f(attr_id,herror)
        else
            call h5dopen_f(file_id,name,attr_id,herror)
            if (rank == 0) call H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
            call h5dclose_f(attr_id,herror)
        end if
  
        call MPI_BCAST(value,NN,MPI_LOGICAL,0,COMM_CART,merror)
  
        if (scalar_yes) then
            if (value(1) == 1) then
                scalar = .true.
            else
                scalar = .false.
            end if
        else
            do m = 1, NN
                if (value(m) == 1) then
                    array(m) = .true.
                else
                    array(m) = .false.
                end if
            end do
        end if
  
  
    end subroutine read_hdf_infoLOG
  
  
  
  
  
  
  
  
  
  
  
    subroutine filespace_props(vel_dir,stride,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
        implicit none
  
        integer          , intent(in)  ::  vel_dir
        integer          , intent(in)  ::  stride(1:3)
        integer          , intent(out) ::  S1w, M1w, dim1
        integer          , intent(out) ::  S2w, M2w, dim2
        integer          , intent(out) ::  S3w, M3w, dim3
  
  
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
        if (BC_1L_global /= -1) then
            if (vel_dir == 1) then
                if (BC_1L_global == -2) then
                    S1w = 1
                    if (iB(1,1) > 1 .and. ls1 ==  0) offset_file(1) = offset_file(1) + 1
                else
                    S1w = 0
                    if (iB(1,1) > 1 .and. ls1 ==  0) offset_file(1) = offset_file(1) + 2
                    if (iB(1,1) > 1 .and. ls1 == -1) offset_file(1) = offset_file(1) + 1
                end if
                if (BC_1U_global == -2) then
                    M1w = M1-1
                else
                    M1w = M1
                end if
            else
                S1w = 1
                M1w = (M1-1)/stride(1) + 1
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (BC_2L_global /= -1) then
            if (vel_dir == 2) then
                if (BC_2L_global == -2) then
                    S2w = 1
                    if (iB(2,1) > 1 .and. ls2 ==  0) offset_file(2) = offset_file(2) + 1
                else
                    S2w = 0
                    if (iB(2,1) > 1 .and. ls2 ==  0) offset_file(2) = offset_file(2) + 2
                    if (iB(2,1) > 1 .and. ls2 == -1) offset_file(2) = offset_file(2) + 1
                end if
                if (BC_2U_global == -2) then
                    M2w = M2-1
                else
                    M2w = M2
                end if
            else
                S2w = 1
                M2w = (M2-1)/stride(2) + 1
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (BC_3L_global /= -1) then
            if (vel_dir == 3) then
                if (BC_3L_global == -2) then
                    S3w = 1
                    if (iB(3,1) > 1 .and. ls3 ==  0) offset_file(3) = offset_file(3) + 1
                else
                    S3w = 0
                    if (iB(3,1) > 1 .and. ls3 ==  0) offset_file(3) = offset_file(3) + 2
                    if (iB(3,1) > 1 .and. ls3 == -1) offset_file(3) = offset_file(3) + 1
                end if
                if (BC_3U_global == -2) then
                    M3w = M3-1
                else
                    M3w = M3
                end if
            else
                S3w = 1
                M3w = (M3-1)/stride(3) + 1
            end if
        end if
        !===========================================================================================================
  
        dim1 = M1w-S1w+1
        dim2 = M2w-S2w+1
        dim3 = M3w-S3w+1
  
        dims_file = (/dim1,dim2,dim3/)
  
  
    end subroutine filespace_props
  
  
  
  
  
  
  
  
  
  
  
    subroutine lambda
  
        implicit none
  
        real                   ::  Dvel(1:3,1:3)
        real                   ::  TT(1:6)
  
        real                   ::  PP, QQ, RR
        real                   ::  rho, theta
        real                   ::  eps
  
        real                   ::  temp
        real                   ::  lam(1:3)
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
        integer                ::  m, n
        real                   ::  pi
  
  
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
        call exchange(1,1,vel(b1L,b2L,b3L,1))
        call exchange(2,2,vel(b1L,b2L,b3L,2))
        call exchange(3,3,vel(b1L,b2L,b3L,3))
  
        !===========================================================================================================
  
        do k = S3p, N3p
            do j = S2p, N2p
                do i = S1p, N1p
                    nl(i,j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
                    !pgi$ unroll = n:8
                    do ii = d1L+1, d1U
                        nl(i,j,k,1) = nl(i,j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1)
                    end do
                end do
            end do
        end do
  
        !===========================================================================================================
  
        do k = S3p, N3p
            do j = S2p, N2p
                do i = S1p, N1p
                    nl(i,j,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
                    !pgi$ unroll = n:8
                    do jj = d2L+1, d2U
                        nl(i,j,k,2) = nl(i,j,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
                    end do
                end do
            end do
        end do
  
        !===========================================================================================================
  
        do k = S3p, N3p
            do j = S2p, N2p
                do i = S1p, N1p
                    nl(i,j,k,3) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
                    !pgi$ unroll = n:8
                    do kk = d3L+1, d3U
                        nl(i,j,k,3) = nl(i,j,k,3) + cIwp(kk,k)*vel(i,j,k+kk,3)
                    end do
                end do
            end do
        end do
  
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        !=== lambda ================================================================================================
        !===========================================================================================================
  
        call exchange_all_all(.false.,nl)
  
        do k = S3p, N3p
            do j = S2p, N2p
                do i = S1p, N1p
           
                    Dvel = 0.
           
                    !--- d/dx -----------------------------------------------------------------------------------------
                    !pgi$ unroll = n:8
                    do ii = b1L, b1U
                        Dvel(:,1) = Dvel(:,1) + cp1(ii,i)*nl(i+ii,j,k,:)
                    end do
           
                    !--- d/dy -----------------------------------------------------------------------------------------
                    !pgi$ unroll = n:8
                    do jj = b2L, b2U
                        Dvel(:,2) = Dvel(:,2) + cp2(jj,j)*nl(i,j+jj,k,:)
                    end do
           
                    !--- d/dz -----------------------------------------------------------------------------------------
                    !pgi$ unroll = n:8
                    do kk = b3L, b3U
                        Dvel(:,3) = Dvel(:,3) + cp3(kk,k)*nl(i,j,k+kk,:)
                    end do
           
           
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
                    if (rho <= eps) then
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
                    else
                        rho = SQRT(rho)
                        QQ  = ((PP/3.)**3 - QQ*(PP/6.) + RR/2.)/rho**3

                        if (ABS(QQ) < 1.) then
                            theta = ACOS(QQ)/3.
                        else
                            if (QQ > 0.) then
                                theta = 0.
                            else
                                theta = pi/3.
                            end if
                        end if
              
                        lam(1) = COS(theta           )
                        lam(2) = COS(theta + 2.*pi/3.)
                        lam(3) = COS(theta + 4.*pi/3.)
              
              
                        !--- sortieren ---------------------------------------------------------------------------------
                        do m = 1, 3
                            !pgi$ unroll = n:8
                            do n = m+1, 3
                                if (lam(m) > lam(n)) then
                                    temp   = lam(n)
                                    lam(n) = lam(m)
                                    lam(m) = temp
                                end if
                            end do
                        end do
              
                        ! Faktor 1/2, da bei TT(1:6) mit 2 durchmultipliziert wurde ...
                        res(i,j,k) = rho*lam(2) + PP/6.
                    end if
           
                end do
            end do
        end do
  
  
    end subroutine lambda
  
  
  
end module cmod_inout
