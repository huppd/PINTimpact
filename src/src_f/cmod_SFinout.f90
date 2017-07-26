!********************************************************************************
!* IMPACT                                                                       *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.*
!* Mai 2005 - Dec 2011                                                          *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.*
!* Apr 2012                                                                     *
!********************************************************************************

module cmod_SFinout

  use iso_c_binding
  use mpi
  use HDF5

  implicit none

  integer         ::  herror
  integer         ::  merror

contains

  subroutine openH5F(  ) bind (c,name='openH5F')

    implicit none

    call h5open_f(herror)

  end subroutine openH5F


  subroutine closeH5F(  ) bind (c,name='closeH5F')

    implicit none

    call h5close_f(herror)

  end subroutine closeH5F


  subroutine write_hdf_2D(    &
      rank,                   &
      COMM_CART,              &
      M,                      &
      BC_L_global,            &
      BC_U_global,            &
      N,                      &
      bL,bU,                  &
      SS,                     &
      NN,                     &
      ls,                     &
      NB,                     &
      iB,                     &
      iShift,                 &
      ftype,                  &
      filecount,              &
      namelen,                &
      phi,                    &
      y1p,                    &
      y2p,                    &
      Re,                     &
      alpha2                  &
      ) bind (c,name='write_hdf5_2D')

    implicit none

    integer(c_int), intent(in) :: rank

    integer(c_int), intent(in) :: COMM_Cart

    integer(c_int), intent(in) :: M(1:3)

    integer(c_int), intent(in) :: BC_L_global(1:3)
    integer(c_int), intent(in) :: BC_U_global(1:3)

    integer(c_int), intent(in) :: N(1:3)

    integer(c_int), intent(in) :: bL(1:3)
    integer(c_int), intent(in) :: bU(1:3)

    integer(c_int), intent(in) :: SS(1:3)

    integer(c_int), intent(in) :: NN(1:3)

    integer(c_int), intent(in) :: ls(1:3)

    integer(c_int), intent(in) :: NB(1:3)

    integer(c_int), intent(in) :: iB(1:3)

    integer(c_int), intent(in) :: iShift(1:3)

    integer(c_int), intent(in) :: ftype

    integer(c_int), intent(in) :: filecount

    integer(c_int), intent(in) :: namelen

    real(c_double), intent(in) :: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

    real(c_double), intent(in) :: y1p( 1:M(1) )
    real(c_double), intent(in) :: y2p( 1:M(2) )
    !        real(c_double), intent(in) :: y3p( 1:M(3) )

    real(c_double), intent(in) :: Re
    real(c_double), intent(in) :: alpha2


    real(c_double)     ::  temp(1:N(1),1:N(2))

    integer(c_int)     ::  i,j

    integer(c_int)     ::  S1w, M1w, dim1
    integer(c_int)     ::  S2w, M2w, dim2

    character(len=5)   ::  count_char

    character(len=namelen) :: filename
    character(len=namelen) :: dsetname

    logical            ::  attr_yes



    integer(HID_T)  ::  memtypeREAL, memtypeINT

    integer(HID_T)  ::  plist_id, file_id, dset_id

    integer(HID_T)  ::  filespace, memspace
    !        integer(HID_T)                 ::  file_id, plist_id


    integer( HSIZE_T)        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
    integer(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)


    !============================================================================
    !=== Ausschrieb-Nr. als String fuer File-Namen ==============================
    !============================================================================
    !  call num_to_string(5,write_count,count_char)
    call num_to_string(rank,5,filecount,count_char)
    if( 0==ftype ) filename = 'velX_'//count_char
    if( 1==ftype ) filename = 'velY_'//count_char
    if( 4==ftype ) filename = 'pre_'//count_char
    !        filename = 'pre_'//count_char
    !
    if( 0==ftype ) dsetname = 'velX'
    if( 1==ftype ) dsetname = 'velY'
    if( 4==ftype ) dsetname = 'pre'
    !        dsetname = 'pre'
    !============================================================================

    do j = SS(2), NN(2)
      do i = SS(1), NN(1)
        temp(i,j) = phi(i,j,1)
      end do
    end do

    !============================================================================
    !=== Attribut-Datentyp ======================================================
    !============================================================================
    call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
    call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
    !============================================================================


    !============================================================================
    !=== Dimensionen und Offsets ================================================
    !============================================================================
    ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden
    ! wegen vel_dir nicht global eingefuehrt.

    S1w = 2  + ls(1)
    S2w = 2  + ls(2)

    M1w = M(1) + ls(1)
    M2w = M(2) + ls(2)

    if (BC_L_global(1) /= -1) then
      S1w = 1
      M1w = M(1)
    end if
    if (BC_L_global(2) /= -1) then
      S2w = 1
      M2w = M(2)
    end if
    !============================================================================
    dim1 = M1w-S1w+1
    dim2 = M2w-S2w+1

    dims_file   = (/           dim1,          dim2 /)
    dims_mem    = (/           N(1),          N(2) /)
    dims_data   = (/(NN(1)-SS(1)+1),(NN(2)-SS(2)+1)/)
    offset_mem  = (/        SS(1)-1,        SS(2)-1/)
    offset_file = (/      iShift(1),      iShift(2)/)

    if( .not.(( .false. .and. iB(1)==1    ) .or. ( .false. .and. iB(2)==   1 ) .or. ( .true. .and. iB(3)==    1  ) .or.    &
      &     ( .false. .and. iB(1)==NB(1)) .or. ( .false. .and. iB(2)==NB(2)) .or. ( .false. .and. iB(3)== NB(3) ))) then
    dims_data   = (/0,0/)
    offset_mem  = (/0,0/)
    offset_file = (/0,0/)
  end if
  !==============================================================================


  !==============================================================================
  !=== HDF5 schreiben ===========================================================
  !==============================================================================
  ! Setup file access property list with parallel I/O access:
  call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

  ! Create the file collectively:
  call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  call h5pclose_f(plist_id,herror)
  !------------------------------------------------------------------------------
  ! Create file space / memory space:
  call h5screate_simple_f(2,dims_file,filespace,herror)
  call h5screate_simple_f(2,dims_mem ,memspace ,herror)

  ! Create the dataset:
  call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .true.
  !        call write_hdf_infoINT (1     ,.true. ,attr_yes,'restart'          ,scalar=restart          )
  call write_hdf_infoINT (memtypeINT,file_id,dset_id,1     ,.true. ,attr_yes,'write_count'      ,scalar=filecount      )
  !        call write_hdf_infoREAL(memtypeREAL,1     ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  !        call write_hdf_infoREAL(memtypeREAL,1     ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  !        call write_hdf_infoREAL(memtypeREAL,1     ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  !        call write_hdf_infoREAL(memtypeREAL,1     ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  !        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  !        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  !        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )

  call write_hdf_infoINT (memtypeINT,file_id,dset_id,2     ,.false.,attr_yes,'S1w S2w'          ,array =(/S1w,S2w/)      )
  call write_hdf_infoINT (memtypeINT,file_id,dset_id,2     ,.false.,attr_yes,'M1w M2w'          ,array =(/M1w,M2w/)      )

  !        call write_hdf_infoREAL(memtypeREAL,1     ,.true. ,attr_yes,'time'             ,scalar=time             )
  !        call write_hdf_infoREAL(memtypeREAL,1     ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
  !        call write_hdf_infoINT (1     ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )

  !        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'concentration_yes',scalar=concentration_yes)
  !        call write_hdf_infoINT (1     ,.true. ,attr_yes,'n_conc'           ,scalar=n_conc           )
  !        call write_hdf_infoREAL(3     ,.false.,attr_yes,'gravity'          ,array =gravity          )
  call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,1     ,.true. ,attr_yes,'Re'               ,scalar=Re               )
  call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,1     ,.true. ,attr_yes,'alpha2'               ,scalar=alpha2               )
  !        call write_hdf_infoREAL(memtypeREAL,n_conc,.false.,attr_yes,'Sc'               ,array =Sc               )
  !        call write_hdf_infoREAL(memtypeREAL,n_conc,.false.,attr_yes,'Ric'              ,array =Ric              )
  !        call write_hdf_infoREAL(memtypeREAL,n_conc,.false.,attr_yes,'usc'              ,array =usc              )
  !------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)

  ! Create property list for collective dataset write:
  call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)

  ! Write the dataset collectively:
  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,temp,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)

  call h5pclose_f(plist_id ,herror)
  call h5dclose_f(dset_id  ,herror)
  call h5sclose_f(memspace ,herror)
  call h5sclose_f(filespace,herror)
  !------------------------------------------------------------------------------
  call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
  call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim2,.false.,.false.,'VectorY',array=y2p(S2w:M2w))
  !------------------------------------------------------------------------------
  call h5fclose_f(file_id,herror)
  !==============================================================================


  end subroutine write_hdf_2D



  !> \brief write h5 file
  !! 
  !! \param[in] rank rank
  !! \param[in] COMM_CART comm
  !! \param[in] M global grid size
  !! \param[in] BC_L_global global BC
  !! \param[in] BC_U_global global BC
  !! \param[in] N local grid size
  !! \param[in] bl
  !! \param[in] bu
  !! \param[in] SS
  !! \param[in] NN
  !! \param[in] ls
  !! \param[in] NB
  !! \param[in] iB
  !! \param[in] iShift
  !! \param[in] dir_name
  !! \param[in] vel_dir vel_dir should be ftype+1
  !! \param[in] filecount
  !! \param[in] namelen
  !! \param[in] stride
  !! \param[in] phi
  !! \param[in] y1p
  !! \param[in] y2p
  !! \param[in] y3p
  !! \param[in] y1u
  !! \param[in] y2v
  !! \param[in] y3w
  !! \param[in] Re
  !! \param[in] alpha2
  !! \note the interval bounds, offset and block size are not global because of vel_dir
  subroutine write_hdf_3D(    &
      rank,                   &
      COMM_CART,              &
      M,                      &
      BC_L_global,            &
      BC_U_global,            &
      N,                      &
      bL,bU,                  &
      SS,                     &
      NN,                     &
      ls,                     &
      NB,                     &
      iB,                     &
      iShift,                 &
      dir_name,               &
      vel_dir,                &
      filecount,              &
      namelen,                &
      stride,                 &
      phi,                    &
      y1p,                    &
      y2p,                    &
      y3p,                    &
      y1u,                    &
      y2v,                    &
      y3w,                    &
      Re,                     &
      alpha2                  &
      ) bind (c,name='write_hdf_3D')

    implicit none

    integer(c_int), intent(in) :: rank

    integer(c_int), intent(in) :: COMM_Cart

    integer(c_int), intent(in) :: M(1:3)

    integer(c_int), intent(in) :: BC_L_global(1:3)
    integer(c_int), intent(in) :: BC_U_global(1:3)

    integer(c_int), intent(in) :: N(1:3)

    integer(c_int), intent(in) :: bL(1:3)
    integer(c_int), intent(in) :: bU(1:3)

    integer(c_int), intent(in) :: SS(1:3)

    integer(c_int), intent(in) :: NN(1:3)

    integer(c_int), intent(in) :: ls(1:3)

    integer(c_int), intent(in) :: NB(1:3)

    integer(c_int), intent(in) :: iB(1:3)

    integer(c_int), intent(in) :: iShift(1:3)

    integer(c_int), intent(in) :: dir_name

    integer(c_int), intent(in) :: vel_dir

    integer(c_int), intent(in) :: filecount

    integer(c_int), intent(in) :: namelen

    real(c_double), intent(in) :: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

    real(c_double), intent(in) :: y1p( 1:M(1) )
    real(c_double), intent(in) :: y2p( 1:M(2) )
    real(c_double), intent(in) :: y3p( 1:M(3) )

    real(c_double), intent(in) :: y1u( 0:M(1) )
    real(c_double), intent(in) :: y2v( 0:M(2) )
    real(c_double), intent(in) :: y3w( 0:M(3) )

    real(c_double), intent(in) :: Re
    real(c_double), intent(in) :: alpha2


    real(c_double)     ::  temp(1:N(1),1:N(2))

    integer(c_int)     ::  i,j


    character(len=5)   ::  count_char

    character(len=namelen) :: filename
    character(len=namelen) :: dsetname

    logical            ::  attr_yes

    integer(HID_T)  ::  memtypeREAL, memtypeINT

    integer(HID_T)  ::  plist_id, file_id, dset_id

    integer(HID_T)  ::  filespace, memspace
    !        integer(HID_T)                 ::  file_id, plist_id


    integer( HSIZE_T)        ::  dims_file  (1:3)
    integer(HSSIZE_T)        ::  offset_file(1:3)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(c_int), intent(in)    ::  stride(1:3)


    integer(HSIZE_T )      ::  stride_mem(1:3)
    integer(HSIZE_T )      ::  dims_mem  (1:3), dims_data(1:3)
    integer(HSSIZE_T)      ::  offset_mem(1:3)

    integer(c_int)         ::  S1w, M1w, dim1
    integer(c_int)         ::  S2w, M2w, dim2
    integer(c_int)         ::  S3w, M3w, dim3



    !====================================================================================
    !=== Ausschrieb-Nr. als String fuer File-Namen ======================================
    !====================================================================================
    !  call num_to_string(5,write_count,count_char)
    call num_to_string(rank,5,filecount,count_char)
    if( 1==dir_name ) filename = 'velX_'//count_char
    if( 2==dir_name ) filename = 'velY_'//count_char
    if( 3==dir_name ) filename = 'velZ_'//count_char
    if( 5==dir_name ) filename = 'pre_'//count_char
    !        filename = 'pre_'//count_char
    !
    if( 1==dir_name ) dsetname = 'velX'
    if( 2==dir_name ) dsetname = 'velY'
    if( 3==dir_name ) dsetname = 'velZ'
    if( 5==dir_name ) dsetname = 'pre'
    !====================================================================================

    ! Nur fuer Druckgitter! (sonst wird es richtig kompliziert, siehe Multigrid
    ! fuer Helmholtz-Problem!).
    if (vel_dir /= 5) then
      if (stride(vel_dir) /= 1) then
        if (rank == 0) write(*,*) 'ERROR! Cannot write velocity field with stride /= 1 in the corresponding direction!'
        call MPI_FINALIZE(merror)
        stop
      end if
    end if
    if ((MOD((N(1)-1),stride(1)) /= 0) .or. (MOD((N(2)-1),stride(2)) /= 0) .or. (MOD((N(3)-1),stride(3)) /= 0)) then
      if (rank == 0) write(*,*) 'ERROR! Cannot write field with this stride!'
      call MPI_FINALIZE(merror)
      stop
    end if


    !====================================================================================
    !=== Attribut-Datentyp ==============================================================
    !====================================================================================
    call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
    call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
    !====================================================================================


    !====================================================================================
    !=== Dimensionen und Offsets ========================================================
    !====================================================================================
    call filespace_props( &
      vel_dir,            &
      stride,             &
      iShift,             &
      iB,                 &
      ls,                 &
      BC_L_global,        &
      BC_U_global,        &
      M,                  &
      S1w,S2w,S3w,        &
      M1w,M2w,M3w,        &
      dim1,dim2,dim3,     &
      dims_file,          &
      offset_file )

    dims_mem   = (/(N(1)+bU(1)-bL(1)+1),(N(2)+bU(2)-bL(2)+1),(N(3)+bU(3)-bL(3)+1)/)
    dims_data  = (/((NN(1)-SS(1))/stride(1)+1),((NN(2)-SS(2))/stride(2)+1),((NN(3)-SS(3))/stride(3)+1)/)
    offset_mem = (/(SS(1)-bL(1)),(SS(2)-bL(2)),(SS(3)-bL(3))/)
    stride_mem = (/stride(1),stride(2),stride(3)/)
    !====================================================================================


    !====================================================================================
    !=== HDF5 schreiben =================================================================
    !====================================================================================
    ! Setup file access property list with parallel I/O access:
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
    call h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

    ! Create the file collectively:
    call h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
    call h5pclose_f(plist_id,herror)
    !------------------------------------------------------------------------------------
    ! Create file space / memory space:
    call h5screate_simple_f(3,dims_file,filespace,herror)
    call h5screate_simple_f(3,dims_mem ,memspace ,herror)

    ! Create the dataset:
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
    !------------------------------------------------------------------------------------
    ! Write the attributes:
    attr_yes = .true.
    !        call write_hdf_infoINT (1     ,.true. ,attr_yes,'restart'          ,scalar=restart          )
    !        call write_hdf_infoINT (1     ,.true. ,attr_yes,'write_count'      ,scalar=write_count      )
    !        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
    !        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
    !        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
    !        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
    !        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
    !        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
    !        call write_hdf_infoLOG (1     ,.true. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )

    call write_hdf_infoINT (memtypeINT,file_id,dset_id,3     ,.false.,attr_yes,'M1  M2  M3 '      ,array =(/M(1), M(2), M(3) /)  )
    call write_hdf_infoINT (memtypeINT,file_id,dset_id,3     ,.false.,attr_yes,'S1w S2w S3w'      ,array =(/S1w,S2w,S3w/)  )
    call write_hdf_infoINT (memtypeINT,file_id,dset_id,3     ,.false.,attr_yes,'M1w M2w M3w'      ,array =(/M1w,M2w,M3w/)  )
    call write_hdf_infoINT (memtypeINT,file_id,dset_id,3     ,.false.,attr_yes,'NB1 NB2 NB3'      ,array =(/NB(1),NB(2),NB(3)/)  )
    call write_hdf_infoINT (memtypeINT,file_id,dset_id,3     ,.false.,attr_yes,'ls1 ls2 ls3'      ,array =(/ls(1),ls(2),ls(3)/)  )

    !call write_hdf_infoREAL(1     ,.true. ,attr_yes,'time'             ,scalar=time             )
    !call write_hdf_infoREAL(1     ,.true. ,attr_yes,'dtime'            ,scalar=dtime            )
    !call write_hdf_infoINT (1     ,.true. ,attr_yes,'timestep'         ,scalar=timestep         )

    !call write_hdf_infoLOG (1     ,.true. ,attr_yes,'concentration_yes',scalar=concentration_yes)
    !call write_hdf_infoINT (1     ,.true. ,attr_yes,'n_conc'           ,scalar=n_conc           )
    !call write_hdf_infoREAL(3     ,.false.,attr_yes,'gravity'          ,array =gravity          )

    !        call write_hdf_infoREAL(1     ,.true. ,attr_yes,'Re'               ,scalar=Re               )
    call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,1     ,.true. ,attr_yes,'Re'               ,scalar=Re               )


    !call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Sc'               ,array =Sc               )
    !call write_hdf_infoREAL(n_conc,.false.,attr_yes,'Ric'              ,array =Ric              )
    !call write_hdf_infoREAL(n_conc,.false.,attr_yes,'usc'              ,array =usc              )

    !call write_hdf_infoREAL(1     ,.true. ,attr_yes,'CFL   '           ,scalar=CFL              )
    !call write_hdf_infoREAL(1     ,.true. ,attr_yes,'thetaL'           ,scalar=thetaL           )
    !call write_hdf_infoREAL(1     ,.true. ,attr_yes,'epsU  '           ,scalar=epsU             )

    !call write_hdf_infoLOG (1     ,.true. ,attr_yes,'Euler_yes'        ,scalar=Euler_yes        )
    !call write_hdf_infoLOG (1     ,.true. ,attr_yes,'twostep_yes'      ,scalar=twostep_yes      )
    !call write_hdf_infoLOG (1     ,.true. ,attr_yes,'mapping_yes'      ,scalar=mapping_yes      )
    !call write_hdf_infoLOG (1     ,.true. ,attr_yes,'upwind_yes'       ,scalar=upwind_yes       )
    !call write_hdf_infoLOG (1     ,.true. ,attr_yes,'upwind_conc_yes'  ,scalar=upwind_conc_yes  )

    call write_hdf_infoINT (memtypeINT,file_id,dset_id,3,.false.,attr_yes,'BC_iL',array=(/BC_L_global(1),BC_L_global(2),BC_L_global(3)/))
    call write_hdf_infoINT (memtypeINT,file_id,dset_id,3,.false.,attr_yes,'BC_iU',array=(/BC_U_global(1),BC_U_global(2),BC_U_global(3)/))

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

    call write_hdf_infoREAL(memtypeREAL,file_id, dset_id, M(1)      ,.false.,attr_yes,'y1p'  ,array=y1p(1:M(1)))
    call write_hdf_infoREAL(memtypeREAL,file_id, dset_id, M(2)      ,.false.,attr_yes,'y2p'  ,array=y2p(1:M(2)))
    call write_hdf_infoREAL(memtypeREAL,file_id, dset_id, M(3)      ,.false.,attr_yes,'y3p'  ,array=y3p(1:M(3)))

    call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,(M(1)+1)   ,.false.,attr_yes,'y1u'  ,array=y1u(0:M(1)))
    call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,(M(2)+1)   ,.false.,attr_yes,'y2v'  ,array=y2v(0:M(2)))
    call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,(M(3)+1)   ,.false.,attr_yes,'y3w'  ,array=y3w(0:M(3)))
    !------------------------------------------------------------------------------------
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
    !------------------------------------------------------------------------------------
    if (vel_dir == 1) then
      call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim1,.false.,.false.,'VectorX',array=y1u(S1w::stride(1)))
    else
      call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim1,.false.,.false.,'VectorX',array=y1p(S1w::stride(1)))
    end if
    if (vel_dir == 2) then
      call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim2,.false.,.false.,'VectorY',array=y2v(S2w::stride(2)))
    else
      call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim2,.false.,.false.,'VectorY',array=y2p(S2w::stride(2)))
    end if
    if (vel_dir == 3) then
      call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim3,.false.,.false.,'VectorZ',array=y3w(S3w::stride(3)))
    else
      call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim3,.false.,.false.,'VectorZ',array=y3p(S3w::stride(3)))
    end if
    !------------------------------------------------------------------------------------
    call h5fclose_f(file_id,herror)
    !====================================================================================


  end subroutine write_hdf_3D



  subroutine write_hdf_infoREAL(memtypeREAL,file_id,dset_id,NN,scalar_yes,attr_yes,name,array,scalar)

    implicit none

    integer(HID_T)   , intent(in)  ::  memtypeREAL

    integer(HID_T)   , intent(in)  ::  file_id, dset_id

    integer          , intent(in)  ::  NN
    logical          , intent(in)  ::  scalar_yes
    logical          , intent(in)  ::  attr_yes
    character(*)     , intent(in)  ::  name
    real(c_double)   , optional, intent(in)  ::  array(1:NN)
    real(c_double)   , optional, intent(in)  ::  scalar

    integer(HID_T)                 ::  memspace
    integer(HID_T)                 ::  attr_id
    integer(HSIZE_T)               ::  dim_mem(1)
    real(c_double)                 ::  value(1:NN)


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











  subroutine write_hdf_infoINT(memtypeINT,file_id,dset_id,NN,scalar_yes,attr_yes,name,array,scalar)

    implicit none

    integer(HID_T)   , intent(in)  ::  memtypeINT

    integer(HID_T)   , intent(in)  ::  file_id, dset_id

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











  subroutine write_hdf_infoLOG(memtypeINT,file_id,dset_id,NN,scalar_yes,attr_yes,name,array,scalar)

    implicit none

    integer(HID_T)   , intent(in)  ::  memtypeINT

    integer(HID_T)   , intent(in)  ::  file_id, dset_id

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



  subroutine num_to_string(rank,n_digits,num,num_char)

    implicit none

    integer     , intent(in ) ::  rank
    integer     , intent(in ) ::  n_digits
    integer     , intent(in ) ::  num
    character(len=n_digits), intent(out) ::  num_char

    !        integer                   ::  tenthousands, thousands, hundreds, tens, ones


    !====================================================================================
    if (num < 0) then
      if (rank == 0) write(*,*) 'ERROR! Cannot convert negative integers to a string.'
      call MPI_FINALIZE(merror)
      stop
      !==================================================================================
    else if (n_digits >= 6) then
      if (rank == 0) write(*,*) 'ERROR! Cannot convert integers > 99999 to a string.'
      call MPI_FINALIZE(merror)
      stop
      !==================================================================================
    else if (n_digits == 1) then
      write(num_char,'(i1.1)') num
      !==================================================================================
    else if (n_digits == 2) then
      write(num_char,'(i2.2)') num
      !==================================================================================
    else if (n_digits == 3) then
      write(num_char,'(i3.3)') num
      !==================================================================================
    else if (n_digits == 4) then
      write(num_char,'(i4.4)') num
      !==================================================================================
    else if (n_digits == 5) then
      !tenthousands =  num / 10000
      !thousands    = (num - 10000 * tenthousands) / 1000
      !hundreds     = (num - 10000 * tenthousands - 1000 * thousands) / 100
      !tens         = (num - 10000 * tenthousands - 1000 * thousands - 100 * hundreds) / 10
      !ones         =  num - 10000 * tenthousands - 1000 * thousands - 100 * hundreds - tens*10
      !
      !num_char(1:1) = CHAR(ICHAR('0')+tenthousands)
      !num_char(2:2) = CHAR(ICHAR('0')+thousands)
      !num_char(3:3) = CHAR(ICHAR('0')+hundreds)
      !num_char(4:4) = CHAR(ICHAR('0')+tens)
      !num_char(5:5) = CHAR(ICHAR('0')+ones)
      write(num_char,'(i5.5)') num
    end if
    !====================================================================================


  end subroutine num_to_string



  !> \note: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen
  !!          "vel_dir" nicht global eingefuehrt, so dass diese Routine immer
  !!          wieder von Neuem aufgerufen werden muss.                                   
  !!        - Routine dient in erster Linie der kuerzeren Schreibweise /
  !!          besseren Uebersicht                      
  !!        - Default bei Intervallgrenzen sind periodische / Nachbarblock-RB                                     
  subroutine filespace_props( &
      vel_dir,                &
      stride,                 &
      iShift,                 &
      iB,                     &
      ls,                     &
      BC_L_global,            &
      BC_U_global,            &
      M,                      &
      S1w,S2w,S3w,            &
      M1w,M2w,M3w,            &
      dim1,dim2,dim3,         &
      dims_file,              &
      offset_file )

    implicit none

    integer          , intent(in)  ::  vel_dir
    integer          , intent(in)  ::  stride(1:3)
    integer          , intent(in)  ::  iShift(1:3)
    integer          , intent(in)  ::  iB(1:3)
    integer          , intent(in)  ::  ls(1:3)
    integer          , intent(in)  ::  BC_L_global(1:3)
    integer          , intent(in)  ::  BC_U_global(1:3)


    integer          , intent(in)  ::  M(1:3)

    integer          , intent(out) ::  S1w, M1w, dim1
    integer          , intent(out) ::  S2w, M2w, dim2
    integer          , intent(out) ::  S3w, M3w, dim3
    integer( HSIZE_T), intent(out) ::  dims_file(1:3)
    integer(HSSIZE_T), intent(out) ::  offset_file(1:3)


    S1w = 2 + ls(1)
    S2w = 2 + ls(2)
    S3w = 2 + ls(3)

    M1w = (M(1)-1)/stride(1) + 1 + ls(1)
    M2w = (M(2)-1)/stride(2) + 1 + ls(2)
    M3w = (M(3)-1)/stride(3) + 1 + ls(3)

    offset_file = (/iShift(1)/stride(1),iShift(2)/stride(2),iShift(3)/stride(3)/)


    !====================================================================================
    if (BC_L_global(1) /= -1) then
      if (vel_dir == 1) then
        if (BC_L_global(1) == -2) then
          S1w = 1
          if (iB(1) > 1 .and. ls(1) ==  0) offset_file(1) = offset_file(1) + 1
        else
          S1w = 0
          if (iB(1) > 1 .and. ls(1) ==  0) offset_file(1) = offset_file(1) + 2
          if (iB(1) > 1 .and. ls(1) == -1) offset_file(1) = offset_file(1) + 1
        end if
        if (BC_U_global(1) == -2) then
          M1w = M(1)-1
        else
          M1w = M(1)
        end if
      else
        S1w = 1
        M1w = (M(1)-1)/stride(1) + 1
      end if
    end if
    !------------------------------------------------------------------------------------
    if (BC_L_global(2) /= -1) then
      if (vel_dir == 2) then
        if (BC_L_global(2) == -2) then
          S2w = 1
          if (iB(2) > 1 .and. ls(2) ==  0) offset_file(2) = offset_file(2) + 1
        else
          S2w = 0
          if (iB(2) > 1 .and. ls(2) ==  0) offset_file(2) = offset_file(2) + 2
          if (iB(2) > 1 .and. ls(2) == -1) offset_file(2) = offset_file(2) + 1
        end if
        if (BC_U_global(2) == -2) then
          M2w = M(2)-1
        else
          M2w = M(2)
        end if
      else
        S2w = 1
        M2w = (M(2)-1)/stride(2) + 1
      end if
    end if
    !------------------------------------------------------------------------------------
    if (BC_L_global(3) /= -1) then
      if (vel_dir == 3) then
        if (BC_L_global(3) == -2) then
          S3w = 1
          if (iB(3) > 1 .and. ls(3) ==  0) offset_file(3) = offset_file(3) + 1
        else
          S3w = 0
          if (iB(3) > 1 .and. ls(3) ==  0) offset_file(3) = offset_file(3) + 2
          if (iB(3) > 1 .and. ls(3) == -1) offset_file(3) = offset_file(3) + 1
        end if
        if (BC_U_global(3) == -2) then
          M3w = M(3)-1
        else
          M3w = M(3)
        end if
      else
        S3w = 1
        M3w = (M(3)-1)/stride(3) + 1
      end if
    end if
    !====================================================================================

    dim1 = M1w-S1w+1
    dim2 = M2w-S2w+1
    dim3 = M3w-S3w+1

    dims_file = (/dim1,dim2,dim3/)


  end subroutine filespace_props



  !subroutine read_restart(    &
      !rank,                   &
      !COMM_CART,              &
      !MM,                     &
      !BC_L_global,            &
      !BC_U_global,            &
      !N,                      &
      !bL,bU,                  &
      !SS,                     &
      !NN,                     &
      !ls,                     &
      !NB,                     &
      !iB,                     &
      !iShift,                 &
      !dir_name,               &
      !vel_dir,                &
      !filecount,              &
      !namelen,                &
      !stride,                 &
      !phi,                    &
      !M,                      &
      !y1p,                    &
      !y2p,                    &
      !y3p,                    &
      !y1u,                    &
      !y2v,                    &
      !y3w,                    &
      !Re,                     &
      !alpha2                  &
      !) bind (c,name='read_restart')

    !implicit none

    !integer(c_int), intent(in) :: rank

    !integer(c_int), intent(in) :: COMM_Cart

    !integer(c_int), intent(in) :: MM(1:3)

    !integer(c_int), intent(in) :: BC_L_global(1:3)
    !integer(c_int), intent(in) :: BC_U_global(1:3)

    !integer(c_int), intent(in) :: N(1:3)

    !integer(c_int), intent(in) :: bL(1:3)
    !integer(c_int), intent(in) :: bU(1:3)

    !integer(c_int), intent(in) :: SS(1:3)

    !integer(c_int), intent(in) :: NN(1:3)

    !integer(c_int), intent(in) :: ls(1:3)

    !integer(c_int), intent(in) :: NB(1:3)

    !integer(c_int), intent(in) :: iB(1:3)

    !integer(c_int), intent(in) :: iShift(1:3)

    !integer(c_int), intent(in) :: dir_name

    !integer(c_int), intent(in) :: vel_dir

    !integer(c_int), intent(in) :: filecount

    !integer(c_int), intent(in) :: namelen

    !integer(c_int), intent(in)    ::  stride(1:3)

    !real(c_double), intent(inout) :: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

    !integer(c_int), intent(in) :: M(1:3)

    !real(c_double), intent(in) :: y1p( 1:M(1) )
    !real(c_double), intent(in) :: y2p( 1:M(2) )
    !real(c_double), intent(in) :: y3p( 1:M(3) )

    !real(c_double), intent(in) :: y1u( 0:M(1) )
    !real(c_double), intent(in) :: y2v( 0:M(2) )
    !real(c_double), intent(in) :: y3w( 0:M(3) )

    !real(c_double), intent(in) :: Re
    !real(c_double), intent(in) :: alpha2

    !!integer(c_int)         :: m

    !integer(c_int)         :: i, ii
    !integer(c_int)         :: j, jj
    !integer(c_int)         :: k, kk

    !real(c_double)         :: old_dtime_out_scal
    !real(c_double)         :: old_dtime_out_vect

    !character(len=1)       :: conc_number


    !if (rank == 0) write(*,'(a,i3,a)') 'reading data for restart',restart,' ...'
    !if (rank == 0) write(*,*)

    !!------------------------------------------------------------------------------
    !!  call num_to_string(5,write_count,count_char)
    !call num_to_string(rank,5,filecount,count_char)
    !if( 1==dir_name ) filename = 'velX_'//count_char
    !if( 2==dir_name ) filename = 'velY_'//count_char
    !if( 3==dir_name ) filename = 'velZ_'//count_char
    !if( 5==dir_name ) filename = 'pre_'//count_char
    !!        filename = 'pre_'//count_char
    !!
    !if( 1==dir_name ) dsetname = 'velX'
    !if( 2==dir_name ) dsetname = 'velY'
    !if( 3==dir_name ) dsetname = 'velZ'
    !if( 5==dir_name ) dsetname = 'pre'

    !call read2_hdf( filename, dsetname,SS(1),SS(2),SS(3),NN(1),NN(2),NN(3),1,vel(bl(1),bl(2),bl(3),vel_dir))
    !!call read2_hdf('velX_restart'//restart_char,'velX_restart',S11B,S21B,S31B,N11B,N21B,N31B,1,vel(bl(1),bl(2),bl(3),1))
    !!call read2_hdf('velY_restart'//restart_char,'velY_restart',S12B,S22B,S32B,N12B,N22B,N32B,2,vel(bl(1),bl(2),bl(3),2))
    !!if (dimens == 3) call read2_hdf('velZ_restart'//restart_char,'velZ_restart',S13B,S23B,S33B,N13B,N23B,N33B,3,vel(bl(1),bl(2),bl(3),3))
    !!------------------------------------------------------------------------------
    !!call read2_hdf('pre_restart'//restart_char,'pre_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,pre)
    !!------------------------------------------------------------------------------

  !end subroutine read_restart


  !> \note - Kann im Gegensatz zu "read_hdf" das eingelesene Feld beliebg im
  !!         Raum anordnen!
  subroutine read2_hdf( &
      rank,             &
      COMM_CART,        &
      BC_L_global,      &
      BC_U_global,      &
      filename,         &
      dsetname,         &
      N,                &
      bL,bU,            &
      SS,               &
      NN,               &
      ls,               &
      iB,               &
      iShift,           &
      vel_dir,          &
      filecount,        &
      phi )

    implicit none

    integer(c_int), intent(in) :: rank

    integer(c_int), intent(in) :: COMM_Cart

    integer(c_int), intent(in) :: BC_L_global(1:3)
    integer(c_int), intent(in) :: BC_U_global(1:3)

    character(*), intent(in)   :: filename
    character(*), intent(in)   :: dsetname

    integer(c_int), intent(in) :: N(1:3)

    integer(c_int), intent(in) :: bL(1:3)
    integer(c_int), intent(in) :: bU(1:3)

    integer(c_int), intent(in) :: SS(1:3)

    integer(c_int), intent(in) :: NN(1:3)

    integer(c_int), intent(in) :: ls(1:3)

    integer(c_int), intent(in) :: iB(1:3)

    integer(c_int), intent(in) :: iShift(1:3)

    integer(c_int), intent(in) :: vel_dir

    integer(c_int), intent(out):: filecount

    real(c_double), intent(out):: phi(bl(1):(N(1)+bu(1)),bl(2):(N(2)+bu(2)),bl(3):(N(3)+bu(3)))

    integer(HSIZE_T )          :: dims_mem(1:3), dims_data(1:3)
    integer(HSSIZE_T)          :: offset_mem(1:3)

    integer(c_int)             :: S1w, M1w, dim1
    integer(c_int)             :: S2w, M2w, dim2
    integer(c_int)             :: S3w, M3w, dim3

    integer(c_int)             :: S1r, N1r, i0, iGrid
    integer(c_int)             :: S2r, N2r, j0, jGrid
    integer(c_int)             :: S3r, N3r, k0, kGrid

    integer(c_int)             :: Siw(1:3), Miw(1:3)

    logical                    :: attr_yes


    integer(HID_T)  :: memtypeREAL, memtypeINT

    integer(HID_T)  :: plist_id, file_id, dset_id

    integer(HID_T)  :: filespace, memspace  

    integer( HSIZE_T)        ::  dims_file  (1:3)
    integer(HSSIZE_T)        ::  offset_file(1:3)

    !-----------------------------------------------------------------------------------!
    ! Hinweis: phi(1:10,15,12) waere somit                                              !
    !       h5dump -d /dsetname -s "11,14,0" -c "1,1,10" filename.h5                    !
    !-----------------------------------------------------------------------------------!


    !====================================================================================
    !=== File oeffnen ===================================================================
    !====================================================================================
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
    !====================================================================================


    !====================================================================================
    !=== Attribute lesen ================================================================
    !====================================================================================
    call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
    call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
    !------------------------------------------------------------------------------------
    ! Open the dataset:
    call h5dopen_f(file_id,dsetname,dset_id,herror)

    if (herror == -1) then
      if (rank == 0) write(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
      call MPI_FINALIZE(merror)
      stop
    end if
    !------------------------------------------------------------------------------------
    ! Read the attributes:
    attr_yes = .true.
    !call read_hdf_infoREAL( rank, 1,.true. ,attr_yes,'time'          ,scalar=time          )
    !call read_hdf_infoREAL( rank, 1,.true. ,attr_yes,'dtime'         ,scalar=dtime         )
    !call read_hdf_infoREAL( rank, 1,.true. ,attr_yes,'time_out_vect' ,scalar=time_out_vect )
    !call read_hdf_infoREAL( rank, 1,.true. ,attr_yes,'time_out_scal' ,scalar=time_out_scal )
    !call read_hdf_infoINT( rank, COMM_CART, memtypeINT, file_id, dset_id, 1,.true. ,attr_yes,'timestep'      ,scalar=timestep      )
    call read_hdf_infoINT( rank, COMM_CART, memtypeINT, file_id, dset_id, 1,.true. ,attr_yes,'write_count'   ,scalar=filecount   )
    !call read_hdf_infoLOG (1,.true. ,attr_yes,'write_out_vect',scalar=write_out_vect)
    !call read_hdf_infoLOG (1,.true. ,attr_yes,'write_out_scal',scalar=write_out_scal)
    !call read_hdf_infoLOG (1,.true. ,attr_yes,'new_dtime'     ,scalar=new_dtime     )
    call read_hdf_infoINT( rank, COMM_CART, memtypeINT, file_id, dset_id, 3,.false.,attr_yes,'S1w S2w S3w'   ,array =Siw           )
    call read_hdf_infoINT( rank, COMM_CART, memtypeINT, file_id, dset_id, 3,.false.,attr_yes,'M1w M2w M3w'   ,array =Miw           )
    !call read_hdf_infoINT( rank, COMM_CART, memtypeINT, file_id, dset_id, 1,.true. ,attr_yes,'n_conc'        ,scalar=n_conc_old    )
    !------------------------------------------------------------------------------------
    call h5dclose_f(dset_id,herror)
    !====================================================================================


    !====================================================================================
    !=== Dimensionen und Offsets ========================================================
    !====================================================================================

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
    dims_mem  = (/(N(1)+bu(1)-bl(1)+1),(N(2)+bu(2)-bl(2)+1),(N(3)+bu(3)-bl(3)+1)/)

    !------------------------------------------------------------------------------------
    iGrid = 0
    jGrid = 0
    kGrid = 0
    if (vel_dir == 1 .and. ib(1) > 1) then
      if (BC_L_global(1)  == -2 .and. ls(1) ==  0) iGrid = 1
      if (BC_L_global(1) > 0 .and. ls(1) ==  0) iGrid = 2
      if (BC_L_global(1) > 0 .and. ls(1) == -1) iGrid = 1
    end if
    if (vel_dir == 2 .and. ib(2) > 1) then
      if (BC_L_global(2)  == -2 .and. ls(2) ==  0) jGrid = 1
      if (BC_L_global(2) > 0 .and. ls(2) ==  0) jGrid = 2
      if (BC_L_global(2) > 0 .and. ls(2) == -1) jGrid = 1
    end if
    if (vel_dir == 3 .and. ib(3) > 1) then
      if (BC_L_global(3)  == -2 .and. ls(3) ==  0) kGrid = 1
      if (BC_L_global(3) > 0 .and. ls(3) ==  0) kGrid = 2
      if (BC_L_global(3) > 0 .and. ls(3) == -1) kGrid = 1
    end if
    !------------------------------------------------------------------------------------
    if ((S1w+i0) <= (NN(1)+iShift(1)) .and. (M1w+i0) >= (SS(1)+iShift(1))) then
      if ((S1w+i0) >= (SS(1)+iShift(1))) then
        S1r = S1w+i0-iShift(1)
      else
        S1r = SS(1)
      end if
      if ((M1w+i0) <= (NN(1)+iShift(1))) then
        N1r = M1w+i0-iShift(1)
      else
        N1r = NN(1)
      end if

      dims_data  (1) = N1r-S1r+1
      offset_mem (1) = S1r-bl(1)
      offset_file(1) = iShift(1)+iGrid-i0

      if (offset_file(1) < 0   ) offset_file(1) = 0
      if (offset_file(1) > dim1) offset_file(1) = dim1
    else
      dims_data  (1) = 0
      offset_mem (1) = 0
      offset_file(1) = 0
    end if
    !------------------------------------------------------------------------------------
    if ((S2w+j0) <= (NN(2)+iShift(2)) .and. (M2w+j0) >= (SS(2)+iShift(2))) then
      if ((S2w+j0) >= (SS(2)+iShift(2))) then
        S2r = S2w+j0-iShift(2)
      else
        S2r = SS(2)
      end if
      if ((M2w+j0) <= (NN(2)+iShift(2))) then
        N2r = M2w+j0-iShift(2)
      else
        N2r = NN(2)
      end if

      dims_data  (2) = N2r-S2r+1
      offset_mem (2) = S2r-bl(2)
      offset_file(2) = iShift(2)+jGrid-j0

      if (offset_file(2) < 0   ) offset_file(2) = 0
      if (offset_file(2) > dim2) offset_file(2) = dim2
    else
      dims_data  (2) = 0
      offset_mem (2) = 0
      offset_file(2) = 0
    end if
    !------------------------------------------------------------------------------------
    if ((S3w+k0) <= (NN(3)+iShift(3)) .and. (M3w+k0) >= (SS(3)+iShift(3))) then
      if ((S3w+k0) >= (SS(3)+iShift(3))) then
        S3r = S3w+k0-iShift(3)
      else
        S3r = SS(3)
      end if
      if ((M3w+k0) <= (NN(3)+iShift(3))) then
        N3r = M3w+k0-iShift(3)
      else
        N3r = NN(3)
      end if

      dims_data  (3) = N3r-S3r+1
      offset_mem (3) = S3r-bl(3)
      offset_file(3) = iShift(3)+kGrid-k0

      if (offset_file(3) < 0   ) offset_file(3) = 0
      if (offset_file(3) > dim3) offset_file(3) = dim3
    else
      dims_data  (3) = 0
      offset_mem (3) = 0
      offset_file(3) = 0
    end if
    !====================================================================================



    !====================================================================================
    !=== Feld lesen =====================================================================
    !====================================================================================
    phi = 0.
    !------------------------------------------------------------------------------------
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
    !====================================================================================


    !====================================================================================
    !=== close file =====================================================================
    !====================================================================================
    call h5fclose_f(file_id,herror)
    !====================================================================================

  end subroutine read2_hdf



  subroutine read_hdf_infoREAL( &
      rank,                     &
      COMM_CART,                &
      file_id,                  &
      dset_id,                  &
      memtypeREAL,              &
      NN,scalar_yes,attr_yes,name,array,scalar)

    implicit none

    integer(c_int), intent(in)  :: rank 

    integer(c_int), intent(in) :: COMM_Cart

    integer(HID_T), intent(in)  ::  file_id , dset_id

    integer(HID_T), intent(inout):: memtypeREAL

    integer(c_int), intent(in)  :: NN
    logical       , intent(in)  :: scalar_yes
    logical       , intent(in)  :: attr_yes
    character(*)  , intent(in)  :: name
    real(c_double), optional, intent(out) ::  array(1:NN)
    real(c_double), optional, intent(out) ::  scalar

    integer(HID_T)                 ::  attr_id
    integer(HSIZE_T)               ::  dim_mem(1)
    real(c_double)                 ::  value(1:NN)


    dim_mem = (/NN/)

    if( attr_yes )then
      call h5aopen_name_f(dset_id,name,attr_id,herror)
      if (rank == 0) call h5aread_f(attr_id,memtypeREAL,value,dim_mem,herror)
      call h5aclose_f(attr_id,herror)
    else
      call h5dopen_f(file_id,name,attr_id,herror)
      if (rank == 0) call H5dread_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
      call h5dclose_f(attr_id,herror)
    end if

    call MPI_BCAST(value,NN,MPI_REAL8,0,COMM_CART,merror)

    if( scalar_yes ) then
      scalar = value(1)
    else
      array  = value
    end if

  end subroutine read_hdf_infoREAL



  subroutine read_hdf_infoINT(  &
      rank,                     &
      COMM_CART,                &
      memtypeINT,               &
      file_id,                  &
      dset_id,                  &
      NN,scalar_yes,attr_yes,name,array,scalar)

    implicit none

    integer(c_int), intent(in)  :: rank 

    integer(c_int), intent(in) :: COMM_Cart

    integer(HID_T), intent(inout):: memtypeINT

    integer(HID_T), intent(in)  ::  file_id, dset_id

    integer(c_int), intent(in)  :: NN
    logical       , intent(in)  :: scalar_yes
    logical       , intent(in)  :: attr_yes
    character(*)  , intent(in)  :: name
    integer(c_int), optional, intent(out) ::  array(1:NN)
    integer(c_int), optional, intent(out) ::  scalar

    integer(HID_T)                 ::  attr_id
    integer(HSIZE_T)               ::  dim_mem(1)
    integer(c_int)                 ::  value(1:NN)


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


end module cmod_SFinout
