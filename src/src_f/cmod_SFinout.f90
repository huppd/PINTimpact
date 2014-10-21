!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* Apr 2012                                                                                                  *
!*************************************************************************************************************

module cmod_SFinout
  
  
    use HDF5
  
    use mpi

    !    use mod_lib, only: num_to_string

    use iso_c_binding

    private
  
    public SF_write_hdf5_2D
    public  openH5F
    public closeH5F

    integer         ::  herror
    integer         ::  merror


contains

    subroutine openH5F(  ) bind (c,name='openH5F')

        implicit none

!        integer(c_int), intent(in) :: herrorout

        call h5open_f(herror)

    end subroutine openH5F

    subroutine closeH5F(  ) bind (c,name='closeH5F')

        implicit none

!        integer(c_int), intent(in) :: herrorout

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
        !        y3p,                    &
        Re,                     &
        alpha2                  &
        )   bind (c,name='write_hdf5_2D')

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


        !===========================================================================================================
        !=== Ausschrieb-Nr. als String fuer File-Namen =============================================================
        !===========================================================================================================
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
        !===========================================================================================================

        do j = SS(2), NN(2)
            do i = SS(1), NN(1)
                temp(i,j) = phi(i,j,1)
            end do
        end do

        !===========================================================================================================
        !=== Attribut-Datentyp =====================================================================================
        !===========================================================================================================
        call h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
        call h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
        !===========================================================================================================


        !===========================================================================================================
        !=== Dimensionen und Offsets ===============================================================================
        !===========================================================================================================
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
        !                eingefuehrt.                                                                              !
        !----------------------------------------------------------------------------------------------------------!
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
        !===========================================================================================================
        dim1 = M1w-S1w+1
        dim2 = M2w-S2w+1

        dims_file   = (/           dim1,          dim2 /)
        dims_mem    = (/           N(1),          N(2) /)
        dims_data   = (/(NN(1)-SS(1)+1),(NN(2)-SS(2)+1)/)
        offset_mem  = (/        SS(1)-1,        SS(2)-1/)
        offset_file = (/      iShift(1),      iShift(2)/)

        if (.not.(( .false. .and. iB(1)==1    ) .or. ( .false. .and. iB(2)==   1 ) .or. ( .true. .and. iB(3)==    1  ) .or.    &
            &     ( .false. .and. iB(1)==NB(1)) .or. ( .false. .and. iB(2)==NB(2)) .or. ( .false. .and. iB(3)== NB(3) ))) then
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
        !-----------------------------------------------------------------------------------------------------------
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
        !-----------------------------------------------------------------------------------------------------------
        call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim1,.false.,.false.,'VectorX',array=y1p(S1w:M1w))
        call write_hdf_infoREAL(memtypeREAL,file_id, dset_id,dim2,.false.,.false.,'VectorY',array=y2p(S2w:M2w))
        !-----------------------------------------------------------------------------------------------------------
        call h5fclose_f(file_id,herror)
    !===========================================================================================================


    end subroutine write_hdf_2D
  
  
  
    subroutine write_hdf_infoREAL(memtypeREAL,file_id,dset_id,NN,scalar_yes,attr_yes,name,array,scalar)

        implicit none

        integer(HID_T)   , intent(in)  ::  memtypeREAL

        integer(HID_T)   , intent(in)  ::  file_id, dset_id

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

        integer                   ::  tenthousands, thousands, hundreds, tens, ones


        !===========================================================================================================
        if (num < 0) then
            if (rank == 0) write(*,*) 'ERROR! Cannot convert negative integers to a string.'
            call MPI_FINALIZE(merror)
            stop
        !===========================================================================================================
        else if (n_digits >= 6) then
            if (rank == 0) write(*,*) 'ERROR! Cannot convert integers > 99999 to a string.'
            call MPI_FINALIZE(merror)
            stop
        !===========================================================================================================
        else if (n_digits == 1) then
            write(num_char,'(i1.1)') num
        !===========================================================================================================
        else if (n_digits == 2) then
            write(num_char,'(i2.2)') num
        !===========================================================================================================
        else if (n_digits == 3) then
            write(num_char,'(i3.3)') num
        !===========================================================================================================
        else if (n_digits == 4) then
            write(num_char,'(i4.4)') num
        !===========================================================================================================
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
    !===========================================================================================================


    end subroutine num_to_string


  
  
  
end module cmod_SFinout
