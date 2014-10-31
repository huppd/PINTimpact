!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* Apr 2012                                                                                                  *
!*************************************************************************************************************

!> \obsolete
module cmod_write_info
  
    use HDF5
  
    use iso_c_binding


    private

  
    public write_hdf_infoReal
    public write_hdf_infoINT
    public write_hdf_infoLOG
  
    integer                ::  herror
  
    integer(HID_T)         ::  file_id, dset_id!, plist_id
    integer(HID_T)         ::  filespace, memspace
    integer(HID_T)         ::  memtypeREAL, memtypeINT
    integer(HSIZE_T )      ::  dims_file  (1:3) ! wird global eingefuehrt, um Probleme bei Uebergabe als Argument nach "CALL filespace_props" zu vermeiden.
    integer(HSSIZE_T)      ::  offset_file(1:3)

    integer                ::  i, j, k ! TEST!!! wo brauchts das global?
  
  
contains

  
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
  
  
  
  
  
  
  
  
  
  
  
!    subroutine read_hdf_infoREAL(NN,scalar_yes,attr_yes,name,array,scalar)
!
!        implicit none
!
!        integer          , intent(in)  ::  NN
!        logical          , intent(in)  ::  scalar_yes
!        logical          , intent(in)  ::  attr_yes
!        character(*)     , intent(in)  ::  name
!        real   , optional, intent(out) ::  array(1:NN)
!        real   , optional, intent(out) ::  scalar
!
!        integer(HID_T)                 ::  attr_id
!        integer(HSIZE_T)               ::  dim_mem(1)
!        real                           ::  value(1:NN)
!
!
!        dim_mem = (/NN/)
!
!        if (attr_yes) then
!            call h5aopen_name_f(dset_id,name,attr_id,herror)
!            if (rank == 0) call h5aread_f(attr_id,memtypeREAL,value,dim_mem,herror)
!            call h5aclose_f(attr_id,herror)
!        else
!            call h5dopen_f(file_id,name,attr_id,herror)
!            if (rank == 0) call H5dread_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
!            call h5dclose_f(attr_id,herror)
!        end if
!
!        call MPI_BCAST(value,NN,MPI_REAL8,0,COMM_CART,merror)
!
!        if (scalar_yes) then
!            scalar = value(1)
!        else
!            array  = value
!        end if
!
!
!    end subroutine read_hdf_infoREAL
!
!
!
!
!
!
!
!
!
!
!
!    subroutine read_hdf_infoINT(NN,scalar_yes,attr_yes,name,array,scalar)
!
!        implicit none
!
!        integer          , intent(in)  ::  NN
!        logical          , intent(in)  ::  scalar_yes
!        logical          , intent(in)  ::  attr_yes
!        character(*)     , intent(in)  ::  name
!        integer, optional, intent(out) ::  array(1:NN)
!        integer, optional, intent(out) ::  scalar
!
!        integer(HID_T)                 ::  attr_id
!        integer(HSIZE_T)               ::  dim_mem(1)
!        integer                        ::  value(1:NN)
!
!
!        dim_mem = (/NN/)
!
!        if (attr_yes) then
!            call h5aopen_name_f(dset_id,name,attr_id,herror)
!            if (rank == 0) call h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
!            call h5aclose_f(attr_id,herror)
!        else
!            call h5dopen_f(file_id,name,attr_id,herror)
!            if (rank == 0) call H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
!            call h5dclose_f(attr_id,herror)
!        end if
!
!        call MPI_BCAST(value,NN,MPI_INTEGER,0,COMM_CART,merror)
!
!        if (scalar_yes) then
!            scalar = value(1)
!        else
!            array  = value
!        end if
!
!
!    end subroutine read_hdf_infoINT
!
!
!
!
!
!
!
!
!
!
!
!    subroutine read_hdf_infoLOG(NN,scalar_yes,attr_yes,name,array,scalar)
!
!        implicit none
!
!        integer          , intent(in)  ::  NN
!        logical          , intent(in)  ::  scalar_yes
!        logical          , intent(in)  ::  attr_yes
!        character(*)     , intent(in)  ::  name
!        logical, optional, intent(out) ::  array(1:NN)
!        logical, optional, intent(out) ::  scalar
!
!        integer(HID_T)                 ::  attr_id
!        integer(HSIZE_T)               ::  dim_mem(1)
!        integer                        ::  m
!        integer                        ::  value(1:NN)
!
!
!        dim_mem = (/NN/)
!
!        if (attr_yes) then
!            call h5aopen_name_f(dset_id,name,attr_id,herror)
!            if (rank == 0) call h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
!            call h5aclose_f(attr_id,herror)
!        else
!            call h5dopen_f(file_id,name,attr_id,herror)
!            if (rank == 0) call H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
!            call h5dclose_f(attr_id,herror)
!        end if
!
!        call MPI_BCAST(value,NN,MPI_LOGICAL,0,COMM_CART,merror)
!
!        if (scalar_yes) then
!            if (value(1) == 1) then
!                scalar = .true.
!            else
!                scalar = .false.
!            end if
!        else
!            do m = 1, NN
!                if (value(m) == 1) then
!                    array(m) = .true.
!                else
!                    array(m) = .false.
!                end if
!            end do
!        end if
!
!
!    end subroutine read_hdf_infoLOG
  
  
  
  
  
  
  
end module cmod_write_info
