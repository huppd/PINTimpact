!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module
module cmod_TransferOp

  use iso_c_binding

  implicit none

contains


  !>  \brief transferes a
  subroutine OP_Transfer( &
      N,                  &
      bLI,bUI,            &
      SI,NI,              &
      bLO,bUO,            &
      SO,NO,              &
      phiIN,              &
      phiOUT ) bind (c,name='OP_Transfer')


    implicit none


    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bLI(3)
    integer(c_int), intent(in)  :: bUI(3)

    integer(c_int), intent(in)  :: SI(3)
    integer(c_int), intent(in)  :: NI(3)

    integer(c_int), intent(in)  :: bLO(3)
    integer(c_int), intent(in)  :: bUO(3)

    integer(c_int), intent(in)  :: SO(3)
    integer(c_int), intent(in)  :: NO(3)

    real(c_double), intent(in)  :: phiIN(bLI(1):(N(1)+bUI(1)),bLI(2):(N(2)+bUI(2)),bLI(3):(N(3)+bUI(3)))

    real(c_double), intent(out) :: phiOUT(bLO(1):(N(1)+bUO(1)),bLO(2):(N(2)+bUO(2)),bLO(3):(N(3)+bUO(3)))


    phiOUT( SO(1):NO(1), SO(2):NO(2), SO(3):NO(3) ) = &
      phiIN ( SI(1):NI(1), SI(2):NI(2), SI(3):NI(3) )


  end subroutine OP_Transfer



end module cmod_TransferOp
