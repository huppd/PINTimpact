!> IMPACT
!! \author Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)
!! \date Mai 2005 - Dec 2011


!> \brief module that exchanges ghost layers
module cmod_exchange

  use iso_c_binding
  use mpi

  implicit none

  integer :: status(MPI_STATUS_SIZE) ! TEST!!! Warum steht das eigentlich hier
  ! lokal? gleiches gilt auch fuer die anderen Module ...

contains


  !> \note  Pseudo-CALL fuer cmod_exchange. Darf auf keinen Fall vom Compiler
  !! geinlined werden!
  subroutine pseudocall(phi)

    implicit none

    real(c_double)   , intent(in)    ::  phi(:,:,:)

  end subroutine pseudocall



  !> \brief exchanges ghost layer
  !!
  !! \note: - [Sij,Nij] ist Untermenge von [Sip,Nip]
  !!        - [Sip,Nip] ist Untermenge von [1,Ni]
  !!        - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,
  !!          was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.
  !!        - Bei Spiegelung muss die zur Spiegelungsebene orthogonale
  !!          Geschwindigkeitskomponente auch im Vorzeichen gespiegelt werden,
  !!          da die Spiegelungsebene nicht durchströmt werden kann.
  !!        - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler
  !!          ansonsten keine Vektorisierung bzw. kein Prefetch einbaut.
  !!
  subroutine F_exchange(  &
      dimens,             &
      COMM,               &
      rankL, rankU,       &
      N,                  &
      bL, bU,             &
      BCL, BCU,           &
      Sp, Np,             &
      SS, NN,             &
      dir,                &
      vel_dir,            &
      phi )               &
      bind ( c, name='F_exchange' )

    implicit none

    integer(c_int), intent(in   ) :: dimens

    integer(c_int), intent(in   ) :: COMM

    integer(c_int), intent(in)    :: rankL(3), rankU(3)

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: bL(3), bU(3)

    integer(c_int), intent(in)    :: BCL(3), BCU(3)

    integer(c_int), intent(in   ) :: Sp(3), Np(3)

    integer(c_int), intent(in   ) :: SS(3), NN(3)

    integer(c_int), intent(in   ) :: dir
    integer(c_int), intent(in   ) :: vel_dir


    real(c_double), intent(inout) ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double)                ::  ghost1LR( 1:bU(1),SS(2):NN(2),SS(3):NN(3))
    real(c_double)                ::  ghost1LS( 1:bU(1),SS(2):NN(2),SS(3):NN(3))
    real(c_double)                ::  ghost1UR(bL(1):-1,SS(2):NN(2),SS(3):NN(3))
    real(c_double)                ::  ghost1US(bL(1):-1,SS(2):NN(2),SS(3):NN(3))

    real(c_double)                ::  ghost2LR(SS(1):NN(1), 1:bU(2),SS(3):NN(3))
    real(c_double)                ::  ghost2LS(SS(1):NN(1), 1:bU(2),SS(3):NN(3))
    real(c_double)                ::  ghost2UR(SS(1):NN(1),bL(2):-1,SS(3):NN(3))
    real(c_double)                ::  ghost2US(SS(1):NN(1),bL(2):-1,SS(3):NN(3))

    real(c_double)                ::  ghost3LR(SS(1):NN(1),SS(2):NN(2), 1:bU(3))
    real(c_double)                ::  ghost3LS(SS(1):NN(1),SS(2):NN(2), 1:bU(3))
    real(c_double)                ::  ghost3UR(SS(1):NN(1),SS(2):NN(2),bL(3):-1)
    real(c_double)                ::  ghost3US(SS(1):NN(1),SS(2):NN(2),bL(3):-1)

    integer(c_int)                ::  length1L, length1U
    integer(c_int)                ::  length2L, length2U
    integer(c_int)                ::  length3L, length3U

    integer(c_int)                ::  req1L, req1U
    integer(c_int)                ::  req2L, req2U
    integer(c_int)                ::  req3L, req3U

    integer(c_int)                :: merror

    integer(c_int)                ::  i, j, k

    !====================================================================================
    if (dir == 1) then

      ! derterming length of message
      length1L = (NN(2)-SS(2)+1)*(NN(3)-SS(3)+1)*ABS(bL(1))
      length1U = (NN(2)-SS(2)+1)*(NN(3)-SS(3)+1)*ABS(bU(1))

      ! reciving message
      if (BCL(1) == 0) call MPI_IRECV(ghost1UR,length1L,MPI_REAL8,rankL(1),1,COMM,req1L,merror)
      if (BCU(1) == 0) call MPI_IRECV(ghost1LR,length1U,MPI_REAL8,rankU(1),2,COMM,req1U,merror)

      ! copying into sending array
      if (BCU(1) == 0) then
        do i = bL(1), -1
          ghost1US(i,SS(2):NN(2),SS(3):NN(3)) = phi((Np(1)+1+i),SS(2):NN(2),SS(3):NN(3))
        end do
      end if
      if (BCL(1) == 0) then
        do i = 1, bU(1)
          ghost1LS(i,SS(2):NN(2),SS(3):NN(3)) = phi((Sp(1)-1+i),SS(2):NN(2),SS(3):NN(3))
        end do
      end if

      ! sending array
      if (BCU(1) == 0) call MPI_SEND(ghost1US,length1L,MPI_REAL8,rankU(1),1,COMM,merror)
      if (BCL(1) == 0) call MPI_SEND(ghost1LS,length1U,MPI_REAL8,rankL(1),2,COMM,merror)

      ! wait for ariving message
      if (BCL(1) == 0) call MPI_WAIT(req1L,status,merror)
      if (BCU(1) == 0) call MPI_WAIT(req1U,status,merror)

      if (BCL(1) == 0) call pseudocall(ghost1UR) ! Soll den Compiler daran
      ! hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
      if (BCU(1) == 0) call pseudocall(ghost1LR)

      ! copying recevied message
      if (BCL(1) == 0) then
        do i = bL(1), -1
          phi((Sp(1)+i),SS(2):NN(2),SS(3):NN(3)) = ghost1UR(i,SS(2):NN(2),SS(3):NN(3))
        end do
      end if
      if (BCU(1) == 0) then
        do i = 1, bU(1)
          phi((Np(1)+i),SS(2):NN(2),SS(3):NN(3)) = ghost1LR(i,SS(2):NN(2),SS(3):NN(3))
        end do
      end if

      !----------------------------------------------------------------------------------

      if (BCL(1) == -1) then
        do i = bL(1), -1
          phi((Sp(1)+i),SS(2):NN(2),SS(3):NN(3)) = phi((Np(1)+1+i),SS(2):NN(2),SS(3):NN(3))
        end do
        do i = 1, bU(1)
          phi((Np(1)+i),SS(2):NN(2),SS(3):NN(3)) = phi((Sp(1)-1+i),SS(2):NN(2),SS(3):NN(3))
        end do
      end if

      !----------------------------------------------------------------------------------

      if (vel_dir == dir) then
        if (BCL(1) > 0) phi( bL(1)  : -1     ,SS(2):NN(2),SS(3):NN(3)) = 0.
        if (BCU(1) > 0) phi((N(1)+1):(N(1)+bU(1)),SS(2):NN(2),SS(3):NN(3)) = 0.

        if (BCL(1)  == -2) phi( bL(1)  :  0     ,SS(2):NN(2),SS(3):NN(3)) = 0.
        if (BCU(1)  == -2) phi( N(1)   :(N(1)+bU(1)),SS(2):NN(2),SS(3):NN(3)) = 0.
      else
        if (BCL(1) > 0 .or. BCL(1) == -2) phi( bL(1)  : 0      ,SS(2):NN(2),SS(3):NN(3)) = 0.
        if (BCU(1) > 0 .or. BCU(1) == -2) phi((N(1)+1):(N(1)+bU(1)),SS(2):NN(2),SS(3):NN(3)) = 0.
      end if

    end if

    !====================================================================================

    if (dir == 2) then

      length2L = (NN(1)-SS(1)+1)*(NN(3)-SS(3)+1)*ABS(bL(2))
      length2U = (NN(1)-SS(1)+1)*(NN(3)-SS(3)+1)*ABS(bU(2))

      if (BCL(2) == 0) call MPI_IRECV(ghost2UR,length2L,MPI_REAL8,rankL(2),3,COMM,req2L,merror)
      if (BCU(2) == 0) call MPI_IRECV(ghost2LR,length2U,MPI_REAL8,rankU(2),4,COMM,req2U,merror)

      if (BCU(2) == 0) then
        do j = bL(2), -1
          ghost2US(SS(1):NN(1),j,SS(3):NN(3)) = phi(SS(1):NN(1),(Np(2)+1+j),SS(3):NN(3))
        end do
      end if
      if (BCL(2) == 0) then
        do j = 1, bU(2)
          ghost2LS(SS(1):NN(1),j,SS(3):NN(3)) = phi(SS(1):NN(1),(Sp(2)-1+j),SS(3):NN(3))
        end do
      end if

      if (BCU(2) == 0) call MPI_SEND(ghost2US,length2L,MPI_REAL8,rankU(2),3,COMM,merror)
      if (BCL(2) == 0) call MPI_SEND(ghost2LS,length2U,MPI_REAL8,rankL(2),4,COMM,merror)

      if (BCL(2) == 0) call MPI_WAIT(req2L,status,merror)
      if (BCU(2) == 0) call MPI_WAIT(req2U,status,merror)

      ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
      if (BCL(2) == 0) call pseudocall(ghost2UR)
      if (BCU(2) == 0) call pseudocall(ghost2LR)

      if (BCL(2) == 0) then
        do j = bL(2), -1
          phi(SS(1):NN(1),(Sp(2)+j),SS(3):NN(3)) = ghost2UR(SS(1):NN(1),j,SS(3):NN(3))
        end do
      end if
      if (BCU(2) == 0) then
        do j = 1, bU(2)
          phi(SS(1):NN(1),(Np(2)+j),SS(3):NN(3)) = ghost2LR(SS(1):NN(1),j,SS(3):NN(3))
        end do
      end if

      !----------------------------------------------------------------------------------

      if (BCL(2) == -1) then
        do j = bL(2), -1
          phi(SS(1):NN(1),(Sp(2)+j),SS(3):NN(3)) = phi(SS(1):NN(1),(Np(2)+1+j),SS(3):NN(3))
        end do
        do j = 1, bU(2)
          phi(SS(1):NN(1),(Np(2)+j),SS(3):NN(3)) = phi(SS(1):NN(1),(Sp(2)-1+j),SS(3):NN(3))
        end do
      end if

      !----------------------------------------------------------------------------------

      if (vel_dir == dir) then
        if (BCL(2) > 0) phi(SS(1):NN(1), bL(2)  : -1     ,SS(3):NN(3)) = 0.
        if (BCU(2) > 0) phi(SS(1):NN(1),(N(2)+1):(N(2)+bU(2)),SS(3):NN(3)) = 0.

        if (BCL(2)  == -2) phi(SS(1):NN(1), bL(2)  :  0     ,SS(3):NN(3)) = 0.
        if (BCU(2)  == -2) phi(SS(1):NN(1), N(2)   :(N(2)+bU(2)),SS(3):NN(3)) = 0.
      else
        if (BCL(2) > 0 .or. BCL(2) == -2) phi(SS(1):NN(1), bL(2)  : 0      ,SS(3):NN(3)) = 0.
        if (BCU(2) > 0 .or. BCU(2) == -2) phi(SS(1):NN(1),(N(2)+1):(N(2)+bU(2)),SS(3):NN(3)) = 0.
      end if

    end if

    !====================================================================================

    if (dir == 3 .and. dimens == 3) then

      length3L = (NN(1)-SS(1)+1)*(NN(2)-SS(2)+1)*ABS(bL(3))
      length3U = (NN(1)-SS(1)+1)*(NN(2)-SS(2)+1)*ABS(bU(3))

      if (BCL(3) == 0) call MPI_IRECV(ghost3UR,length3L,MPI_REAL8,rankL(3),5,COMM,req3L,merror)
      if (BCU(3) == 0) call MPI_IRECV(ghost3LR,length3U,MPI_REAL8,rankU(3),6,COMM,req3U,merror)

      if (BCU(3) == 0) then
        do k = bL(3), -1
          ghost3US(SS(1):NN(1),SS(2):NN(2),k) = phi(SS(1):NN(1),SS(2):NN(2),(Np(3)+1+k))
        end do
      end if
      if (BCL(3) == 0) then
        do k = 1, bU(3)
          ghost3LS(SS(1):NN(1),SS(2):NN(2),k) = phi(SS(1):NN(1),SS(2):NN(2),(Sp(3)-1+k))
        end do
      end if

      if (BCU(3) == 0) call MPI_SEND(ghost3US,length3L,MPI_REAL8,rankU(3),5,COMM,merror)
      if (BCL(3) == 0) call MPI_SEND(ghost3LS,length3U,MPI_REAL8,rankL(3),6,COMM,merror)

      if (BCL(3) == 0) call MPI_WAIT(req3L,status,merror)
      if (BCU(3) == 0) call MPI_WAIT(req3U,status,merror)

      ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
      if (BCL(3) == 0) call pseudocall(ghost3UR)
      if (BCU(3) == 0) call pseudocall(ghost3LR)

      if (BCL(3) == 0) then
        do k = bL(3), -1
          phi(SS(1):NN(1),SS(2):NN(2),(Sp(3)+k)) = ghost3UR(SS(1):NN(1),SS(2):NN(2),k)
        end do
      end if
      if (BCU(3) == 0) then
        do k = 1, bU(3)
          phi(SS(1):NN(1),SS(2):NN(2),(Np(3)+k)) = ghost3LR(SS(1):NN(1),SS(2):NN(2),k)
        end do
      end if

      !----------------------------------------------------------------------------------

      if (BCL(3) == -1) then
        do k = bL(3), -1
          phi(SS(1):NN(1),SS(2):NN(2),(Sp(3)+k)) = phi(SS(1):NN(1),SS(2):NN(2),(Np(3)+1+k))
        end do
        do k = 1, bU(3)
          phi(SS(1):NN(1),SS(2):NN(2),(Np(3)+k)) = phi(SS(1):NN(1),SS(2):NN(2),(Sp(3)-1+k))
        end do
      end if

      !----------------------------------------------------------------------------------

      if (vel_dir == dir) then
        if (BCL(3)  >  0) phi(SS(1):NN(1),SS(2):NN(2), bL(3)  : -1         ) = 0.
        if (BCU(3)  >  0) phi(SS(1):NN(1),SS(2):NN(2),(N(3)+1):(N(3)+bU(3))) = 0.

        if (BCL(3) == -2) phi(SS(1):NN(1),SS(2):NN(2), bL(3)  :  0         ) = 0.
        if (BCU(3) == -2) phi(SS(1):NN(1),SS(2):NN(2), N(3)   :(N(3)+bU(3))) = 0.
      else
        if (BCL(3) > 0 .or. BCL(3) == -2) phi(SS(1):NN(1),SS(2):NN(2), bL(3)  : 0          ) = 0.
        if (BCU(3) > 0 .or. BCU(3) == -2) phi(SS(1):NN(1),SS(2):NN(2),(N(3)+1):(N(3)+bU(3))) = 0.
      end if

    end if

    !====================================================================================

  end subroutine F_exchange

end module cmod_exchange
