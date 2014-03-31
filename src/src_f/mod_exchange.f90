!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module that deals with the MPI compunication?
module mod_exchange
  
  
  use mod_dims
  use mod_vars
  
  
  private
  
  public exchange, exchange2, exchange_all_all, exchange_relax, exchange_part
  
  INCLUDE 'mpif.h'
  
  integer                ::  status(MPI_STATUS_SIZE) ! TEST!!! Warum steht das eigentlich hier lokal? gleiches gilt auch fuer die anderen Module ...
  
  contains
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  subroutine exchange(dir,vel_dir,phi)
  
  implicit none
  
  integer, intent(in)    ::  dir
  integer, intent(in)    ::  vel_dir
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  real                   ::  ghost1LR( 1:b1U,1:N2,1:N3) ! TEST!!! Feld global anlegen und all_all-Variante loeschen? Performance ...
  real                   ::  ghost1LS( 1:b1U,1:N2,1:N3)
  real                   ::  ghost1UR(b1L:-1,1:N2,1:N3)
  real                   ::  ghost1US(b1L:-1,1:N2,1:N3)
  
  real                   ::  ghost2LR(1:N1, 1:b2U,1:N3)
  real                   ::  ghost2LS(1:N1, 1:b2U,1:N3)
  real                   ::  ghost2UR(1:N1,b2L:-1,1:N3)
  real                   ::  ghost2US(1:N1,b2L:-1,1:N3)
  
  real                   ::  ghost3LR(1:N1,1:N2, 1:b3U)
  real                   ::  ghost3LS(1:N1,1:N2, 1:b3U)
  real                   ::  ghost3UR(1:N1,1:N2,b3L:-1)
  real                   ::  ghost3US(1:N1,1:N2,b3L:-1)
  
  integer                ::  length1L, length1U
  integer                ::  length2L, length2U
  integer                ::  length3L, length3U
  
  integer                ::  i, j, k
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
  !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
  !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
  !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
  !                kein Prefetch einbaut.                                                                               !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  if (dir == 1) then
     
     length1L = N2*N3*ABS(b1L)
     length1U = N2*N3*ABS(b1U)
     
     if (BC_1L == 0) call MPI_IRECV(ghost1UR,length1L,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
     if (BC_1U == 0) call MPI_IRECV(ghost1LR,length1U,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
     
     if (BC_1U == 0) then
        do i = b1L, -1
           ghost1US(i,1:N2,1:N3) = phi((N1p+1+i),1:N2,1:N3)
        end do
     end if
     if (BC_1L == 0) then
        do i = 1, b1U
           ghost1LS(i,1:N2,1:N3) = phi((S1p-1+i),1:N2,1:N3)
        end do
     end if
     
     if (BC_1U == 0) call MPI_SEND(ghost1US,length1L,MPI_REAL8,rank1U,1,COMM_CART,merror)
     if (BC_1L == 0) call MPI_SEND(ghost1LS,length1U,MPI_REAL8,rank1L,2,COMM_CART,merror)
     
     if (BC_1L == 0) call MPI_WAIT(req1L,status,merror)
     if (BC_1U == 0) call MPI_WAIT(req1U,status,merror)
     
     if (BC_1L == 0) call pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_1U == 0) call pseudocall(ghost1LR)
     
     if (BC_1L == 0) then
        do i = b1L, -1
           phi((S1p+i),1:N2,1:N3) = ghost1UR(i,1:N2,1:N3)
        end do
     end if
     if (BC_1U == 0) then
        do i = 1, b1U
           phi((N1p+i),1:N2,1:N3) = ghost1LR(i,1:N2,1:N3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (BC_1L == -1) then
        do i = b1L, -1
           phi((S1p+i),1:N2,1:N3) = phi((N1p+1+i),1:N2,1:N3)
        end do
        do i = 1, b1U
           phi((N1p+i),1:N2,1:N3) = phi((S1p-1+i),1:N2,1:N3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_1L > 0) phi( b1L  : -1     ,1:N2,1:N3) = 0.
        if (BC_1U > 0) phi((N1+1):(N1+b1U),1:N2,1:N3) = 0.
        
        if (BC_1L  == -2) phi( b1L  :  0     ,1:N2,1:N3) = 0.
        if (BC_1U  == -2) phi( N1   :(N1+b1U),1:N2,1:N3) = 0.
     else
        if (BC_1L > 0 .or. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3) = 0.
        if (BC_1U > 0 .or. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3) = 0.
     end if
     
  end if
  
  !======================================================================================================================
  
  if (dir == 2) then
     
     length2L = N1*N3*ABS(b2L)
     length2U = N1*N3*ABS(b2U)
     
     if (BC_2L == 0) call MPI_IRECV(ghost2UR,length2L,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
     if (BC_2U == 0) call MPI_IRECV(ghost2LR,length2U,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
     
     if (BC_2U == 0) then
        do j = b2L, -1
           ghost2US(1:N1,j,1:N3) = phi(1:N1,(N2p+1+j),1:N3)
        end do
     end if
     if (BC_2L == 0) then
        do j = 1, b2U
           ghost2LS(1:N1,j,1:N3) = phi(1:N1,(S2p-1+j),1:N3)
        end do
     end if
     
     if (BC_2U == 0) call MPI_SEND(ghost2US,length2L,MPI_REAL8,rank2U,3,COMM_CART,merror)
     if (BC_2L == 0) call MPI_SEND(ghost2LS,length2U,MPI_REAL8,rank2L,4,COMM_CART,merror)
     
     if (BC_2L == 0) call MPI_WAIT(req2L,status,merror)
     if (BC_2U == 0) call MPI_WAIT(req2U,status,merror)
     
     if (BC_2L == 0) call pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_2U == 0) call pseudocall(ghost2LR)
     
     if (BC_2L == 0) then
        do j = b2L, -1
           phi(1:N1,(S2p+j),1:N3) = ghost2UR(1:N1,j,1:N3)
        end do
     end if
     if (BC_2U == 0) then
        do j = 1, b2U
           phi(1:N1,(N2p+j),1:N3) = ghost2LR(1:N1,j,1:N3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (BC_2L == -1) then
        do j = b2L, -1
           phi(1:N1,(S2p+j),1:N3) = phi(1:N1,(N2p+1+j),1:N3)
        end do
        do j = 1, b2U
           phi(1:N1,(N2p+j),1:N3) = phi(1:N1,(S2p-1+j),1:N3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_2L > 0) phi(1:N1, b2L  : -1     ,1:N3) = 0.
        if (BC_2U > 0) phi(1:N1,(N2+1):(N2+b2U),1:N3) = 0.
        
        if (BC_2L  == -2) phi(1:N1, b2L  :  0     ,1:N3) = 0.
        if (BC_2U  == -2) phi(1:N1, N2   :(N2+b2U),1:N3) = 0.
     else
        if (BC_2L > 0 .or. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3) = 0.
        if (BC_2U > 0 .or. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3) = 0.
     end if
     
  end if
  
  !======================================================================================================================
  
  if (dir == 3 .and. dimens == 3) then
     
     length3L = N1*N2*ABS(b3L)
     length3U = N1*N2*ABS(b3U)
     
     if (BC_3L == 0) call MPI_IRECV(ghost3UR,length3L,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
     if (BC_3U == 0) call MPI_IRECV(ghost3LR,length3U,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
     
     if (BC_3U == 0) then
        do k = b3L, -1
           ghost3US(1:N1,1:N2,k) = phi(1:N1,1:N2,(N3p+1+k))
        end do
     end if
     if (BC_3L == 0) then
        do k = 1, b3U
           ghost3LS(1:N1,1:N2,k) = phi(1:N1,1:N2,(S3p-1+k))
        end do
     end if
     
     if (BC_3U == 0) call MPI_SEND(ghost3US,length3L,MPI_REAL8,rank3U,5,COMM_CART,merror)
     if (BC_3L == 0) call MPI_SEND(ghost3LS,length3U,MPI_REAL8,rank3L,6,COMM_CART,merror)
     
     if (BC_3L == 0) call MPI_WAIT(req3L,status,merror)
     if (BC_3U == 0) call MPI_WAIT(req3U,status,merror)
     
     if (BC_3L == 0) call pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_3U == 0) call pseudocall(ghost3LR)
     
     if (BC_3L == 0) then
        do k = b3L, -1
           phi(1:N1,1:N2,(S3p+k)) = ghost3UR(1:N1,1:N2,k)
        end do
     end if
     if (BC_3U == 0) then
        do k = 1, b3U
           phi(1:N1,1:N2,(N3p+k)) = ghost3LR(1:N1,1:N2,k)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (BC_3L == -1) then
        do k = b3L, -1
           phi(1:N1,1:N2,(S3p+k)) = phi(1:N1,1:N2,(N3p+1+k))
        end do
        do k = 1, b3U
           phi(1:N1,1:N2,(N3p+k)) = phi(1:N1,1:N2,(S3p-1+k))
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_3L > 0) phi(1:N1,1:N2, b3L  : -1     ) = 0.
        if (BC_3U > 0) phi(1:N1,1:N2,(N3+1):(N3+b3U)) = 0.
        
        if (BC_3L  == -2) phi(1:N1,1:N2, b3L  :  0     ) = 0.
        if (BC_3U  == -2) phi(1:N1,1:N2, N3   :(N3+b3U)) = 0.
     else
        if (BC_3L > 0 .or. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ) = 0.
        if (BC_3U > 0 .or. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U)) = 0.
     end if
     
  end if
  
  !======================================================================================================================
  
  
  end subroutine exchange
  
  
  
  
  
  
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  subroutine exchange2(dir,vel_dir,SS1,SS2,SS3,NN1,NN2,NN3,phi) ! Anmerkung: Routine exakt identisch zu "exchange", allerdings mit variablen Intervallgrenzen
  
  implicit none
  
  integer, intent(in   ) ::  dir
  integer, intent(in   ) ::  vel_dir
  integer, intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  real                   ::  ghost1LR( 1:b1U,SS2:NN2,SS3:NN3)
  real                   ::  ghost1LS( 1:b1U,SS2:NN2,SS3:NN3)
  real                   ::  ghost1UR(b1L:-1,SS2:NN2,SS3:NN3)
  real                   ::  ghost1US(b1L:-1,SS2:NN2,SS3:NN3)
  
  real                   ::  ghost2LR(SS1:NN1, 1:b2U,SS3:NN3)
  real                   ::  ghost2LS(SS1:NN1, 1:b2U,SS3:NN3)
  real                   ::  ghost2UR(SS1:NN1,b2L:-1,SS3:NN3)
  real                   ::  ghost2US(SS1:NN1,b2L:-1,SS3:NN3)
  
  real                   ::  ghost3LR(SS1:NN1,SS2:NN2, 1:b3U)
  real                   ::  ghost3LS(SS1:NN1,SS2:NN2, 1:b3U)
  real                   ::  ghost3UR(SS1:NN1,SS2:NN2,b3L:-1)
  real                   ::  ghost3US(SS1:NN1,SS2:NN2,b3L:-1)
  
  integer                ::  length1L, length1U
  integer                ::  length2L, length2U
  integer                ::  length3L, length3U
  
  integer                ::  i, j, k
  integer                ::  dummy
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
  !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
  !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
  !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchstr�mt werden kann!.                !
  !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
  !                kein Prefetch einbaut.                                                                               !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  if (dir == 1) then
     
     length1L = (NN2-SS2+1)*(NN3-SS3+1)*ABS(b1L)
     length1U = (NN2-SS2+1)*(NN3-SS3+1)*ABS(b1U)
     
     if (BC_1L == 0) call MPI_IRECV(ghost1UR,length1L,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
     if (BC_1U == 0) call MPI_IRECV(ghost1LR,length1U,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
     
     if (BC_1U == 0) then
        do i = b1L, -1
           ghost1US(i,SS2:NN2,SS3:NN3) = phi((N1p+1+i),SS2:NN2,SS3:NN3)
        end do
     end if
     if (BC_1L == 0) then
        do i = 1, b1U
           ghost1LS(i,SS2:NN2,SS3:NN3) = phi((S1p-1+i),SS2:NN2,SS3:NN3)
        end do
     end if
     
     if (BC_1U == 0) call MPI_SEND(ghost1US,length1L,MPI_REAL8,rank1U,1,COMM_CART,merror)
     if (BC_1L == 0) call MPI_SEND(ghost1LS,length1U,MPI_REAL8,rank1L,2,COMM_CART,merror)
     
     if (BC_1L == 0) call MPI_WAIT(req1L,status,merror)
     if (BC_1U == 0) call MPI_WAIT(req1U,status,merror)
     
     if (BC_1L == 0) call pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_1U == 0) call pseudocall(ghost1LR)
     
     if (BC_1L == 0) then
        do i = b1L, -1
           phi((S1p+i),SS2:NN2,SS3:NN3) = ghost1UR(i,SS2:NN2,SS3:NN3)
        end do
     end if
     if (BC_1U == 0) then
        do i = 1, b1U
           phi((N1p+i),SS2:NN2,SS3:NN3) = ghost1LR(i,SS2:NN2,SS3:NN3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (BC_1L == -1) then
        do i = b1L, -1
           phi((S1p+i),SS2:NN2,SS3:NN3) = phi((N1p+1+i),SS2:NN2,SS3:NN3)
        end do
        do i = 1, b1U
           phi((N1p+i),SS2:NN2,SS3:NN3) = phi((S1p-1+i),SS2:NN2,SS3:NN3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_1L > 0) phi( b1L  : -1     ,SS2:NN2,SS3:NN3) = 0.
        if (BC_1U > 0) phi((N1+1):(N1+b1U),SS2:NN2,SS3:NN3) = 0.
        
        if (BC_1L  == -2) phi( b1L  :  0     ,SS2:NN2,SS3:NN3) = 0.
        if (BC_1U  == -2) phi( N1   :(N1+b1U),SS2:NN2,SS3:NN3) = 0.
     else
        if (BC_1L > 0 .or. BC_1L == -2) phi( b1L  : 0      ,SS2:NN2,SS3:NN3) = 0.
        if (BC_1U > 0 .or. BC_1U == -2) phi((N1+1):(N1+b1U),SS2:NN2,SS3:NN3) = 0.
     end if
     
  end if
  
  !======================================================================================================================
  
  if (dir == 2) then
     
     length2L = (NN1-SS1+1)*(NN3-SS3+1)*ABS(b2L)
     length2U = (NN1-SS1+1)*(NN3-SS3+1)*ABS(b2U)
     
     if (BC_2L == 0) call MPI_IRECV(ghost2UR,length2L,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
     if (BC_2U == 0) call MPI_IRECV(ghost2LR,length2U,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
     
     if (BC_2U == 0) then
        do j = b2L, -1
           ghost2US(SS1:NN1,j,SS3:NN3) = phi(SS1:NN1,(N2p+1+j),SS3:NN3)
        end do
     end if
     if (BC_2L == 0) then
        do j = 1, b2U
           ghost2LS(SS1:NN1,j,SS3:NN3) = phi(SS1:NN1,(S2p-1+j),SS3:NN3)
        end do
     end if
     
     if (BC_2U == 0) call MPI_SEND(ghost2US,length2L,MPI_REAL8,rank2U,3,COMM_CART,merror)
     if (BC_2L == 0) call MPI_SEND(ghost2LS,length2U,MPI_REAL8,rank2L,4,COMM_CART,merror)
     
     if (BC_2L == 0) call MPI_WAIT(req2L,status,merror)
     if (BC_2U == 0) call MPI_WAIT(req2U,status,merror)
     
     if (BC_2L == 0) call pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_2U == 0) call pseudocall(ghost2LR)
     
     if (BC_2L == 0) then
        do j = b2L, -1
           phi(SS1:NN1,(S2p+j),SS3:NN3) = ghost2UR(SS1:NN1,j,SS3:NN3)
        end do
     end if
     if (BC_2U == 0) then
        do j = 1, b2U
           phi(SS1:NN1,(N2p+j),SS3:NN3) = ghost2LR(SS1:NN1,j,SS3:NN3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (BC_2L == -1) then
        do j = b2L, -1
           phi(SS1:NN1,(S2p+j),SS3:NN3) = phi(SS1:NN1,(N2p+1+j),SS3:NN3)
        end do
        do j = 1, b2U
           phi(SS1:NN1,(N2p+j),SS3:NN3) = phi(SS1:NN1,(S2p-1+j),SS3:NN3)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_2L > 0) phi(SS1:NN1, b2L  : -1     ,SS3:NN3) = 0.
        if (BC_2U > 0) phi(SS1:NN1,(N2+1):(N2+b2U),SS3:NN3) = 0.
        
        if (BC_2L  == -2) phi(SS1:NN1, b2L  :  0     ,SS3:NN3) = 0.
        if (BC_2U  == -2) phi(SS1:NN1, N2   :(N2+b2U),SS3:NN3) = 0.
     else
        if (BC_2L > 0 .or. BC_2L == -2) phi(SS1:NN1, b2L  : 0      ,SS3:NN3) = 0.
        if (BC_2U > 0 .or. BC_2U == -2) phi(SS1:NN1,(N2+1):(N2+b2U),SS3:NN3) = 0.
     end if
     
  end if
  
  !======================================================================================================================
  
  if (dir == 3 .and. dimens == 3) then
     
     length3L = (NN1-SS1+1)*(NN2-SS2+1)*ABS(b3L)
     length3U = (NN1-SS1+1)*(NN2-SS2+1)*ABS(b3U)
     
     if (BC_3L == 0) call MPI_IRECV(ghost3UR,length3L,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
     if (BC_3U == 0) call MPI_IRECV(ghost3LR,length3U,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
     
     if (BC_3U == 0) then
        do k = b3L, -1
           ghost3US(SS1:NN1,SS2:NN2,k) = phi(SS1:NN1,SS2:NN2,(N3p+1+k))
        end do
     end if
     if (BC_3L == 0) then
        do k = 1, b3U
           ghost3LS(SS1:NN1,SS2:NN2,k) = phi(SS1:NN1,SS2:NN2,(S3p-1+k))
        end do
     end if
     
     if (BC_3U == 0) call MPI_SEND(ghost3US,length3L,MPI_REAL8,rank3U,5,COMM_CART,merror)
     if (BC_3L == 0) call MPI_SEND(ghost3LS,length3U,MPI_REAL8,rank3L,6,COMM_CART,merror)
     
     if (BC_3L == 0) call MPI_WAIT(req3L,status,merror)
     if (BC_3U == 0) call MPI_WAIT(req3U,status,merror)
     
     if (BC_3L == 0) call pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_3U == 0) call pseudocall(ghost3LR)
     
     if (BC_3L == 0) then
        do k = b3L, -1
           phi(SS1:NN1,SS2:NN2,(S3p+k)) = ghost3UR(SS1:NN1,SS2:NN2,k)
        end do
     end if
     if (BC_3U == 0) then
        do k = 1, b3U
           phi(SS1:NN1,SS2:NN2,(N3p+k)) = ghost3LR(SS1:NN1,SS2:NN2,k)
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (BC_3L == -1) then
        do k = b3L, -1
           phi(SS1:NN1,SS2:NN2,(S3p+k)) = phi(SS1:NN1,SS2:NN2,(N3p+1+k))
        end do
        do k = 1, b3U
           phi(SS1:NN1,SS2:NN2,(N3p+k)) = phi(SS1:NN1,SS2:NN2,(S3p-1+k))
        end do
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_3L > 0) phi(SS1:NN1,SS2:NN2, b3L  : -1     ) = 0.
        if (BC_3U > 0) phi(SS1:NN1,SS2:NN2,(N3+1):(N3+b3U)) = 0.
        
        if (BC_3L  == -2) phi(SS1:NN1,SS2:NN2, b3L  :  0     ) = 0.
        if (BC_3U  == -2) phi(SS1:NN1,SS2:NN2, N3   :(N3+b3U)) = 0.
     else
        if (BC_3L > 0 .or. BC_3L == -2) phi(SS1:NN1,SS2:NN2, b3L  : 0      ) = 0.
        if (BC_3U > 0 .or. BC_3U == -2) phi(SS1:NN1,SS2:NN2,(N3+1):(N3+b3U)) = 0.
     end if
     
  end if
  
  !======================================================================================================================
  
  
  end subroutine exchange2
  
  
  
  
  
  
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  subroutine exchange_all_all(vel_yes,phi)
  
  implicit none
  
  logical, intent(in)    ::  vel_yes
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  real                   ::  ghost1LR( 1:b1U,1:N2,1:N3,1:3)
  real                   ::  ghost1LS( 1:b1U,1:N2,1:N3,1:3)
  real                   ::  ghost1UR(b1L:-1,1:N2,1:N3,1:3)
  real                   ::  ghost1US(b1L:-1,1:N2,1:N3,1:3)
  
  real                   ::  ghost2LR(1:N1, 1:b2U,1:N3,1:3)
  real                   ::  ghost2LS(1:N1, 1:b2U,1:N3,1:3)
  real                   ::  ghost2UR(1:N1,b2L:-1,1:N3,1:3)
  real                   ::  ghost2US(1:N1,b2L:-1,1:N3,1:3)
  
  real                   ::  ghost3LR(1:N1,1:N2, 1:b3U,1:3)
  real                   ::  ghost3LS(1:N1,1:N2, 1:b3U,1:3)
  real                   ::  ghost3UR(1:N1,1:N2,b3L:-1,1:3)
  real                   ::  ghost3US(1:N1,1:N2,b3L:-1,1:3)
  
  integer                ::  length1L, length1U
  integer                ::  length2L, length2U
  integer                ::  length3L, length3U
  
  integer                ::  i, j, k
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
  !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
  !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
  !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
  !                kein Prefetch einbaut.                                                                               !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  length1L = N2*N3*ABS(b1L)*3
  length1U = N2*N3*ABS(b1U)*3
  
  if (BC_1L == 0) call MPI_IRECV(ghost1UR,length1L,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
  if (BC_1U == 0) call MPI_IRECV(ghost1LR,length1U,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
  
  if (BC_1U == 0) then
     do i = b1L, -1
        ghost1US(i,1:N2,1:N3,:) = phi((N1p+1+i),1:N2,1:N3,:) ! TEST!!! Hier koennte man auch 1:dimens schreiben!!!
     end do
  end if
  if (BC_1L == 0) then
     do i = 1, b1U
        ghost1LS(i,1:N2,1:N3,:) = phi((S1p-1+i),1:N2,1:N3,:)
     end do
  end if
  
  if (BC_1U == 0) call MPI_SEND(ghost1US,length1L,MPI_REAL8,rank1U,1,COMM_CART,merror)
  if (BC_1L == 0) call MPI_SEND(ghost1LS,length1U,MPI_REAL8,rank1L,2,COMM_CART,merror)
  
  if (BC_1L == 0) call MPI_WAIT(req1L,status,merror)
  if (BC_1U == 0) call MPI_WAIT(req1U,status,merror)
  
  if (BC_1L == 0) call pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  if (BC_1U == 0) call pseudocall(ghost1LR)
  
  if (BC_1L == 0) then
     do i = b1L, -1
        phi((S1p+i),1:N2,1:N3,:) = ghost1UR(i,1:N2,1:N3,:)
     end do
  end if
  if (BC_1U == 0) then
     do i = 1, b1U
        phi((N1p+i),1:N2,1:N3,:) = ghost1LR(i,1:N2,1:N3,:)
     end do
  end if
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (BC_1L == -1) then
     do i = b1L, -1
        phi((S1p+i),1:N2,1:N3,:) = phi((N1p+1+i),1:N2,1:N3,:)
     end do
     do i = 1, b1U
        phi((N1p+i),1:N2,1:N3,:) = phi((S1p-1+i),1:N2,1:N3,:)
     end do
  end if
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (vel_yes) then
     if (BC_1L > 0) phi( b1L  : -1     ,1:N2,1:N3,1) = 0.
     if (BC_1U > 0) phi((N1+1):(N1+b1U),1:N2,1:N3,1) = 0.
     
     if (BC_1L  == -2) phi( b1L  :  0     ,1:N2,1:N3,1) = 0.
     if (BC_1U  == -2) phi( N1   :(N1+b1U),1:N2,1:N3,1) = 0.
  else
     if (BC_1L > 0 .or. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3,1) = 0.
     if (BC_1U > 0 .or. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3,1) = 0.
  end if
  
  if (BC_1L > 0 .or. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3,2) = 0.
  if (BC_1L > 0 .or. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3,3) = 0.
  if (BC_1U > 0 .or. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3,2) = 0.
  if (BC_1U > 0 .or. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3,3) = 0.
  
  !======================================================================================================================
  
  length2L = N1*N3*ABS(b2L)*3
  length2U = N1*N3*ABS(b2U)*3
  
  if (BC_2L == 0) call MPI_IRECV(ghost2UR,length2L,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
  if (BC_2U == 0) call MPI_IRECV(ghost2LR,length2U,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
  
  if (BC_2U == 0) then
     do j = b2L, -1
        ghost2US(1:N1,j,1:N3,:) = phi(1:N1,(N2p+1+j),1:N3,:)
     end do
  end if
  if (BC_2L == 0) then
     do j = 1, b2U
        ghost2LS(1:N1,j,1:N3,:) = phi(1:N1,(S2p-1+j),1:N3,:)
     end do
  end if
  
  if (BC_2U == 0) call MPI_SEND(ghost2US,length2L,MPI_REAL8,rank2U,3,COMM_CART,merror)
  if (BC_2L == 0) call MPI_SEND(ghost2LS,length2U,MPI_REAL8,rank2L,4,COMM_CART,merror)
  
  if (BC_2L == 0) call MPI_WAIT(req2L,status,merror)
  if (BC_2U == 0) call MPI_WAIT(req2U,status,merror)
  
  if (BC_2L == 0) call pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  if (BC_2U == 0) call pseudocall(ghost2LR)
  
  if (BC_2L == 0) then
     do j = b2L, -1
        phi(1:N1,(S2p+j),1:N3,:) = ghost2UR(1:N1,j,1:N3,:)
     end do
  end if
  if (BC_2U == 0) then
     do j = 1, b2U
        phi(1:N1,(N2p+j),1:N3,:) = ghost2LR(1:N1,j,1:N3,:)
     end do
  end if
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (BC_2L == -1) then
     do j = b2L, -1
        phi(1:N1,(S2p+j),1:N3,:) = phi(1:N1,(N2p+1+j),1:N3,:)
     end do
     do j = 1, b2U
        phi(1:N1,(N2p+j),1:N3,:) = phi(1:N1,(S2p-1+j),1:N3,:)
     end do
  end if
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (vel_yes) then
     if (BC_2L > 0) phi(1:N1, b2L  : -1     ,1:N3,2) = 0.
     if (BC_2U > 0) phi(1:N1,(N2+1):(N2+b2U),1:N3,2) = 0.
     
     if (BC_2L  == -2) phi(1:N1, b2L  :  0     ,1:N3,2) = 0.
     if (BC_2U  == -2) phi(1:N1, N2   :(N2+b2U),1:N3,2) = 0.
  else
     if (BC_2L > 0 .or. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3,2) = 0.
     if (BC_2U > 0 .or. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3,2) = 0.
  end if
  
  if (BC_2L > 0 .or. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3,1) = 0.
  if (BC_2L > 0 .or. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3,3) = 0.
  if (BC_2U > 0 .or. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3,1) = 0.
  if (BC_2U > 0 .or. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3,3) = 0.
  
  !======================================================================================================================
  
  if (dimens == 3) then
  
  length3L = N1*N2*ABS(b3L)*3
  length3U = N1*N2*ABS(b3U)*3
  
  if (BC_3L == 0) call MPI_IRECV(ghost3UR,length3L,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
  if (BC_3U == 0) call MPI_IRECV(ghost3LR,length3U,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
  
  if (BC_3U == 0) then
     do k = b3L, -1
        ghost3US(1:N1,1:N2,k,:) = phi(1:N1,1:N2,(N3p+1+k),:)
     end do
  end if
  if (BC_3L == 0) then
     do k = 1, b3U
        ghost3LS(1:N1,1:N2,k,:) = phi(1:N1,1:N2,(S3p-1+k),:)
     end do
  end if
  
  if (BC_3U == 0) call MPI_SEND(ghost3US,length3L,MPI_REAL8,rank3U,5,COMM_CART,merror)
  if (BC_3L == 0) call MPI_SEND(ghost3LS,length3U,MPI_REAL8,rank3L,6,COMM_CART,merror)
  
  if (BC_3L == 0) call MPI_WAIT(req3L,status,merror)
  if (BC_3U == 0) call MPI_WAIT(req3U,status,merror)
  
  if (BC_3L == 0) call pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  if (BC_3U == 0) call pseudocall(ghost3LR)
  
  if (BC_3L == 0) then
     do k = b3L, -1
        phi(1:N1,1:N2,(S3p+k),:) = ghost3UR(1:N1,1:N2,k,:)
     end do
  end if
  if (BC_3U == 0) then
     do k = 1, b3U
        phi(1:N1,1:N2,(N3p+k),:) = ghost3LR(1:N1,1:N2,k,:)
     end do
  end if
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (BC_3L == -1) then
     do k = b3L, -1
        phi(1:N1,1:N2,(S3p+k),:) = phi(1:N1,1:N2,(N3p+1+k),:)
     end do
     do k = 1, b3U
        phi(1:N1,1:N2,(N3p+k),:) = phi(1:N1,1:N2,(S3p-1+k),:)
     end do
  end if
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (vel_yes) then
     if (BC_3L > 0) phi(1:N1,1:N2, b3L  : -1     ,3) = 0.
     if (BC_3U > 0) phi(1:N1,1:N2,(N3+1):(N3+b3U),3) = 0.
     
     if (BC_3L  == -2) phi(1:N1,1:N2, b3L  :  0     ,3) = 0.
     if (BC_3U  == -2) phi(1:N1,1:N2, N3   :(N3+b3U),3) = 0.
  else
     if (BC_3L > 0 .or. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ,3) = 0.
     if (BC_3U > 0 .or. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U),3) = 0.
  end if
  
  if (BC_3L > 0 .or. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ,1) = 0.
  if (BC_3L > 0 .or. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ,2) = 0.
  if (BC_3U > 0 .or. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U),1) = 0.
  if (BC_3U > 0 .or. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U),2) = 0.
  
  end if
  !======================================================================================================================
  
  
  end subroutine exchange_all_all
  
  
  
  
  
  
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  subroutine exchange_relax(g,fb1,fb2,fb3,vel_dir,mirror_yes,phi) ! TEST!!! aufraeumen und Variablen substituieren ...
  
  implicit none
  
  !*************************************************************************************************
  integer, intent(in   ) ::  g
  real   , intent(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  integer                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  integer                ::  rank1L, rank1U, rank2L, rank2U, rank3L, rank3U
  !*************************************************************************************************
  
  integer, intent(in)    ::  fb1
  integer, intent(in)    ::  fb2
  integer, intent(in)    ::  fb3
  
  integer, intent(in)    ::  vel_dir
  logical, intent(in)    ::  mirror_yes
  
  real                   ::  ghost1LR(1:NN(2,g),1:NN(3,g)) ! TEST!!! Arrays sollten frueher oder spaeter global allociert werden ...
  real                   ::  ghost1LS(1:NN(2,g),1:NN(3,g))
  real                   ::  ghost1UR(1:NN(2,g),1:NN(3,g))
  real                   ::  ghost1US(1:NN(2,g),1:NN(3,g))
  
  real                   ::  ghost2LR(1:NN(1,g),1:NN(3,g))
  real                   ::  ghost2LS(1:NN(1,g),1:NN(3,g))
  real                   ::  ghost2UR(1:NN(1,g),1:NN(3,g))
  real                   ::  ghost2US(1:NN(1,g),1:NN(3,g))
  
  real                   ::  ghost3LR(1:NN(1,g),1:NN(2,g))
  real                   ::  ghost3LS(1:NN(1,g),1:NN(2,g))
  real                   ::  ghost3UR(1:NN(1,g),1:NN(2,g))
  real                   ::  ghost3US(1:NN(1,g),1:NN(2,g))
  
  integer                ::  length1
  integer                ::  length2
  integer                ::  length3
  
  integer                ::  i, N1, S1R, N1R
  integer                ::  j, N2, S2R, N2R
  integer                ::  k, N3, S3R, N3R
  
  
  !*****************************************************************************************
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)
  
  rank1L = ngb(1,1,g)
  rank1U = ngb(2,1,g)
  rank2L = ngb(1,2,g)
  rank2U = ngb(2,2,g)
  rank3L = ngb(1,3,g)
  rank3U = ngb(2,3,g)
  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  !*****************************************************************************************
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Intervallgrenzen sind der Einfachheit halber analog zu exchange gewählt.                             !
  !              - Im Fall von Symmetrie-RB muss IMMER entsprechend oben und/oder unten gespiegelt werden.              !
  !              - Die Spiegelung wird anstelle einer Koeffizienten-Manipulation vorgenommen, da ansonsten auch in den  !
  !                Grob-Gitter-Routinen die Punkte ausserhalb der betrachteten Domain zu Null gesetzt werden müssten,   !
  !                was im Endeffekt (vermutlich) der gleiche Aufwand wäre.                                              !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Die "vel_dir"-Variante wird ausschliesslich von allen Helmholtz-Operationen auf dem feinsten Gitter  !
  !                benutzt, INKL. Interpolation.                                                                        !
  !              - Bei Spiegelung (mirror_yes bzw. vel_dir) wird das Nullsetzen des Randes einer expliziten Behandlung  !
  !                in den aufrufenden Routinen vorgezogen, da die Komplexität und der Programmieraufwand im Vergleich   !
  !                zum Geschwindigkeitsgewinn zu gross erscheinen (im Gegensatz zum Poisson-Problem mit Neumann-RB, bei !
  !                dem am Rand kein Laplace-3d-Stencil sondern nur eine 1d-Ableitung vorkommt (!)).                     !
  !              - Die "mirror_yes"-Variante wird ausschliesslich von allen Helmholtz- und Poisson-Operationen benutzt, !
  !                OHNE Interpolation.                                                                                  !
  !              - In der ersten Richtungen kann der Compiler offenbar nicht vektorisieren und in der zweiten zudem     !
  !                nicht unrollen. Der Grund dafür ist nicht klar.                                                      !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  length1 = N2*N3
  
  if (fb1 >= 0 .and. BC_1L ==  0) call MPI_IRECV(ghost1UR,length1,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
  if (fb1 <= 0 .and. BC_1U ==  0) call MPI_IRECV(ghost1LR,length1,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
  
  if (fb1 >= 0 .and. BC_1U ==  0) ghost1US(1:N2,1:N3) = phi(N1R,1:N2,1:N3)
  if (fb1 <= 0 .and. BC_1L ==  0) ghost1LS(1:N2,1:N3) = phi(S1R,1:N2,1:N3)
  
  if (fb1 >= 0 .and. BC_1U ==  0) call MPI_SEND(ghost1US,length1,MPI_REAL8,rank1U,1,COMM_CART,merror)
  if (fb1 <= 0 .and. BC_1L ==  0) call MPI_SEND(ghost1LS,length1,MPI_REAL8,rank1L,2,COMM_CART,merror)
  
  if (fb1 >= 0 .and. BC_1L ==  0) call MPI_WAIT(req1L,status,merror)
  if (fb1 <= 0 .and. BC_1U ==  0) call MPI_WAIT(req1U,status,merror)
  
  if (fb1 >= 0 .and. BC_1L ==  0) call pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  if (fb1 <= 0 .and. BC_1U ==  0) call pseudocall(ghost1LR)
  
  if (fb1 >= 0 .and. BC_1L ==  0) phi(S1R-1,1:N2,1:N3) = ghost1UR(1:N2,1:N3)
  if (fb1 <= 0 .and. BC_1U ==  0) phi(N1R+1,1:N2,1:N3) = ghost1LR(1:N2,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (fb1 >= 0 .and. BC_1L == -1) phi(S1R-1,1:N2,1:N3) = phi(N1R,1:N2,1:N3)
  if (fb1 <= 0 .and. BC_1L == -1) phi(N1R+1,1:N2,1:N3) = phi(S1R,1:N2,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (vel_dir == 1) then
     if (BC_1L == -2) phi(0   ,1:N2,1:N3) = 0.
     if (BC_1U == -2) phi(N1  ,1:N2,1:N3) = 0.
  else if (mirror_yes) then
     if (BC_1L == -2) phi(0   ,1:N2,1:N3) = 0.
     if (BC_1U == -2) phi(N1+1,1:N2,1:N3) = 0.
  end if
  
  !======================================================================================================================
  
  length2 = N1*N3
  
  if (fb2 >= 0 .and. BC_2L ==  0) call MPI_IRECV(ghost2UR,length2,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
  if (fb2 <= 0 .and. BC_2U ==  0) call MPI_IRECV(ghost2LR,length2,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
  
  if (fb2 >= 0 .and. BC_2U ==  0) ghost2US(1:N1,1:N3) = phi(1:N1,N2R,1:N3)
  if (fb2 <= 0 .and. BC_2L ==  0) ghost2LS(1:N1,1:N3) = phi(1:N1,S2R,1:N3)
  
  if (fb2 >= 0 .and. BC_2U ==  0) call MPI_SEND(ghost2US,length2,MPI_REAL8,rank2U,3,COMM_CART,merror)
  if (fb2 <= 0 .and. BC_2L ==  0) call MPI_SEND(ghost2LS,length2,MPI_REAL8,rank2L,4,COMM_CART,merror)
  
  if (fb2 >= 0 .and. BC_2L ==  0) call MPI_WAIT(req2L,status,merror)
  if (fb2 <= 0 .and. BC_2U ==  0) call MPI_WAIT(req2U,status,merror)
  
  if (fb2 >= 0 .and. BC_2L ==  0) call pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  if (fb2 <= 0 .and. BC_2U ==  0) call pseudocall(ghost2LR)
  
  if (fb2 >= 0 .and. BC_2L ==  0) phi(1:N1,S2R-1,1:N3) = ghost2UR(1:N1,1:N3)
  if (fb2 <= 0 .and. BC_2U ==  0) phi(1:N1,N2R+1,1:N3) = ghost2LR(1:N1,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (fb2 >= 0 .and. BC_2L == -1) phi(1:N1,S2R-1,1:N3) = phi(1:N1,N2R,1:N3)
  if (fb2 <= 0 .and. BC_2L == -1) phi(1:N1,N2R+1,1:N3) = phi(1:N1,S2R,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (vel_dir == 2) then
     if (BC_2L == -2) phi(1:N1,0   ,1:N3) = 0.
     if (BC_2U == -2) phi(1:N1,N2  ,1:N3) = 0.
  else if (mirror_yes) then
     if (BC_2L == -2) phi(1:N1,0   ,1:N3) = 0.
     if (BC_2U == -2) phi(1:N1,N2+1,1:N3) = 0.
  end if
  
  !======================================================================================================================
  if (dimens == 3) then
  
  length3 = N1*N2
  
  if (fb3 >= 0 .and. BC_3L ==  0) call MPI_IRECV(ghost3UR,length3,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
  if (fb3 <= 0 .and. BC_3U ==  0) call MPI_IRECV(ghost3LR,length3,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
  
  if (fb3 >= 0 .and. BC_3U ==  0) ghost3US(1:N1,1:N2) = phi(1:N1,1:N2,N3R)
  if (fb3 <= 0 .and. BC_3L ==  0) ghost3LS(1:N1,1:N2) = phi(1:N1,1:N2,S3R)
  
  if (fb3 >= 0 .and. BC_3U ==  0) call MPI_SEND(ghost3US,length3,MPI_REAL8,rank3U,5,COMM_CART,merror)
  if (fb3 <= 0 .and. BC_3L ==  0) call MPI_SEND(ghost3LS,length3,MPI_REAL8,rank3L,6,COMM_CART,merror)
  
  if (fb3 >= 0 .and. BC_3L ==  0) call MPI_WAIT(req3L,status,merror)
  if (fb3 <= 0 .and. BC_3U ==  0) call MPI_WAIT(req3U,status,merror)
  
  if (fb3 >= 0 .and. BC_3L ==  0) call pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  if (fb3 <= 0 .and. BC_3U ==  0) call pseudocall(ghost3LR)
  
  if (fb3 >= 0 .and. BC_3L ==  0) phi(1:N1,1:N2,S3R-1) = ghost3UR(1:N1,1:N2)
  if (fb3 <= 0 .and. BC_3U ==  0) phi(1:N1,1:N2,N3R+1) = ghost3LR(1:N1,1:N2)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (fb3 >= 0 .and. BC_3L == -1) phi(1:N1,1:N2,S3R-1) = phi(1:N1,1:N2,N3R)
  if (fb3 <= 0 .and. BC_3L == -1) phi(1:N1,1:N2,N3R+1) = phi(1:N1,1:N2,S3R)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  if (vel_dir == 3) then
     if (BC_3L == -2) phi(1:N1,1:N2,0   ) = 0.
     if (BC_3U == -2) phi(1:N1,1:N2,N3  ) = 0.
  else if (mirror_yes) then
     if (BC_3L == -2) phi(1:N1,1:N2,0   ) = 0.
     if (BC_3U == -2) phi(1:N1,1:N2,N3+1) = 0.
  end if
  
  end if
  !======================================================================================================================
  
  
  end subroutine exchange_relax
  
  
  
  
  
  
  
  
  
  
  
  subroutine exchange_part(dir,vel_dir,ss,phi)
  
  implicit none
  
  integer, intent(in   ) ::  dir
  integer, intent(in   ) ::  vel_dir
  integer, intent(in   ) ::  ss
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  real                   ::  ghost1LR(0:(N2+1),0:(N3+1))
  real                   ::  ghost1LS(0:(N2+1),0:(N3+1))
  real                   ::  ghost1UR(0:(N2+1),0:(N3+1))
  real                   ::  ghost1US(0:(N2+1),0:(N3+1))
  
  real                   ::  ghost2LR(0:(N1+1),0:(N3+1))
  real                   ::  ghost2LS(0:(N1+1),0:(N3+1))
  real                   ::  ghost2UR(0:(N1+1),0:(N3+1))
  real                   ::  ghost2US(0:(N1+1),0:(N3+1))
  
  real                   ::  ghost3LR(0:(N1+1),0:(N2+1))
  real                   ::  ghost3LS(0:(N1+1),0:(N2+1))
  real                   ::  ghost3UR(0:(N1+1),0:(N2+1))
  real                   ::  ghost3US(0:(N1+1),0:(N2+1))
  
  integer                ::  length
  integer                ::  i, j, k
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
  !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
  !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
  !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
  !                kein Prefetch einbaut.                                                                               !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  if (dir == 1) then
     
     length = (N2+2)*(N3+2)
     
     if (BC_1L == 0) call MPI_IRECV(ghost1UR,length,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
     if (BC_1U == 0) call MPI_IRECV(ghost1LR,length,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
     
     if (BC_1U == 0) ghost1US(0:(N2+1),0:(N3+1)) = phi(N1p+ss,0:(N2+1),0:(N3+1))
     if (BC_1L == 0) ghost1LS(0:(N2+1),0:(N3+1)) = phi(S1p-ss,0:(N2+1),0:(N3+1))
     
     if (BC_1U == 0) call MPI_SEND(ghost1US,length,MPI_REAL8,rank1U,1,COMM_CART,merror)
     if (BC_1L == 0) call MPI_SEND(ghost1LS,length,MPI_REAL8,rank1L,2,COMM_CART,merror)
     
     if (BC_1L == 0) call MPI_WAIT(req1L,status,merror)
     if (BC_1U == 0) call MPI_WAIT(req1U,status,merror)
     
     if (BC_1L == 0) call pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_1U == 0) call pseudocall(ghost1LR)
     
     if (ss == 1) then
        if (BC_1L == 0) phi(S1p  ,0:(N2+1),0:(N3+1)) = ghost1UR(0:(N2+1),0:(N3+1)) + phi(S1p,0:(N2+1),0:(N3+1))
        if (BC_1U == 0) phi(N1p  ,0:(N2+1),0:(N3+1)) = ghost1LR(0:(N2+1),0:(N3+1)) + phi(N1p,0:(N2+1),0:(N3+1))
     else
        if (BC_1L == 0) phi(S1p-1,0:(N2+1),0:(N3+1)) = ghost1UR(0:(N2+1),0:(N3+1))
        if (BC_1U == 0) phi(N1p+1,0:(N2+1),0:(N3+1)) = ghost1LR(0:(N2+1),0:(N3+1))
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (ss == 1) then
        if (BC_1L == -1) then
           phi(S1p  ,0:(N2+1),0:(N3+1)) = phi(S1p,0:(N2+1),0:(N3+1)) + phi(N1p+1,0:(N2+1),0:(N3+1))
           phi(N1p  ,0:(N2+1),0:(N3+1)) = phi(N1p,0:(N2+1),0:(N3+1)) + phi(S1p-1,0:(N2+1),0:(N3+1))
        end if
     else
        if (BC_1L == -1) then
           phi(S1p-1,0:(N2+1),0:(N3+1)) = phi(N1p,0:(N2+1),0:(N3+1))
           phi(N1p+1,0:(N2+1),0:(N3+1)) = phi(S1p,0:(N2+1),0:(N3+1))
        end if
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_1L > 0) phi(  -1,0:(N2+1),0:(N3+1)) = 0.
        if (BC_1U > 0) phi(N1+1,0:(N2+1),0:(N3+1)) = 0.
        
        if (BC_1L  == -2) phi(   0,0:(N2+1),0:(N3+1)) = -phi(   1,0:(N2+1),0:(N3+1))
        if (BC_1U  == -2) phi(N1  ,0:(N2+1),0:(N3+1)) = -phi(N1-1,0:(N2+1),0:(N3+1))
     else
        if (BC_1L > 0) phi(   0,0:(N2+1),0:(N3+1)) = 0.
        if (BC_1U > 0) phi(N1+1,0:(N2+1),0:(N3+1)) = 0.
        
        !IF (ss == 1) THEN
        !   IF (BC_1L == -2) phi(1 ,0:(N2+1),0:(N3+1)) = 2.*phi(1 ,0:(N2+1),0:(N3+1)) ! TEST!!! feedback forces
        !   IF (BC_1U == -2) phi(N1,0:(N2+1),0:(N3+1)) = 2.*phi(N1,0:(N2+1),0:(N3+1))
        !ELSE
           if (BC_1L == -2) phi(   0,0:(N2+1),0:(N3+1)) = 0. ! TEST!!! Fluid
           if (BC_1U == -2) phi(N1+1,0:(N2+1),0:(N3+1)) = 0.
        !END IF
     end if
     
  end if
  
  !======================================================================================================================
  
  if (dir == 2) then
     
     length = (N1+2)*(N3+2)
     
     if (BC_2L == 0) call MPI_IRECV(ghost2UR,length,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
     if (BC_2U == 0) call MPI_IRECV(ghost2LR,length,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
     
     if (BC_2U == 0) ghost2US(0:(N1+1),0:(N3+1)) = phi(0:(N1+1),N2p+ss,0:(N3+1))
     if (BC_2L == 0) ghost2LS(0:(N1+1),0:(N3+1)) = phi(0:(N1+1),S2p-ss,0:(N3+1))
     
     if (BC_2U == 0) call MPI_SEND(ghost2US,length,MPI_REAL8,rank2U,3,COMM_CART,merror)
     if (BC_2L == 0) call MPI_SEND(ghost2LS,length,MPI_REAL8,rank2L,4,COMM_CART,merror)
     
     if (BC_2L == 0) call MPI_WAIT(req2L,status,merror)
     if (BC_2U == 0) call MPI_WAIT(req2U,status,merror)
     
     if (BC_2L == 0) call pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_2U == 0) call pseudocall(ghost2LR)
     
     if (ss == 1) then
        if (BC_2L == 0) phi(0:(N1+1),S2p  ,0:(N3+1)) = ghost2UR(0:(N1+1),0:(N3+1)) + phi(0:(N1+1),S2p,0:(N3+1))
        if (BC_2U == 0) phi(0:(N1+1),N2p  ,0:(N3+1)) = ghost2LR(0:(N1+1),0:(N3+1)) + phi(0:(N1+1),N2p,0:(N3+1))
     else
        if (BC_2L == 0) phi(0:(N1+1),S2p-1,0:(N3+1)) = ghost2UR(0:(N1+1),0:(N3+1))
        if (BC_2U == 0) phi(0:(N1+1),N2p+1,0:(N3+1)) = ghost2LR(0:(N1+1),0:(N3+1))
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (ss == 1) then
        if (BC_2L == -1) then
           phi(0:(N1+1),S2p  ,0:(N3+1)) = phi(0:(N1+1),S2p,0:(N3+1)) + phi(0:(N1+1),N2p+1,0:(N3+1))
           phi(0:(N1+1),N2p  ,0:(N3+1)) = phi(0:(N1+1),N2p,0:(N3+1)) + phi(0:(N1+1),S2p-1,0:(N3+1))
        end if
     else
        if (BC_2L == -1) then
           phi(0:(N1+1),S2p-1,0:(N3+1)) = phi(0:(N1+1),N2p,0:(N3+1))
           phi(0:(N1+1),N2p+1,0:(N3+1)) = phi(0:(N1+1),S2p,0:(N3+1))
        end if
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_2L > 0) phi(0:(N1+1),  -1,0:(N3+1)) = 0.
        if (BC_2U > 0) phi(0:(N1+1),N2+1,0:(N3+1)) = 0.
        
        if (BC_2L  == -2) phi(0:(N1+1),   0,0:(N3+1)) = -phi(0:(N1+1),   1,0:(N3+1))
        if (BC_2U  == -2) phi(0:(N1+1),N2  ,0:(N3+1)) = -phi(0:(N1+1),N2-1,0:(N3+1))
     else
        if (BC_2L > 0) phi(0:(N1+1),   0,0:(N3+1)) = 0.
        if (BC_2U > 0) phi(0:(N1+1),N2+1,0:(N3+1)) = 0.
        
        !IF (ss == 1) THEN
        !   IF (BC_2L == -2) phi(0:(N1+1),1 ,0:(N3+1)) = 2.*phi(0:(N1+1),1 ,0:(N3+1))
        !   IF (BC_2U == -2) phi(0:(N1+1),N2,0:(N3+1)) = 2.*phi(0:(N1+1),N2,0:(N3+1))
        !ELSE
           if (BC_2L == -2) phi(0:(N1+1),   0,0:(N3+1)) = 0.
           if (BC_2U == -2) phi(0:(N1+1),N2+1,0:(N3+1)) = 0.
        !END IF
     end if
     
  end if
  
  !======================================================================================================================
  
  if (dir == 3 .and. dimens == 3) then
     
     length = (N1+2)*(N2+2)
     
     if (BC_3L == 0) call MPI_IRECV(ghost3UR,length,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
     if (BC_3U == 0) call MPI_IRECV(ghost3LR,length,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
     
     if (BC_3U == 0) ghost3US(0:(N1+1),0:(N2+1)) = phi(0:(N1+1),0:(N2+1),N3p+ss)
     if (BC_3L == 0) ghost3LS(0:(N1+1),0:(N2+1)) = phi(0:(N1+1),0:(N2+1),S3p-ss)
     
     if (BC_3U == 0) call MPI_SEND(ghost3US,length,MPI_REAL8,rank3U,5,COMM_CART,merror)
     if (BC_3L == 0) call MPI_SEND(ghost3LS,length,MPI_REAL8,rank3L,6,COMM_CART,merror)
     
     if (BC_3L == 0) call MPI_WAIT(req3L,status,merror)
     if (BC_3U == 0) call MPI_WAIT(req3U,status,merror)
     
     if (BC_3L == 0) call pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     if (BC_3U == 0) call pseudocall(ghost3LR)
     
     if (ss == 1) then
        if (BC_3L == 0) phi(0:(N1+1),0:(N2+1),S3p  ) = ghost3UR(0:(N1+1),0:(N2+1)) + phi(0:(N1+1),0:(N2+1),S3p)
        if (BC_3U == 0) phi(0:(N1+1),0:(N2+1),N3p  ) = ghost3LR(0:(N1+1),0:(N2+1)) + phi(0:(N1+1),0:(N2+1),N3p)
     else
        if (BC_3L == 0) phi(0:(N1+1),0:(N2+1),S3p-1) = ghost3UR(0:(N1+1),0:(N2+1))
        if (BC_3U == 0) phi(0:(N1+1),0:(N2+1),N3p+1) = ghost3LR(0:(N1+1),0:(N2+1))
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (ss == 1) then
        if (BC_3L == -1) then
           phi(0:(N1+1),0:(N2+1),S3p  ) = phi(0:(N1+1),0:(N2+1),S3p) + phi(0:(N1+1),0:(N2+1),N3p+1)
           phi(0:(N1+1),0:(N2+1),N3p  ) = phi(0:(N1+1),0:(N2+1),N3p) + phi(0:(N1+1),0:(N2+1),S3p-1)
        end if
     else
        if (BC_3L == -1) then
           phi(0:(N1+1),0:(N2+1),S3p-1) = phi(0:(N1+1),0:(N2+1),N3p)
           phi(0:(N1+1),0:(N2+1),N3p+1) = phi(0:(N1+1),0:(N2+1),S3p)
        end if
     end if
     
     !-------------------------------------------------------------------------------------------------------------------
     
     if (vel_dir == dir) then
        if (BC_3L > 0) phi(0:(N1+1),0:(N2+1),  -1) = 0.
        if (BC_3U > 0) phi(0:(N1+1),0:(N2+1),N3+1) = 0.
        
        if (BC_3L  == -2) phi(0:(N1+1),0:(N2+1),   0) = -phi(0:(N1+1),0:(N2+1),   1)
        if (BC_3U  == -2) phi(0:(N1+1),0:(N2+1),N3  ) = -phi(0:(N1+1),0:(N2+1),N3-1)
     else
        if (BC_3L > 0) phi(0:(N1+1),0:(N2+1),   0) = 0.
        if (BC_3U > 0) phi(0:(N1+1),0:(N2+1),N3+1) = 0.
        
        !IF (ss == 1) THEN
        !   IF (BC_3L == -2) phi(0:(N1+1),0:(N2+1),1 ) = 2.*phi(0:(N1+1),0:(N2+1),1 )
        !   IF (BC_3U == -2) phi(0:(N1+1),0:(N2+1),N3) = 2.*phi(0:(N1+1),0:(N2+1),N3)
        !ELSE
           if (BC_3L == -2) phi(0:(N1+1),0:(N2+1),   0) = 0.
           if (BC_3U == -2) phi(0:(N1+1),0:(N2+1),N3+1) = 0.
        !END IF
     end if
     
  end if
  
  !======================================================================================================================
  
  
  end subroutine exchange_part
  
  
  
  
end module mod_exchange
