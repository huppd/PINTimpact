!> \brief module providing functions to initiliaze and apply RestrictionOp
module cmod_RestrictionOp

    use iso_c_binding

contains

    subroutine MG_getCR(    &
        S,                  &
        N,                  &
        BC_L, BC_U,         &
        cR ) bind(c,name='MG_getCR')

        implicit none

        integer               ::  i!, ii!, iimax
!        integer               ::  j, jj, jjmax
!        integer               ::  k, kk, kkmax

        integer(c_int), intent(in)  :: S
        integer(c_int), intent(in)  :: N
        integer(c_int), intent(in)  :: BC_L, BC_U

        real(c_double), intent(out)   :: cR(-1:1,S:N)


        cR = 0.

        !===========================================================================================================
        !=== Restriktion, linienweise, 1d ==========================================================================
        !===========================================================================================================

!        iimax = N

        do i = S, N
            cR(-1,i) = 1./4.
            cR( 0,i) = 2./4.
            cR( 1,i) = 1./4.
        end do

        if( BC_L > 0 ) then
            cR(-1,S) = 0.
            cR( 0,S) = 1.
            cR( 1,S) = 0.

!            cR(-1,2) = 0. ! TEST!!! Sollte evtl. noch ergaenzt werden ...
!            cR( 0,2) = 1.
!            cR( 1,2) = 0.
        end if
        if (BC_L == -2) then
            cR( 1,S) = cR( 1,S) + cR(-1,S)
            cR(-1,S) = 0.
        end if

        if (BC_U > 0) then
            cR(-1,N) = 0.
            cR( 0,N) = 1.
            cR( 1,N) = 0.

!            cR(-1,N-1) = 0. ! TEST!!!
!            cR( 0,N-1) = 1.
!            cR( 1,N-1) = 0.
        end if

        if (BC_U == -2) then
            cR(-1,N) = cR( 1,N) + cR(-1,N)
            cR( 1,N) = 0.
        end if

    end subroutine MG_getCR



    subroutine MG_restrict( &
        dimens,             &
        Nf,                 &
        bLf,bUf,            &
        SSf,NNf,            &
        Nc,                 &
        bLc,bUc,            &
        SSc,NNc,            &
        cR1,cR2,cR3,        &
        phif,               &
        phic ) bind(c,name='MG_restrict')

        implicit none

        integer(c_int), intent(in)     :: dimens

        integer(c_int), intent(in)     :: Nf(3)

        integer(c_int), intent(in)     :: bLf(3)
        integer(c_int), intent(in)     :: bUf(3)

        integer(c_int), intent(in)     :: SSf(3)
        integer(c_int), intent(in)     :: NNf(3)

        integer(c_int), intent(in)     :: Nc(3)

        integer(c_int), intent(in)     :: bLc(3)
        integer(c_int), intent(in)     :: bUc(3)

        integer(c_int), intent(in)     :: SSc(3)
        integer(c_int), intent(in)     :: NNc(3)

        real(c_double),  intent(in)    :: cR1 ( -1:1, SSc(1):NNc(1) )
        real(c_double),  intent(in)    :: cR2 ( -1:1, SSc(2):NNc(2) )
        real(c_double),  intent(in)    :: cR3 ( -1:1, SSc(3):NNc(3) )

        real(c_double),  intent(in)    :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

        real(c_double),  intent(out)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


        integer                ::  i, ii!, di, imax, iimax
        integer                ::  j, jj!, dj, jmax, jjmax
        integer                ::  k, kk!, dk, kmax, kkmax
        integer                         :: dd(3)

                !        integer                ::  sizsg(1:3), offsg(1:3), dispg
                !        real   , allocatable   ::  sendbuf(:,:,:)
                !        real                   ::  recvbuf(1:NN(1,g+1)*NN(2,g+1)*NN(3,g+1))


                !----------------------------------------------------------------------------------------------------------!
                ! Anmerkungen: - für allgemeine di, dj, dk geeignet!                                                       !
                !              - überlappende Schicht in Blöcken wird (der Einfachheit halber) ebenfalls ausgetauscht, ist !
                !                aber im Prinzip redundant (genauer: coarse(S1R:N1R,S2R:N2R,S3R:N3R) = ...).               !
                !              - Motivation für diese kurze Routine ist die Möglichkeit, auch Varianten wie Full-Weighting !
                !                etc. ggf. einzubauen, ansonsten könnte sie auch eingespaart werden.                       !
                !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
                !                nicht gebraucht (erleichtert die Programmierung), so dass eine Initialisierung notwendig  !
                !                ist. Dies wiederum bedingt die INTENT(inout)-Deklaration.                                 !
                !----------------------------------------------------------------------------------------------------------!


                !        imax    = NNf(1)
                !        jmax    = NNf(2)
                !        kmax    = NNf(3)

                !        iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
                !        jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
                !        kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1

                !        di      = (NN(1,g)-1)/(iimax-1)
                !        dj      = (NN(2,g)-1)/(jjmax-1)
                !        dk      = (NN(3,g)-1)/(kkmax-1)

        dd = 1
        do i=1,dimens
            dd(i) = ( NNf(i)-SSf(i) )/( NNc(i)-SSc(i) )
        end do


        !        if (BC(1,1,g) > 0 .and. BC(1,2,g) > 0) fine1(  1      ,  1      ,1:NN(3,g)) = (fine1(2        ,1      ,1:NN(3,g)) + fine1(1      ,2        ,1:NN(3,g)))/2.
        !        if (BC(1,1,g) > 0 .and. BC(2,2,g) > 0) fine1(  1      ,  NN(2,g),1:NN(3,g)) = (fine1(2        ,NN(2,g),1:NN(3,g)) + fine1(1      ,NN(2,g)-1,1:NN(3,g)))/2.
        !        if (BC(2,1,g) > 0 .and. BC(1,2,g) > 0) fine1(  NN(1,g),  1      ,1:NN(3,g)) = (fine1(NN(1,g)-1,1      ,1:NN(3,g)) + fine1(NN(1,g),2        ,1:NN(3,g)))/2.
        !        if (BC(2,1,g) > 0 .and. BC(2,2,g) > 0) fine1(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine1(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine1(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
        !
        !        if (BC(1,1,g) > 0 .and. BC(1,3,g) > 0) fine1(  1      ,1:NN(2,g),  1      ) = (fine1(2        ,1:NN(2,g),1      ) + fine1(1      ,1:NN(2,g),2        ))/2.
        !        if (BC(1,1,g) > 0 .and. BC(2,3,g) > 0) fine1(  1      ,1:NN(2,g),  NN(3,g)) = (fine1(2        ,1:NN(2,g),NN(3,g)) + fine1(1      ,1:NN(2,g),NN(3,g)-1))/2.
        !        if (BC(2,1,g) > 0 .and. BC(1,3,g) > 0) fine1(  NN(1,g),1:NN(2,g),  1      ) = (fine1(NN(1,g)-1,1:NN(2,g),1      ) + fine1(NN(1,g),1:NN(2,g),2        ))/2.
        !        if (BC(2,1,g) > 0 .and. BC(2,3,g) > 0) fine1(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine1(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine1(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.
        !
        !        if (BC(1,2,g) > 0 .and. BC(1,3,g) > 0) fine1(1:NN(1,g),  1      ,  1      ) = (fine1(1:NN(1,g),2        ,1      ) + fine1(1:NN(1,g),1      ,2        ))/2.
        !        if (BC(1,2,g) > 0 .and. BC(2,3,g) > 0) fine1(1:NN(1,g),  1      ,  NN(3,g)) = (fine1(1:NN(1,g),2        ,NN(3,g)) + fine1(1:NN(1,g),1      ,NN(3,g)-1))/2.
        !        if (BC(2,2,g) > 0 .and. BC(1,3,g) > 0) fine1(1:NN(1,g),  NN(2,g),  1      ) = (fine1(1:NN(1,g),NN(2,g)-1,1      ) + fine1(1:NN(1,g),NN(2,g),2        ))/2.
        !        if (BC(2,2,g) > 0 .and. BC(2,3,g) > 0) fine1(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine1(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine1(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
        !
        !
        !        if (BC(1,1,g) > 0 .and. BC(1,2,g) > 0) fine2(  1      ,  1      ,1:NN(3,g)) = (fine2(2        ,1      ,1:NN(3,g)) + fine2(1      ,2        ,1:NN(3,g)))/2.
        !        if (BC(1,1,g) > 0 .and. BC(2,2,g) > 0) fine2(  1      ,  NN(2,g),1:NN(3,g)) = (fine2(2        ,NN(2,g),1:NN(3,g)) + fine2(1      ,NN(2,g)-1,1:NN(3,g)))/2.
        !        if (BC(2,1,g) > 0 .and. BC(1,2,g) > 0) fine2(  NN(1,g),  1      ,1:NN(3,g)) = (fine2(NN(1,g)-1,1      ,1:NN(3,g)) + fine2(NN(1,g),2        ,1:NN(3,g)))/2.
        !        if (BC(2,1,g) > 0 .and. BC(2,2,g) > 0) fine2(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine2(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine2(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
        !
        !        if (BC(1,1,g) > 0 .and. BC(1,3,g) > 0) fine2(  1      ,1:NN(2,g),  1      ) = (fine2(2        ,1:NN(2,g),1      ) + fine2(1      ,1:NN(2,g),2        ))/2.
        !        if (BC(1,1,g) > 0 .and. BC(2,3,g) > 0) fine2(  1      ,1:NN(2,g),  NN(3,g)) = (fine2(2        ,1:NN(2,g),NN(3,g)) + fine2(1      ,1:NN(2,g),NN(3,g)-1))/2.
        !        if (BC(2,1,g) > 0 .and. BC(1,3,g) > 0) fine2(  NN(1,g),1:NN(2,g),  1      ) = (fine2(NN(1,g)-1,1:NN(2,g),1      ) + fine2(NN(1,g),1:NN(2,g),2        ))/2.
        !        if (BC(2,1,g) > 0 .and. BC(2,3,g) > 0) fine2(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine2(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine2(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.
        !
        !        if (BC(1,2,g) > 0 .and. BC(1,3,g) > 0) fine2(1:NN(1,g),  1      ,  1      ) = (fine2(1:NN(1,g),2        ,1      ) + fine2(1:NN(1,g),1      ,2        ))/2.
        !        if (BC(1,2,g) > 0 .and. BC(2,3,g) > 0) fine2(1:NN(1,g),  1      ,  NN(3,g)) = (fine2(1:NN(1,g),2        ,NN(3,g)) + fine2(1:NN(1,g),1      ,NN(3,g)-1))/2.
        !        if (BC(2,2,g) > 0 .and. BC(1,3,g) > 0) fine2(1:NN(1,g),  NN(2,g),  1      ) = (fine2(1:NN(1,g),NN(2,g)-1,1      ) + fine2(1:NN(1,g),NN(2,g),2        ))/2.
        !        if (BC(2,2,g) > 0 .and. BC(2,3,g) > 0) fine2(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine2(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine2(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
        !

        !        if (ls1 ==  0 .and. (BC(1,1,g) == 0 .or. BC(1,1,g) == -1)) iimax = iimax-1
        !        if (ls1 == -1 .and. (BC(2,1,g) == 0 .or. BC(2,1,g) == -1)) iimax = iimax-1
        !
        !        if (ls2 ==  0 .and. (BC(1,2,g) == 0 .or. BC(1,2,g) == -1)) jjmax = jjmax-1
        !        if (ls2 == -1 .and. (BC(2,2,g) == 0 .or. BC(2,2,g) == -1)) jjmax = jjmax-1
        !
        !        if (ls3 ==  0 .and. (BC(1,3,g) == 0 .or. BC(1,3,g) == -1)) kkmax = kkmax-1
        !        if (ls3 == -1 .and. (BC(2,3,g) == 0 .or. BC(2,3,g) == -1)) kkmax = kkmax-1


        !        if (1 == 2) then ! TEST!!!
        !            if (add_yes) then
        !                !pgi$ unroll = n:8
        !                coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk) - fine2(1:imax:di,1:jmax:dj,1:kmax:dk)
        !            else
        !                coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk)
        !            end if
        !
        !        else

        !            if (add_yes) then ! TEST!!! etwas seriöser einbauen ...
        !
        !                call exchange_relax(g,0,0,0,0,.true.,fine1)
        !                call exchange_relax(g,0,0,0,0,.true.,fine2)
        !
        !                if (dimens == 3) then
        !                    do kk = 1, kkmax
        !                        k = dk*(kk-1)+1
        !                        do jj = 1, jjmax
        !                            j = dj*(jj-1)+1
        !                            do ii = 1, iimax
        !                                i = di*(ii-1)+1
        !                                coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
        !                                    &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
        !                                    &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
        !                                    &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
        !                                    &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k)) +  &
        !                                    &              cR3(-1,kk,g+1)*(fine1(i,j,k-1)-fine2(i,j,k-1)) +  &
        !                                    &              cR3( 1,kk,g+1)*(fine1(i,j,k+1)-fine2(i,j,k+1))) / 3.
        !                            end do
        !                        end do
        !                    end do
        !                else
        !                    k  = 1
        !                    kk = 1
        !                    do jj = 1, jjmax
        !                        j = dj*(jj-1)+1
        !                        do ii = 1, iimax
        !                            i = di*(ii-1)+1
        !                            coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
        !                                &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
        !                                &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
        !                                &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
        !                                &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k))) / 2.
        !                        end do
        !                    end do
        !                end if
        !
        !            else
        !
        !                call exchange_relax(g,0,0,0,0,.true.,fine1)

        if (dimens == 3) then
            do kk = SSc(3), NNc(3)
                k = dd(3)*(kk-1)+1
                do jj = SSc(2), NNc(2)
                    j = dd(2)*(jj-1)+1
                    do ii = SSc(1), NNc(1)
                        i = dd(1)*(ii-1)+1
                        phic(ii,jj,kk) = ((cR1( 0,ii)+cR2( 0,jj)+cR3( 0,kk))*phif(i  ,j  ,k  ) +  &
                            &              cR1(-1,ii)                       *phif(i-1,j  ,k  ) +  &
                            &              cR1( 1,ii)                       *phif(i+1,j  ,k  ) +  &
                            &                         cR2(-1,jj)            *phif(i  ,j-1,k  ) +  &
                            &                         cR2( 1,jj)            *phif(i  ,j+1,k  ) +  &
                            &                                     cR3(-1,kk)*phif(i  ,j  ,k-1) +  &
                            &                                     cR3( 1,kk)*phif(i  ,j  ,k+1)) / 3.
                    end do
                end do
            end do
        else
            k  = 1
            kk = 1
            do jj = SSc(2), NNc(2)
                j = dd(2)*(jj-1)+1
                do ii = SSc(1), NNc(1)
                    i = dd(1)*(ii-1)+1
                    phic(ii,jj,kk) = ((cR1( 0,ii)+cR2(0,jj))*phif(i  ,j  ,k) +  &
                        &              cR1(-1,ii)           *phif(i-1,j  ,k) +  &
                        &              cR1( 1,ii)           *phif(i+1,j  ,k) +  &
                        &                         cR2(-1,jj)*phif(i  ,j-1,k) +  &
                        &                         cR2( 1,jj)*phif(i  ,j+1,k)) / 2.
                end do
            end do
        end if


    !        end if



    !        if (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) > 1) then
    !
    !            allocate(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
    !
    !            do kk = 1, kkmax
    !                do jj = 1, jjmax
    !                    do ii = 1, iimax
    !                        sendbuf(ii,jj,kk) = coarse(ii,jj,kk)
    !                    end do
    !                end do
    !            end do
    !
    !            call MPI_GATHERv(sendbuf,iimax*jjmax*kkmax,MPI_REAL8,recvbuf,recvR(1,g+1),dispR(1,g+1),MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
    !
    !            deallocate(sendbuf)
    !
    !
    !            if (participate_yes(g+1)) then
    !                do k = 1, n_gather(3,g+1)
    !                    do j = 1, n_gather(2,g+1)
    !                        do i = 1, n_gather(1,g+1)
    !
    !                            sizsg(1:3) = sizsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
    !                            offsg(1:3) = offsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
    !                            dispg      = dispR(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
    !
    !                            do kk = 1, sizsg(3)
    !                                do jj = 1, sizsg(2)
    !                                    do ii = 1, sizsg(1)
    !                                        coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3)) = recvbuf(dispg+ii+(jj-1)*sizsg(1)+(kk-1)*sizsg(1)*sizsg(2))
    !                                    end do
    !                                end do
    !                            end do
    !
    !                        end do
    !                    end do
    !                end do
    !            end if
    !
    !        end if
    !

    end subroutine MG_restrict


end module cmod_RestrictionOp
