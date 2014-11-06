!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module
module cmod_DivGrad2ndOOp
  
    use iso_c_binding
  
contains

    subroutine Op_getCDG(   &
        dimens,             &
        M,                  &
        N,                  &
        BL,                 &
        BU,                 &
        BCL,                &
        BCU,                &
        y1u,                &
        y2v,                &
        y3w,                &
        x1p,                &
        x2p,                &
        x3p,                &
        x1u,                &
        x2v,                &
        x3w,                &
        cdg1,               &
        cdg2,               &
        cdg3                &
        ) bind(c,name='Op_getCDG')

        implicit none

        integer(c_int), intent(in)    :: dimens

        integer(c_int), intent(in)    :: M(3)

        integer(c_int), intent(in)    :: N(3)

        integer(c_int), intent(in)    :: BL(3)
        integer(c_int), intent(in)    :: BU(3)

        integer(c_int), intent(in)    :: BCL(3)
        integer(c_int), intent(in)    :: BCU(3)

        real(c_double), intent(in)    :: y1u(0:M(1))
        real(c_double), intent(in)    :: y2v(0:M(2))
        real(c_double), intent(in)    :: y3w(0:M(3))

        real(c_double), intent(in)    :: x1p(bl(1):N(1)+bu(1))
        real(c_double), intent(in)    :: x2p(bl(2):N(2)+bu(2))
        real(c_double), intent(in)    :: x3p(bl(3):N(3)+bu(3))

        real(c_double), intent(in)    :: x1u(bl(1):N(1)+bu(1))
        real(c_double), intent(in)    :: x2v(bl(2):N(2)+bu(2))
        real(c_double), intent(in)    :: x3w(bl(3):N(3)+bu(3))

        real(c_double), intent(inout) :: cdg1(-1:1,1:N(1))
        real(c_double), intent(inout) :: cdg2(-1:1,1:N(2))
        real(c_double), intent(inout) :: cdg3(-1:1,1:N(3))

        real                          :: cDu1R(-1:0,1:N(1))
        real                          :: cDv2R(-1:0,1:N(2))
        real                          :: cDw3R(-1:0,1:N(3))

        real                          :: cGp1R( 0:1,0:N(1))
        real                          :: cGp2R( 0:1,0:N(2))
        real                          :: cGp3R( 0:1,0:N(3))

        integer                       :: i, j, k
        integer                       :: imax, jmax, kmax

        cdg1 = 0.
        cdg2 = 0.
        cdg3 = 0.

        cDu1R = 0.
        cGp1R = 0.

        !---------------------------------------------------------------------------!
        ! Achtung: Falls nicht periodisch, wird hier am Rand ein Fehler gemacht,    !
        !          der nur über die Sonderbehandlung weiter unten korrigiert wird!  !
        !---------------------------------------------------------------------------!
        !imax = NN(1,g)
        imax = N(1)

        do i = 1, imax
            cDu1R(-1,i) = -1./( x1u(i) - x1u(i-1) )
            cDu1R( 0,i) =  1./( x1u(i) - x1u(i-1) )
        end do

        do i = 1, imax-1
            cGp1R( 0,i) = -1./( x1p(i+1) - x1p(i) )
            cGp1R( 1,i) =  1./( x1p(i+1) - x1p(i) )
        end do

        !if (BC(1,1,g) > 0) then
        if (BCL(1) > 0) then
            cGp1R( :,0   ) =  0.
            cDu1R( :,1   ) =  0.
            !cDu1R( 0,1   ) =  1./(x1u (1  ) - x1u (0  ))
            cDu1R( 0,1   ) =  1./(y1u (1  ) - y1u (0  )) ! TEST!!! Der Schoenheit halber, s.u.
        else
            cGp1R( 0,0   ) = -1./(x1p(1) - x1p(0))
            cGp1R( 1,0   ) =  1./(x1p(1) - x1p(0))
        end if

        !if (BC(2,1,g) > 0) then
        if (BCU(1) > 0) then
            cGp1R( :,imax) =  0.
            cDu1R( :,imax) =  0.
            !cDu1R(-1,imax) = -1./(x1u (N(1)      ) - x1u (N(1)-1  )) ! TEST!!! Das geht in die Hose ...
            cDu1R(-1,imax) = -1./(y1u (M(1)      ) - y1u (M(1)-1  ))
        else
            cGp1R( 0,imax) = -1./( x1p(imax+1) - x1p(imax) )
            cGp1R( 1,imax) =  1./( x1p(imax+1) - x1p(imax) )
        end if
        !--------------------------------------------------------------------------------------------------------
        do i = 1, imax
            cdg1(-1,i) = cDu1R(-1,i)*cGp1R(0,i-1)
            cdg1( 0,i) = cDu1R(-1,i)*cGp1R(1,i-1) + cDu1R(0,i)*cGp1R(0,i)
            cdg1( 1,i) =                            cDu1R(0,i)*cGp1R(1,i)
        end do

        ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
        !if (BC(1,1,g) > 0) cdg1( :,1   ,g) = 2.*cdg1(:,1   ) ! TEST!!! Ist das wirklich so optimal?
        !if (BC(2,1,g) > 0) cdg1( :,imax,g) = 2.*cdg1(:,imax)
        if (BCL(1) > 0) cdg1( :,1   ) = 2.*cdg1(:,1   ) ! TEST!!! Ist das wirklich so optimal?
        if (BCU(1) > 0) cdg1( :,imax) = 2.*cdg1(:,imax)

        !if (BC(1,1,g) == -2) then
        if (BCL(1) == -2) then
            cdg1( 1,1) = cdg1( 1, 1 ) + cdg1( -1, 1 )
            cdg1(-1,1) = 0.
        end if
        !if (BC(2,1,g) == -2) then
        if (BCU(2) == -2) then
            cdg1(-1,imax) = cdg1(-1,imax) + cdg1( 1, imax )
            cdg1( 1,imax) = 0.
        end if
        !========================================================================================================
        !jmax = NN(2,g)
        jmax = N(2)

        cDv2R = 0.
        cGp2R = 0.

        do j = 1, jmax
            cDv2R(-1,j) = -1./( x2v(j) - x2v(j-1) )
            cDv2R( 0,j) =  1./( x2v(j) - x2v(j-1) )
        end do

        do j = 1, jmax-1
            cGp2R( 0,j) = -1./( x2p(j+1) - x2p(j) )
            cGp2R( 1,j) =  1./( x2p(j+1) - x2p(j) )
        end do

        !if (BC(1,2,g) > 0) then
        if (BCL(2) > 0) then
            cGp2R( :,0   ) =  0.
            cDv2R( :,1   ) =  0.
            !cDv2R( 0,1   ) =  1./(x2v (1  ) - x2v (0    ))
            cDv2R( 0,1   ) =  1./(y2v (1  ) - y2v (0    )) ! TEST!!! Der Schoenheit halber, s.u. ! TEST!!! Das Ergebnis ist nicht 100% korrekt. Warum???
        else
            cGp2R( 0,0   ) = -1./( x2p(1) - x2p(0) )
            cGp2R( 1,0   ) =  1./( x2p(1) - x2p(0) )
        end if

        !if (BC(2,2,g) > 0) then
        if (BCU(2) > 0) then
            cGp2R( :,jmax) =  0.
            cDv2R( :,jmax) =  0.
            !cDv2R(-1,jmax) = -1./( x2v (N(2)      ) - x2v (N(2)-1  ) ) ! TEST!!! Das geht in die Hose ...
            cDv2R(-1,jmax) = -1./( y2v (M(2)      ) - y2v (M(2)-1  ) )
        else
            cGp2R( 0,jmax) = -1./(x2p(jmax+1) - x2p(jmax) )
            cGp2R( 1,jmax) =  1./(x2p(jmax+1) - x2p(jmax) )
        end if
        !--------------------------------------------------------------------------------------------------------
        do j = 1, jmax
            cdg2(-1,j) = cDv2R(-1,j)*cGp2R(0,j-1)
            cdg2( 0,j) = cDv2R(-1,j)*cGp2R(1,j-1) + cDv2R(0,j)*cGp2R(0,j)
            cdg2( 1,j) =                            cDv2R(0,j)*cGp2R(1,j)
        end do

        ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
        !if (BC(1,2,g) > 0) cdg2( :,1   ) = 2.*cdg2(:,1   )
        !if (BC(2,2,g) > 0) cdg2( :,jmax) = 2.*cdg2(:,jmax)
        if (BCL(2) > 0) cdg2( :,1   ) = 2.*cdg2(:,1   )
        if (BCU(2) > 0) cdg2( :,jmax) = 2.*cdg2(:,jmax)

        !if (BC(1,2,g) == -2) then
        if (BCL(2) == -2) then
            cdg2( 1,1) = cdg2( 1, 1 ) + cdg2(-1, 1 )
            cdg2(-1,1) = 0.
        end if
        !if (BC(2,2,g) == -2) then
        if (BCU(2) == -2) then
            cdg2(-1,jmax) = cdg2(-1,jmax) + cdg2( 1,jmax)
            cdg2( 1,jmax) = 0.
        end if
        !========================================================================================================
        if (dimens == 3) then

            kmax = N(3)

            cDw3R = 0.
            cGp3R = 0.

            do k = 1, kmax
                cDw3R(-1,k) = -1./( x3w(k) - x3w(k-1) )
                cDw3R( 0,k) =  1./( x3w(k) - x3w(k-1) )
            end do

            do k = 1, kmax-1
                cGp3R( 0,k) = -1./( x3p(k+1) - x3p(k) )
                cGp3R( 1,k) =  1./( x3p(k+1) - x3p(k) )
            end do

            !if (BC(1,3,g) > 0) then
            if (BCL(3) > 0) then
                cGp3R( :,0   ) =  0.
                cDw3R( :,1   ) =  0.
                !cDw3R( 0,1   ) =  1./(x3w (1  ) - x3w (0  ))
                cDw3R( 0,1   ) =  1./(y3w (1  ) - y3w (0  )) ! TEST!!! Der Schoenheit halber, s.u.
            else
                cGp3R( 0,0   ) = -1./(x3p(1) - x3p(0))
                cGp3R( 1,0   ) =  1./(x3p(1) - x3p(0))
            end if

            !if (BC(2,3,g) > 0) then
            if (BCU(3) > 0) then
                cGp3R( :,kmax) =  0.
                cDw3R( :,kmax) =  0.
                !cDw3R(-1,kmax) = -1./(x3w (N3      ) - x3w (N3-1  )) ! TEST!!! Das geht in die Hose ...
                cDw3R(-1,kmax) = -1./(y3w (M(3)      ) - y3w (M(3)-1  ))
            else
                cGp3R( 0,kmax) = -1./(x3p(kmax+1) - x3p(kmax))
                cGp3R( 1,kmax) =  1./(x3p(kmax+1) - x3p(kmax))
            end if
            !--------------------------------------------------------------------------------------------------------
            do k = 1, kmax
                cdg3(-1,k) = cDw3R(-1,k)*cGp3R(0,k-1)
                cdg3( 0,k) = cDw3R(-1,k)*cGp3R(1,k-1) + cDw3R(0,k)*cGp3R(0,k)
                cdg3( 1,k) =                            cDw3R(0,k)*cGp3R(1,k)
            end do

            ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
            !if (BC(1,3,g) > 0) cdg3( :,1   ) = 2.*cdg3(:,1   )
            !if (BC(2,3,g) > 0) cdg3( :,kmax) = 2.*cdg3(:,kmax)
            if (BCL(3) > 0) cdg3( :,1   ) = 2.*cdg3(:,1   )
            if (BCU(3) > 0) cdg3( :,kmax) = 2.*cdg3(:,kmax)

            if (BCL(3) == -2) then
                cdg3( 1,1) = cdg3( 1,1) + cdg3(-1,1)
                cdg3(-1,1) = 0.
            end if
            if (BCU(3) == -2) then
                cdg3(-1,kmax) = cdg3(-1,kmax) + cdg3( 1,kmax)
                cdg3( 1,kmax) = 0.
            end if

        end if
           !========================================================================================================


    end subroutine Op_getCDG




    !>  \brief computes \f$ \mathrm{div = \nabla\cdot\phi } \f$
    subroutine OP_DivGrad2ndOOp(    &
        dimens,                     &
        N,                          &
        bL,bU,                      &
        BCL,BCU,                    &
!        SS,NN,                      &
        cdg1,                       &
        cdg2,                       &
        cdg3,                       &
        phi,                        &
        Lap ) bind (c,name='OP_DivGrad2ndOOp')
  
        implicit none
  
        integer(c_int), intent(in)    :: dimens

        integer(c_int), intent(in)    :: N(3)

        integer(c_int), intent(in)    :: bL(3)
        integer(c_int), intent(in)    :: bU(3)

        integer(c_int), intent(in)    :: BCL(3)
        integer(c_int), intent(in)    :: BCU(3)

        real(c_double), intent(in)    :: cdg1(-1:1,1:N(1))
        real(c_double), intent(in)    :: cdg2(-1:1,1:N(2))
        real(c_double), intent(in)    :: cdg3(-1:1,1:N(3))

        real(c_double), intent(in)    :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(inout) :: Lap (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        integer                ::  i,  S1R, N1R, S11R, N11R
        integer                ::  j,  S2R, N2R, S22R, N22R
        integer                ::  k,  S3R, N3R, S33R, N33R


            !--------------------------------------------------------------
!            SNB(1:2,1,g) = (/2+ls1,NN(1,g)+ls1/)
!            SNB(1:2,2,g) = (/2+ls2,NN(2,g)+ls2/)
!            SNB(1:2,3,g) = (/2+ls3,NN(3,g)+ls3/)

            S1R = 1
            N1R = N(1)-1

            if (BCL(1) >   0) S1R = 1
            if (BCL(1) == -2) S1R = 1
            if (BCU(1) >   0) N1R = N(1)
            if (BCU(1) == -2) N1R = N(1)

            S2R = 1
            N2R = N(2)-1

            if (BCL(2) >   0) S2R = 1
            if (BCL(2) == -2) S2R = 1
            if (BCU(2) >   0) N2R = N(2)
            if (BCU(2) == -2) N2R = N(2)

            S3R = 1
            N3R = N(3)-1

            if (BCL(3) >   0) S3R = 1
            if (BCL(3) == -2) S3R = 1
            if (BCU(3) >   0) N3R = N(3)
            if (BCU(3) == -2) N3R = N(3)
            !--------------------------------------------------------------

            !--------------------------------------------------------------
!            SNF(1:2,1,g) = (/2+ls1,NN(1,g)+ls1/)
!            SNF(1:2,2,g) = (/2+ls2,NN(2,g)+ls2/)
!            SNF(1:2,3,g) = (/2+ls3,NN(3,g)+ls3/)

            S11R = 1
            N11R = N(1)-1

            if (BCL(1) >   0) S11R = 2
            if (BCL(1) == -2) S11R = 1
            if (BCU(1) >   0) N11R = N(1)-1
            if (BCU(1) == -2) N11R = N(1)

            S22R = 1
            N22R = N(2)-1

            if (BCL(2) >   0) S22R = 2
            if (BCL(2) == -2) S2R = 1
            if (BCU(2) >   0) N22R = N(2)-1
            if (BCU(2) == -2) N22R = N(2)

            S33R = 1
            N33R = N(3)-1

            if (BCL(3) >   0) S33R = 2
            if (BCL(3) == -2) S33R = 1
            if (BCU(3) >   0) N33R = N(3)-1
            if (BCU(3) == -2) N33R = N(3)
            !--------------------------------------------------------------

!        S1R  = SNB(1,1,g)
!        S2R  = SNB(1,2,g)
!        S3R  = SNB(1,3,g)
!
!        N1R  = SNB(2,1,g)
!        N2R  = SNB(2,2,g)
!        N3R  = SNB(2,3,g)
!
!        S11R = SNF(1,1,g)
!        S22R = SNF(1,2,g)
!        S33R = SNF(1,3,g)
!
!        N11R = SNF(2,1,g)
!        N22R = SNF(2,2,g)
!        N33R = SNF(2,3,g)
!        !*****************************************************************************************
  
  
        !-----------------------------------------------------------------------------------------------------!
        ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!        !
        !            - Direktes Verbinden mit Restriktion erscheint nicht sehr attraktiv, da Lap auf dem      !
        !              jeweils feineren Gitter ohnehin benötigt wird und coarse nur per MOD( , ) oder IF-     !
        !              Abfrage belegt werden könnte. Zum anderen Wird das redundante feinere Feld auch bei    !
        !              der Interpolation verwendet (ebenfalls zur Vereinfachung der Programmierung) und ist   !
        !              somit ohnehin vorhanden. Geschwindigkeitsmässig ist vermutlich auch nicht mehr viel zu !
        !              holen.                                                                                 !
        !-----------------------------------------------------------------------------------------------------!

  
        !===========================================================================================================
        if (dimens == 3) then

            do k = S33R, N33R
                do j = S22R, N22R
                    !pgi$ unroll = n:8
                    do i = S11R, N11R
                        Lap(i,j,k) =  cdg1(-1,i)*phi(i-1,j  ,k  ) + cdg1(1,i)*phi(i+1,j  ,k  )   &
                                &  +  cdg2(-1,j)*phi(i  ,j-1,k  ) + cdg2(1,j)*phi(i  ,j+1,k  )   &
                                &  +  cdg3(-1,k)*phi(i  ,j  ,k-1) + cdg3(1,k)*phi(i  ,j  ,k+1)   &
                                 &  + (cdg1(0,i) + cdg2(0,j) + cdg3(0,k))*phi(i,j,k)
                    end do
                end do
            end do

        !===========================================================================================================
        else

            do k = S33R, N33R
                do j = S22R, N22R
                    !pgi$ unroll = n:8
                    do i = S11R, N11R
                        Lap(i,j,k) =  cdg1(-1,i)*phi(i-1,j  ,k) + cdg1(1,i)*phi(i+1,j  ,k)   &
                                &  +  cdg2(-1,j)*phi(i  ,j-1,k) + cdg2(1,j)*phi(i  ,j+1,k)   &
                                &  + (cdg1(0,i) + cdg2(0,j))*phi(i,j,k)
                    end do
                end do
            end do

        end if
        !===========================================================================================================


        !===========================================================================================================
        !=== Randbedingungen =======================================================================================
        !===========================================================================================================

!        if (BC_1L > 0) then
        if (BCL(1) > 0) then
            i = 1
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do j = S2R, N2R
                    Lap(i,j,k) = cdg1(0,i)*phi(i,j,k) + cdg1(1,i)*phi(i+1,j,k)
                end do
            end do
        end if

!        if (BC_1U > 0) then
        if (BCU(1) > 0) then
            i = N(1)
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do j = S2R, N2R
                    Lap(i,j,k) = cdg1(-1,i)*phi(i-1,j,k) + cdg1(0,i)*phi(i,j,k)
                end do
            end do
        end if

        !===========================================================================================================

!        if (BC_2L > 0) then
        if (BCL(2) > 0) then
            j = 1
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg2(0,j)*phi(i,j,k) + cdg2(1,j)*phi(i,j+1,k)
                end do
            end do
        end if

!        if (BC_2U > 0) then
        if (BCU(2) > 0) then
            j = N(2)
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg2(-1,j)*phi(i,j-1,k) + cdg2(0,j)*phi(i,j,k)
                end do
            end do
        end if

        !===========================================================================================================

!        if (BC_3L > 0) then
        if (BCL(3) > 0) then
            k = 1
            do j = S2R, N2R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg3(0,k)*phi(i,j,k) + cdg3(1,k)*phi(i,j,k+1)
                end do
            end do
        end if

!        if (BC_3U > 0) then
        if (BCU(3) > 0) then
            k = N(3)
            do j = S2R, N2R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg3(-1,k)*phi(i,j,k-1) + cdg3(0,k)*phi(i,j,k)
                end do
            end do
        end if
  
        !===========================================================================================================


!        if (corner_yes) call handle_corner_Lap(g,Lap)


    end subroutine OP_DivGrad2ndOOp

  
end module cmod_DivGrad2ndOOp
