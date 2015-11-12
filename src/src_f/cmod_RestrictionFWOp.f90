!> \brief module providing functions to initiliaze and apply RestrictionFWOp
module cmod_RestrictionFWOp

  use iso_c_binding

  use mpi

  implicit none

contains



  subroutine MG_restrictFW( &
      dimens,               &
      Nf,                   &
      bLf,bUf,              &
      Nc,                   &
      bLc,bUc,              &
      iimax,                &
      dd,                   &
      cR1,cR2,cR3,          &
      phif,                 &
      phic ) bind(c,name='MG_restrictFW')

    implicit none

    integer(c_int), intent(in)     :: dimens

    integer(c_int), intent(in)     :: Nf(3)

    integer(c_int), intent(in)     :: bLf(3)
    integer(c_int), intent(in)     :: bUf(3)


    integer(c_int), intent(in)     :: Nc(3)

    integer(c_int), intent(in)     :: bLc(3)
    integer(c_int), intent(in)     :: bUc(3)

    integer(c_int), intent(in)     :: iimax(3)

    integer(c_int), intent(in)     :: dd(3)

    real(c_double),  intent(in)    :: cR1 ( -1:1, 1:iimax(1) )
    real(c_double),  intent(in)    :: cR2 ( -1:1, 1:iimax(2) )
    real(c_double),  intent(in)    :: cR3 ( -1:1, 1:iimax(3) )

    real(c_double),  intent(in)    :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

    real(c_double),  intent(out)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


    integer(c_int)                ::  i, ii
    integer(c_int)                ::  j, jj
    integer(c_int)                ::  k, kk


    !----------------------------------------------------------------------------------------------------------!
    ! Anmerkungen: - für allgemeine di, dj, dk geeignet!                                                       !
    !              - überlappende Schicht in Blöcken wird (der Einfachheit halber) ebenfalls ausgetauscht, ist !
    !                aber im Prinzip redundant (genauer: phiC(S1R:N1R,S2R:N2R,S3R:N3R) = ...).               !
    !              - Motivation für diese kurze Routine ist die Möglichkeit, auch Varianten wie Full-Weighting !
    !                etc. ggf. einzubauen, ansonsten könnte sie auch eingespaart werden.                       !
    !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
    !                nicht gebraucht (erleichtert die Programmierung), so dass eine Initialisierung notwendig  !
    !                ist. Dies wiederum bedingt die INTENT(inout)-Deklaration.                                 !
    !----------------------------------------------------------------------------------------------------------!





    if (dimens == 3) then
      do kk = 1, iimax(3)
        k = dd(3)*(kk-1)+1
        do jj = 1, iimax(2) 
          j = dd(2)*(jj-1)+1
          do ii = 1, iimax(1)
            i = dd(1)*(ii-1)+1
            phic(ii,jj,kk) =      &
              &              cR1(-1,ii)*cR2(-1,jj)*cR3(-1,kk)*phif(i-1,j-1,k-1) +  &
              &              cR1( 0,ii)*cR2(-1,jj)*cR3(-1,kk)*phif(i  ,j-1,k-1) +  &
              &              cR1(+1,ii)*cR2(-1,jj)*cR3(-1,kk)*phif(i+1,j-1,k-1) +  &
              !
              &              cR1(-1,ii)*cR2( 0,jj)*cR3(-1,kk)*phif(i-1,j  ,k-1) +  &
              &              cR1( 0,ii)*cR2( 0,jj)*cR3(-1,kk)*phif(i  ,j  ,k-1) +  &
              &              cR1(+1,ii)*cR2( 0,jj)*cR3(-1,kk)*phif(i+1,j  ,k-1) +  &
              !
              &              cR1(-1,ii)*cR2(+1,jj)*cR3(-1,kk)*phif(i-1,j+1,k-1) +  &
              &              cR1( 0,ii)*cR2(+1,jj)*cR3(-1,kk)*phif(i  ,j+1,k-1) +  &
              &              cR1(+1,ii)*cR2(+1,jj)*cR3(-1,kk)*phif(i+1,j+1,k-1) +  &
              !
              !
              &              cR1(-1,ii)*cR2(-1,jj)*cR3( 0,kk)*phif(i-1,j-1,k  ) +  &
              &              cR1( 0,ii)*cR2(-1,jj)*cR3( 0,kk)*phif(i  ,j-1,k  ) +  &
              &              cR1(+1,ii)*cR2(-1,jj)*cR3( 0,kk)*phif(i+1,j-1,k  ) +  &
              !
              &              cR1(-1,ii)*cR2( 0,jj)*cR3( 0,kk)*phif(i-1,j  ,k  ) +  &
              &              cR1( 0,ii)*cR2( 0,jj)*cR3( 0,kk)*phif(i  ,j  ,k  ) +  &
              &              cR1(+1,ii)*cR2( 0,jj)*cR3( 0,kk)*phif(i+1,j  ,k  ) +  &
              !
              &              cR1(-1,ii)*cR2(+1,jj)*cR3( 0,kk)*phif(i-1,j+1,k  ) +  &
              &              cR1( 0,ii)*cR2(+1,jj)*cR3( 0,kk)*phif(i  ,j+1,k  ) +  &
              &              cR1(+1,ii)*cR2(+1,jj)*cR3( 0,kk)*phif(i+1,j+1,k  ) +  &
              !
              !
              &              cR1(-1,ii)*cR2(-1,jj)*cR3(+1,kk)*phif(i-1,j-1,k+1) +  &
              &              cR1( 0,ii)*cR2(-1,jj)*cR3(+1,kk)*phif(i  ,j-1,k+1) +  &
              &              cR1(+1,ii)*cR2(-1,jj)*cR3(+1,kk)*phif(i+1,j-1,k+1) +  &
              !
              &              cR1(-1,ii)*cR2( 0,jj)*cR3(+1,kk)*phif(i-1,j  ,k+1) +  &
              &              cR1( 0,ii)*cR2( 0,jj)*cR3(+1,kk)*phif(i  ,j  ,k+1) +  &
              &              cR1(+1,ii)*cR2( 0,jj)*cR3(+1,kk)*phif(i+1,j  ,k+1) +  &
              !
              &              cR1(-1,ii)*cR2(+1,jj)*cR3(+1,kk)*phif(i-1,j+1,k+1) +  &
              &              cR1( 0,ii)*cR2(+1,jj)*cR3(+1,kk)*phif(i  ,j+1,k+1) +  &
              &              cR1(+1,ii)*cR2(+1,jj)*cR3(+1,kk)*phif(i+1,j+1,k+1)
          end do
        end do
      end do
    else
      k  = 1
      kk = 1
      do jj = 1, iimax(2)
        j = dd(2)*(jj-1)+1
        do ii = 1, iimax(1)
          i = dd(1)*(ii-1)+1
          phic(ii,jj,kk) =    &
            !&              1*1*phif(i-1,j-1,k) +  &
          !&              2*1*phif(i  ,j-1,k) +  &
          !&              1*1*phif(i+1,j-1,k) +  &
          !!
          !&              1*2*phif(i-1,j  ,k) +  &
          !&              2*2*phif(i  ,j  ,k) +  &
          !&              1*2*phif(i+1,j  ,k) +  &
          !!
          !&              cR1(-1,ii)*cR2( 0,jj)*phif(i-1,j+1,k) +  &
          !&              cR1( 0,ii)*cR2( 0,jj)*phif(i  ,j+1,k) +  &
          !&              cR1(+1,ii)*cR2( 0,jj)*phif(i+1,j+1,k)
          &              cR1(-1,ii)*cR2(-1,jj)*phif(i-1,j-1,k) +  &
          &              cR1( 0,ii)*cR2(-1,jj)*phif(i  ,j-1,k) +  &
          &              cR1(+1,ii)*cR2(-1,jj)*phif(i+1,j-1,k) +  &
          !
          &              cR1(-1,ii)*cR2( 0,jj)*phif(i-1,j  ,k) +  &
          &              cR1( 0,ii)*cR2( 0,jj)*phif(i  ,j  ,k) +  &
          &              cR1(+1,ii)*cR2( 0,jj)*phif(i+1,j  ,k) +  &
          !
          &              cR1(-1,ii)*cR2(+1,jj)*phif(i-1,j+1,k) +  &
          &              cR1( 0,ii)*cR2(+1,jj)*phif(i  ,j+1,k) +  &
          &              cR1(+1,ii)*cR2(+1,jj)*phif(i+1,j+1,k)
        end do
      end do
    end if


  end subroutine MG_restrictFW



  !> \todo maybe use SS and NN with boundary instead should be more cleaner
  !! \todo make restrictFWion proper in 3d not just in 1d!!!
  subroutine MG_restrictFWV(    &
      dimens,                   &
      dir,                      &
      Nf,                       &
      bLf,bUf,                  &
      SSf,NNf,                  &
      Nc,                       &
      bLc,bUc,                  &
      SSc,NNc,                  &
      iimax,                    &
      dd,                       &
      cRV,                      &
      cR1,cR2,cR3,              &
      phif,                     &
      phic ) bind (c,name='MG_restrictFWV')

    implicit none

    integer(c_int), intent(in)     :: dimens

    integer(c_int), intent(in)     :: dir

    integer(c_int), intent(in)     :: Nf(1:3)

    integer(c_int), intent(in)     :: bLf(1:3)
    integer(c_int), intent(in)     :: bUf(1:3)

    integer(c_int), intent(in)     :: SSf(1:3)
    integer(c_int), intent(in)     :: NNf(1:3)

    integer(c_int), intent(in)     :: Nc(1:3)

    integer(c_int), intent(in)     :: bLc(1:3)
    integer(c_int), intent(in)     :: bUc(1:3)

    integer(c_int), intent(in)     :: SSc(1:3)
    integer(c_int), intent(in)     :: NNc(1:3)

    integer(c_int), intent(in)     :: iimax(1:3)

    integer(c_int), intent(in)     :: dd(1:3)

    real(c_double),  intent(in)    :: cRV ( 1:2, 0:iimax(dir) )

    real(c_double),  intent(in)    :: cR1 ( -1:1, 1:iimax(1) )
    real(c_double),  intent(in)    :: cR2 ( -1:1, 1:iimax(2) )
    real(c_double),  intent(in)    :: cR3 ( -1:1, 1:iimax(3) )

    real(c_double),  intent(in)    :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

    real(c_double),  intent(out)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


    integer(c_int)               ::  i, ii
    integer(c_int)               ::  j, jj
    integer(c_int)               ::  k, kk


    !----------------------------------------------------------------------------------------------------------!
    ! Anmerkungen: - Null-Setzen am Rand nicht notwendig!                                                      !
    !              - Da nur in Richtung der jeweiligen Geschwindigkeitskomponente gemittelt wird, muss nicht   !
    !                die spezialisierte Helmholtz-Variante aufgerufen werden.                                  !
    !              - Austauschrichtung ist invers zu ex1, ex2, ex3. Bei mehreren Blöcken wird auch der jeweils !
    !                redundante "überlappende" Punkt aufgrund der zu grossen Intervallgrenzen (1:iimax) zwar   !
    !                berechnet, aber aufgrund des Einweg-Austauschs falsch berechnet! Dieses Vorgehen wurde    !
    !                bislang aus übersichtsgründen vorgezogen, macht aber eine Initialisierung notwendig.      !
    !                Dabei werden Intervalle der Form 0:imax anstelle von 1:imax bearbeitet, da hier nur die   !
    !                das feinste Geschwindigkeitsgitter behandelt wird!                                        !
    !              - INTENT(inout) ist bei den feinen Gittern notwendig, da Ghost-Werte ausgetauscht werden    !
    !                müssen.                                                                                   !
    !              - Zuviele Daten werden ausgetauscht; eigentlich müsste in der Grenzfläche nur jeder 4.      !
    !                Punkt behandelt werden (4x zuviel!). Leider etwas unschön, könnte aber durch eine         !
    !                spezialisierte Austauschroutine behandelt werden, da das übergeben von Feldern mit        !
    !                Intervallen von b1L:(iimax+b1U) nur sehr schlecht funktionieren würde (d.h. mit Um-       !
    !                kopieren).                                                                                !
    !----------------------------------------------------------------------------------------------------------!

    ! TEST!!! Test schreiben, um n_gather(:,2) .GT. 1 hier zu vermeiden! Gleiches gilt natürlich für die Interpolation.



    !dd=1
    !do i=1,dimens
    !!if( 0/=NNc(i)-SSc(i) ) dd(i) = ( NNf(i)-SSf(i) )/( NNc(i)-SSc(i) )
    !if( 0/=Nc(i) ) dd(i) = ( Nf(i) )/( Nc(i)-1 )
    !end do
    !!write(*,*) dd


    !===========================================================================================================
    if( 1==dir ) then

      if( 2==dimens ) then

        kk = SSc(3)
        k = SSc(3)
        do jj = SSc(2), iimax(2)
          j = dd(2)*(jj-1)+1
          do ii = SSc(1), iimax(1)
            i = dd(1)*(ii-1)+1
            if( 1==dd(1) ) then
              phic(ii,jj,kk) = phif(i,j,k)
            else
              phic(ii,jj,kk) = &
                & cRV(1,ii)*cR2(-1,jj)*phif(i  ,j-1,k-1)  +  &
                & cRV(2,ii)*cR2(-1,jj)*phif(i+1,j-1,k-1)  +  &
                !                                             
                & cRV(1,ii)*cR2( 0,jj)*phif(i  ,j  ,k-1)  +  &
                & cRV(2,ii)*cR2( 0,jj)*phif(i+1,j  ,k-1)  +  &
                !                                             
                & cRV(1,ii)*cR2(+1,jj)*phif(i  ,j+1,k-1)  +  &
                & cRV(2,ii)*cR2(+1,jj)*phif(i+1,j+1,k-1)
            end if
          end do
        end do

      else

        do kk = SSc(3), iimax(3)
          k = dd(3)*(kk-1)+1
          do jj = SSc(2), iimax(2)
            j = dd(2)*(jj-1)+1
            do ii = SSc(1), iimax(1)
              i = dd(1)*(ii-1)+1
              if( 1==dd(1) ) then
                !phic(ii,jj,kk) = phif(i,j,k)
                phic(ii,jj,kk) = &
                  &  cR2( -1, jj)*cR3(-1,kk)*phif(i,j-1,k-1) + &
                  &  cR2(  0, jj)*cR3(-1,kk)*phif(i,j  ,k-1) + &
                  &  cR2( +1, jj)*cR3(-1,kk)*phif(i,j+1,k-1) + &
                  !                                    
                  &  cR2( -1, jj)*cR3( 0,kk)*phif(i,j-1,k  ) + &
                  &  cR2(  0, jj)*cR3( 0,kk)*phif(i,j  ,k  ) + &
                  &  cR2( +1, jj)*cR3( 0,kk)*phif(i,j+1,k  ) + &
                  !                                    
                  &  cR2( -1, jj)*cR3(+1,kk)*phif(i,j-1,k+1) + &
                  &  cR2(  0, jj)*cR3(+1,kk)*phif(i,j  ,k+1) + &
                  &  cR2( +1, jj)*cR3(+1,kk)*phif(i,j+1,k+1)    
              else
                !phic(ii,jj,kk) = cRV(1,ii)*phif(i,j,k) + cRV(2,ii)*phif(i+1,j,k)
                phic(ii,jj,kk) = &
                  & cRV(1,ii)*cR2(-1,jj)*cR3(-1,kk)*phif(i  ,j-1,k-1) + &
                  & cRV(2,ii)*cR2(-1,jj)*cR3(-1,kk)*phif(i+1,j-1,k-1) + &
                  !
                  & cRV(1,ii)*cR2( 0,jj)*cR3(-1,kk)*phif(i  ,j  ,k-1) + &
                  & cRV(2,ii)*cR2( 0,jj)*cR3(-1,kk)*phif(i+1,j  ,k-1) + &
                  !
                  & cRV(1,ii)*cR2(+1,jj)*cR3(-1,kk)*phif(i  ,j+1,k-1) + &
                  & cRV(2,ii)*cR2(+1,jj)*cR3(-1,kk)*phif(i+1,j+1,k-1) + &
                  !
                  & cRV(1,ii)*cR2(-1,jj)*cR3( 0,kk)*phif(i  ,j-1,k  ) + &
                  & cRV(2,ii)*cR2(-1,jj)*cR3( 0,kk)*phif(i+1,j-1,k  ) + &
                  !
                  & cRV(1,ii)*cR2( 0,jj)*cR3( 0,kk)*phif(i  ,j  ,k  ) + &
                  & cRV(2,ii)*cR2( 0,jj)*cR3( 0,kk)*phif(i+1,j  ,k  ) + &
                  !
                  & cRV(1,ii)*cR2(+1,jj)*cR3( 0,kk)*phif(i  ,j+1,k  ) + &
                  & cRV(2,ii)*cR2(+1,jj)*cR3( 0,kk)*phif(i+1,j+1,k  ) + &
                  !
                  & cRV(1,ii)*cR2(-1,jj)*cR3(+1,kk)*phif(i  ,j-1,k+1) + &
                  & cRV(2,ii)*cR2(-1,jj)*cR3(+1,kk)*phif(i+1,j-1,k+1) + &
                  !
                  & cRV(1,ii)*cR2( 0,jj)*cR3(+1,kk)*phif(i  ,j  ,k+1) + &
                  & cRV(2,ii)*cR2( 0,jj)*cR3(+1,kk)*phif(i+1,j  ,k+1) + &
                  !
                  & cRV(1,ii)*cR2(+1,jj)*cR3(+1,kk)*phif(i  ,j+1,k+1) + &
                  & cRV(2,ii)*cR2(+1,jj)*cR3(+1,kk)*phif(i+1,j+1,k+1)    
              end if
            end do
          end do
        end do

      end if

    end if
    !===========================================================================================================
    if( 2==dir ) then

      if( 2==dimens ) then

        kk = SSc(3)
        k = SSc(3)
        do jj = SSc(2), iimax(2)
          j = dd(2)*(jj-1)+1
          do ii = SSc(1), iimax(1)
            i = dd(1)*(ii-1)+1
            if( 1==dd(2) ) then
              phic(ii,jj,kk) = phif(i,j,k)
            else
              phic(ii,jj,kk) = cRV(1,jj)*phif(i,j,k) + cRV(2,jj)*phif(i,j+1,k)
            end if
          end do
        end do

      else

        do kk = SSc(3), iimax(3)
          k = dd(3)*(kk-1)+1
          do jj = SSc(2), iimax(2)
            j = dd(2)*(jj-1)+1
            do ii = SSc(1), iimax(1)
              i = dd(1)*(ii-1)+1
              if( 1==dd(2) ) then
                !phic(ii,jj,kk) = phif(i,j,k)
                phic(ii,jj,kk) = &
                  &  cR1( -1, ii)*cR3(-1,kk)*phif(i-1,j,k-1) + &
                  &  cR1(  0, ii)*cR3(-1,kk)*phif(i  ,j,k-1) + &
                  &  cR1( +1, ii)*cR3(-1,kk)*phif(i+1,j,k-1) + &
                  !
                  &  cR1( -1, ii)*cR3( 0,kk)*phif(i-1,j,k  ) + &
                  &  cR1(  0, ii)*cR3( 0,kk)*phif(i  ,j,k  ) + &
                  &  cR1( +1, ii)*cR3( 0,kk)*phif(i+1,j,k  ) + &
                  !
                  &  cR1( -1, ii)*cR3(+1,kk)*phif(i-1,j,k+1) + &
                  &  cR1(  0, ii)*cR3(+1,kk)*phif(i  ,j,k+1) + &
                  &  cR1( +1, ii)*cR3(+1,kk)*phif(i+1,j,k+1)    
              else
                !phic(ii,jj,kk) = cRV(1,jj)*phif(i,j,k) + cRV(2,jj)*phif(i,j+1,k)
                phic(ii,jj,kk) = &
                  & cR1(-1,ii)*cRV(1,jj)*cR3(-1,kk)*phif(i-1,j  ,k-1) + &
                  & cR1( 0,ii)*cRV(1,jj)*cR3(-1,kk)*phif(i  ,j  ,k-1) + &
                  & cR1(+1,ii)*cRV(1,jj)*cR3(-1,kk)*phif(i+1,j  ,k-1) + &
                  !
                  & cR1(-1,ii)*cRV(2,jj)*cR3(-1,kk)*phif(i-1,j+1,k-1) + &
                  & cR1( 0,ii)*cRV(2,jj)*cR3(-1,kk)*phif(i  ,j+1,k-1) + &
                  & cR1(+1,ii)*cRV(2,jj)*cR3(-1,kk)*phif(i+1,j+1,k-1) + &
                  !
                  & cR1(-1,ii)*cRV(1,jj)*cR3( 0,kk)*phif(i-1,j  ,k  ) + &
                  & cR1( 0,ii)*cRV(1,jj)*cR3( 0,kk)*phif(i  ,j  ,k  ) + &
                  & cR1(+1,ii)*cRV(1,jj)*cR3( 0,kk)*phif(i+1,j  ,k  ) + &
                  !
                  & cR1(-1,ii)*cRV(2,jj)*cR3( 0,kk)*phif(i-1,j+1,k  ) + &
                  & cR1( 0,ii)*cRV(2,jj)*cR3( 0,kk)*phif(i  ,j+1,k  ) + &
                  & cR1(+1,ii)*cRV(2,jj)*cR3( 0,kk)*phif(i+1,j+1,k  ) + &
                  !
                  & cR1(-1,ii)*cRV(1,jj)*cR3(+1,kk)*phif(i-1,j  ,k+1) + &
                  & cR1( 0,ii)*cRV(1,jj)*cR3(+1,kk)*phif(i  ,j  ,k+1) + &
                  & cR1(+1,ii)*cRV(1,jj)*cR3(+1,kk)*phif(i+1,j  ,k+1) + &
                  !
                  & cR1(-1,ii)*cRV(2,jj)*cR3(+1,kk)*phif(i-1,j+1,k+1) + &
                  & cR1( 0,ii)*cRV(2,jj)*cR3(+1,kk)*phif(i  ,j+1,k+1) + &
                  & cR1(+1,ii)*cRV(2,jj)*cR3(+1,kk)*phif(i+1,j+1,k+1)    
              end if
            end do
          end do
        end do

      end if

    end if
    !===========================================================================================================
    if( 3==dimens .and. dir==3 ) then

      do kk = SSc(3), iimax(3)
        k = dd(3)*(kk-1)+1
        do jj = SSc(2), iimax(2)
          j = dd(2)*(jj-1)+1
          do ii = SSc(1), iimax(1)
            i = dd(1)*(ii-1)+1
            if( 1==dd(3) ) then
              phic(ii,jj,kk) = &
                & cR1(-1,ii)*cR2(-1,jj)*phif(i-1,j-1,k) + &
                & cR1( 0,ii)*cR2(-1,jj)*phif(i  ,j-1,k) + &
                & cR1(+1,ii)*cR2(-1,jj)*phif(i+1,j-1,k) + &
                !
                & cR1(-1,ii)*cR2( 0,jj)*phif(i-1,j  ,k) + &
                & cR1( 0,ii)*cR2( 0,jj)*phif(i  ,j  ,k) + &
                & cR1(+1,ii)*cR2( 0,jj)*phif(i+1,j  ,k) + &
                !
                & cR1(-1,ii)*cR2(+1,jj)*phif(i-1,j+1,k) + &
                & cR1( 0,ii)*cR2(+1,jj)*phif(i  ,j+1,k) + &
                & cR1(+1,ii)*cR2(+1,jj)*phif(i+1,j+1,k)
            else
              phic(ii,jj,kk) = &
                & cR1(-1,ii)*cR2(-1,jj)*cRV(1,kk)*phif(i-1,j-1,k  ) + &
                & cR1( 0,ii)*cR2(-1,jj)*cRV(1,kk)*phif(i  ,j-1,k  ) + &
                & cR1(+1,ii)*cR2(-1,jj)*cRV(1,kk)*phif(i+1,j-1,k  ) + &
                !
                & cR1(-1,ii)*cR2(-1,jj)*cRV(2,kk)*phif(i-1,j-1,k+1) + &
                & cR1( 0,ii)*cR2(-1,jj)*cRV(2,kk)*phif(i  ,j-1,k+1) + &
                & cR1(+1,ii)*cR2(-1,jj)*cRV(2,kk)*phif(i+1,j-1,k+1) + &
                !
                & cR1(-1,ii)*cR2( 0,jj)*cRV(1,kk)*phif(i-1,j  ,k  ) + &
                & cR1( 0,ii)*cR2( 0,jj)*cRV(1,kk)*phif(i  ,j  ,k  ) + &
                & cR1(+1,ii)*cR2( 0,jj)*cRV(1,kk)*phif(i+1,j  ,k  ) + &
                !
                & cR1(-1,ii)*cR2( 0,jj)*cRV(2,kk)*phif(i-1,j  ,k+1) + &
                & cR1( 0,ii)*cR2( 0,jj)*cRV(2,kk)*phif(i  ,j  ,k+1) + &
                & cR1(+1,ii)*cR2( 0,jj)*cRV(2,kk)*phif(i+1,j  ,k+1) + &
                !
                & cR1(-1,ii)*cR2(+1,jj)*cRV(1,kk)*phif(i-1,j+1,k  ) + &
                & cR1( 0,ii)*cR2(+1,jj)*cRV(1,kk)*phif(i  ,j+1,k  ) + &
                & cR1(+1,ii)*cR2(+1,jj)*cRV(1,kk)*phif(i+1,j+1,k  ) + &
                !
                & cR1(-1,ii)*cR2(+1,jj)*cRV(2,kk)*phif(i-1,j+1,k+1) + &
                & cR1( 0,ii)*cR2(+1,jj)*cRV(2,kk)*phif(i  ,j+1,k+1) + &
                & cR1(+1,ii)*cR2(+1,jj)*cRV(2,kk)*phif(i+1,j+1,k+1)    
            end if
          end do
        end do
      end do

    end if
    !===========================================================================================================


  end subroutine MG_restrictFWV


end module cmod_RestrictionFWOp
