!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

module mod_particles
  
  
  use mod_dims
  use mod_vars
  use mod_exchange
  
  
  private
  
  public move_particles, interpol_coeffs
  
  
  INCLUDE 'mpif.h'
  
  contains
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  subroutine move_particles(exch_yes)
  
  implicit none
  
  logical, intent(in)    ::  exch_yes
  
  real                   ::  cp1L, cp1U, cu1L, cu1U
  real                   ::  cp2L, cp2U, cv2L, cv2U
  real                   ::  cp3L, cp3U, cw3L, cw3U
  
  real                   ::  bufferR(1:n_args,1:n_part_max)
  real                   ::  bufferS(1:n_args,1:n_part_max)
  
  real                   ::  mult, dd1
  real                   ::  force(1:3)
  
  integer                ::  i, ip, iu
  integer                ::  j, jp, jv
  integer                ::  k, kp, kw
  
  integer                ::  m, countR, countS
  integer                ::  status(MPI_STATUS_SIZE)
  
  integer                ::  too_many, too_many_global
  
  
  if (n_args < 14) then ! TEST!!! woanders hin ...
     if (rank == 0) write(*,*) 'ERROR! n_args < 14!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  ! particle arguments:
  !   particles(1:3,1:n_part) == particle position
  !   particles(4:6,1:n_part) == particle velocity
  !   particles(7:9,1:n_part) == particle acceleration
  !   particles(10 ,1:n_part) == particle feedback force
  !   particles(11 ,1:n_part) == particle settling velocity
  !   particles(12 ,1:n_part) == particle Stokes number
  !   particles(13 ,1:n_part) == particle number
  !   particles(14 ,1:n_part) == particle group
  !   particles(15:,1:n_part) == other particle properties
  
  
  !===========================================================================================================
  !=== ghost-cell-update =====================================================================================
  !===========================================================================================================
  ! Anmerkung: - "Normales" "exchange..." kann nicht verwendet werden, da sonst Partikel in Ecken/Kanten
  !              falsch gerechnet werden.
  call exchange_part(1,1,0,vel(b1L,b2L,b3L,1)) ! TEST!!! Es wird momentan immer in beide Richtungen kommuniziert, es waere aber jeweils nur eine Richtung notwendig (siehe auch unten)!
  call exchange_part(2,1,0,vel(b1L,b2L,b3L,1))
  call exchange_part(3,1,0,vel(b1L,b2L,b3L,1))
  
  call exchange_part(1,2,0,vel(b1L,b2L,b3L,2))
  call exchange_part(2,2,0,vel(b1L,b2L,b3L,2))
  call exchange_part(3,2,0,vel(b1L,b2L,b3L,2))
  
  call exchange_part(1,3,0,vel(b1L,b2L,b3L,3))
  call exchange_part(2,3,0,vel(b1L,b2L,b3L,3))
  call exchange_part(3,3,0,vel(b1L,b2L,b3L,3))
  !===========================================================================================================
  
  
  if (1 == 1) then ! TEST!!! Steuerung fuer passive Partikel ueberlegen (siehe auch unten) ...
                      pp(0:(N1+1),0:(N2+1),0:(N3+1)) = 0. ! TEST!!! Grenzen ok? Ganz generell, auch in "exchange_part" ...
                      rr(0:(N1+1),0:(N2+1),0:(N3+1)) = 0.
     if (dimens == 3) Ap(0:(N1+1),0:(N2+1),0:(N3+1)) = 0.
  end if
  
  
  ! TEST!!! generell: IF (dimens == 3) herausziehen!
  
  m = 1
  do while (m <= n_part)
     
     !========================================================================================================
     !=== Fluidgeschwindigkeit auf Partikelposition interpolieren ============================================
     !========================================================================================================
     call interpol_coeffs(m,cp1L,cp1U,cu1L,cu1U,cp2L,cp2U,cv2L,cv2U,cp3L,cp3U,cw3L,cw3U,ip,iu,jp,jv,kp,kw)
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
        bufferR(1,m) = cu1L*cp2L*cp3L * vel(iu-1,jp-1,kp-1,1)  &
                &    + cu1L*cp2L*cp3U * vel(iu-1,jp-1,kp  ,1)  &
                &    + cu1L*cp2U*cp3L * vel(iu-1,jp  ,kp-1,1)  &
                &    + cu1L*cp2U*cp3U * vel(iu-1,jp  ,kp  ,1)  &
                &    + cu1U*cp2L*cp3L * vel(iu  ,jp-1,kp-1,1)  &
                &    + cu1U*cp2L*cp3U * vel(iu  ,jp-1,kp  ,1)  &
                &    + cu1U*cp2U*cp3L * vel(iu  ,jp  ,kp-1,1)  &
                &    + cu1U*cp2U*cp3U * vel(iu  ,jp  ,kp  ,1)
        
        bufferR(2,m) = cp1L*cv2L*cp3L * vel(ip-1,jv-1,kp-1,2)  &
                &    + cp1L*cv2L*cp3U * vel(ip-1,jv-1,kp  ,2)  &
                &    + cp1U*cv2L*cp3L * vel(ip  ,jv-1,kp-1,2)  &
                &    + cp1U*cv2L*cp3U * vel(ip  ,jv-1,kp  ,2)  &
                &    + cp1L*cv2U*cp3L * vel(ip-1,jv  ,kp-1,2)  &
                &    + cp1L*cv2U*cp3U * vel(ip-1,jv  ,kp  ,2)  &
                &    + cp1U*cv2U*cp3L * vel(ip  ,jv  ,kp-1,2)  &
                &    + cp1U*cv2U*cp3U * vel(ip  ,jv  ,kp  ,2)
        
        bufferR(3,m) = cp1L*cp2L*cw3L * vel(ip-1,jp-1,kw-1,3)  &
                &    + cp1L*cp2U*cw3L * vel(ip-1,jp  ,kw-1,3)  &
                &    + cp1U*cp2L*cw3L * vel(ip  ,jp-1,kw-1,3)  &
                &    + cp1U*cp2U*cw3L * vel(ip  ,jp  ,kw-1,3)  &
                &    + cp1L*cp2L*cw3U * vel(ip-1,jp-1,kw  ,3)  &
                &    + cp1L*cp2U*cw3U * vel(ip-1,jp  ,kw  ,3)  &
                &    + cp1U*cp2L*cw3U * vel(ip  ,jp-1,kw  ,3)  &
                &    + cp1U*cp2U*cw3U * vel(ip  ,jp  ,kw  ,3)
     else
        kp = 2+ls3 ! == S3p == N3p
        bufferR(1,m) = cu1L*cp2L      * vel(iu-1,jp-1,kp  ,1)  &
                &    + cu1L*cp2U      * vel(iu-1,jp  ,kp  ,1)  &
                &    + cu1U*cp2L      * vel(iu  ,jp-1,kp  ,1)  &
                &    + cu1U*cp2U      * vel(iu  ,jp  ,kp  ,1)
        
        bufferR(2,m) = cp1L*cv2L      * vel(ip-1,jv-1,kp  ,2)  &
                &    + cp1U*cv2L      * vel(ip  ,jv-1,kp  ,2)  &
                &    + cp1L*cv2U      * vel(ip-1,jv  ,kp  ,2)  &
                &    + cp1U*cv2U      * vel(ip  ,jv  ,kp  ,2)
     end if
     !========================================================================================================
     
     
     i = 12
     !=======================================================================================================
     !=== particle time integration =========================================================================
     !=======================================================================================================
     if (particles(i,m) /= 0.) then
        !--- alte RHS umspeichern ---------------------------------------------------------------------------
        if (substep /= 1) then
                            bufferR(4,m) = particles(4,m)
                            bufferR(5,m) = particles(5,m)
           if (dimens == 3) bufferR(6,m) = particles(6,m)
           
                            bufferR(7,m) = particles(7,m)
                            bufferR(8,m) = particles(8,m)
           if (dimens == 3) bufferR(9,m) = particles(9,m)
        end if
        
        !--- relative particle velocity ---------------------------------------------------------------------
                         particles(7,m) = bufferR(1,m) - particles(4,m)
                         particles(8,m) = bufferR(2,m) - particles(5,m)
        if (dimens == 3) particles(9,m) = bufferR(3,m) - particles(6,m)
        
        !--- feedback-forces on fluid -----------------------------------------------------------------------
        j = 10
        k = 11
                         force(1) = -particles(j,m) * particles(7,m) / particles(k,m) ! TEST!!! Kann das nicht zur Initialisierung? (siehe auch unten ...)
                         force(2) = -particles(j,m) * particles(8,m) / particles(k,m)
        if (dimens == 3) force(3) = -particles(j,m) * particles(9,m) / particles(k,m)
        
        !--- particle inertia and settling velocity ---------------------------------------------------------
        j = 11
                         particles(7,m) = (particles(7,m) + particles(j,m)*gravity(1)) / particles(i,m)
                         particles(8,m) = (particles(8,m) + particles(j,m)*gravity(2)) / particles(i,m)
        if (dimens == 3) particles(9,m) = (particles(9,m) + particles(j,m)*gravity(3)) / particles(i,m)
        
        !--- Runge-Kutta time integration (particle position) -----------------------------------------------
        if (substep /= 1) then
           mult = dtime*bRK(substep)
                            particles(1,m) = particles(1,m) + mult*bufferR(4,m)
                            particles(2,m) = particles(2,m) + mult*bufferR(5,m)
           if (dimens == 3) particles(3,m) = particles(3,m) + mult*bufferR(6,m)
        end if
        
        mult = dtime*aRK(substep)
                         particles(1,m) = particles(1,m) + mult*particles(4,m)
                         particles(2,m) = particles(2,m) + mult*particles(5,m)
        if (dimens == 3) particles(3,m) = particles(3,m) + mult*particles(6,m)
        
        !--- Runge-Kutta time integration (particle velocity) -----------------------------------------------
        if (substep /= 1) then
           mult = dtime*bRK(substep)
                            particles(4,m) = particles(4,m) + mult*bufferR(7,m)
                            particles(5,m) = particles(5,m) + mult*bufferR(8,m)
           if (dimens == 3) particles(6,m) = particles(6,m) + mult*bufferR(9,m)
        end if
        
        mult = dtime*aRK(substep)
                         particles(4,m) = particles(4,m) + mult*particles(7,m)
                         particles(5,m) = particles(5,m) + mult*particles(8,m)
        if (dimens == 3) particles(6,m) = particles(6,m) + mult*particles(9,m)
        
     !=======================================================================================================
     else
        !--- alte RHS umspeichern ---------------------------------------------------------------------------
        if (substep /= 1) then
                            bufferR(4,m) = particles(4,m)
                            bufferR(5,m) = particles(5,m)
           if (dimens == 3) bufferR(6,m) = particles(6,m)
        end if
        
        !--- feedback-forces on fluid -----------------------------------------------------------------------
        j = 10
                         force(1) =  particles(j,m)*gravity(1) ! TEST!!! Kann das nicht zur Initialisierung? (siehe auch oben ...)
                         force(2) =  particles(j,m)*gravity(2)
        if (dimens == 3) force(3) =  particles(j,m)*gravity(3)
        
        !--- particle settling velocity ---------------------------------------------------------------------
        j = 11
                         particles(4,m) =  bufferR(1,m) + particles(j,m)*gravity(1)
                         particles(5,m) =  bufferR(2,m) + particles(j,m)*gravity(2)
        if (dimens == 3) particles(6,m) =  bufferR(3,m) + particles(j,m)*gravity(3)
        
        !--- Runge-Kutta time integration (particle position) -----------------------------------------------
        if (substep /= 1) then
           mult = dtime*bRK(substep)
                            particles(1,m) = particles(1,m) + mult*bufferR(4,m)
                            particles(2,m) = particles(2,m) + mult*bufferR(5,m)
           if (dimens == 3) particles(3,m) = particles(3,m) + mult*bufferR(6,m)
        end if
        
        mult = dtime*aRK(substep)
                         particles(1,m) = particles(1,m) + mult*particles(4,m)
                         particles(2,m) = particles(2,m) + mult*particles(5,m)
        if (dimens == 3) particles(3,m) = particles(3,m) + mult*particles(6,m)
        
     end if
     !========================================================================================================
     
     
     !========================================================================================================
     !=== apply feedback-forces on fluid =====================================================================
     !========================================================================================================
     j = 10
     if (particles(j,m) /= 0.) then
        if (dimens == 3) then
           !--------------------------------------------------------------------------------------------------
           ! Drei Alternativen:
           !-------------------
           !cp1L = cp1L/dx1p(ip-1)
           !cp1U = cp1U/dx1p(ip  )
           !cu1L = cu1L/dx1u(iu-1)
           !cu1U = cu1U/dx1u(iu  )
           !
           !cp2L = cp2L/dx2p(jp-1)
           !cp2U = cp2U/dx2p(jp  )
           !cv2L = cv2L/dx2v(jv-1)
           !cv2U = cv2U/dx2v(jv  )
           !
           !cp3L = cp3L/dx3p(kp-1)
           !cp3U = cp3U/dx3p(kp  )
           !cw3L = cw3L/dx3w(kw-1)
           !cw3U = cw3U/dx3w(kw  )
           !-------------------
           !force(1) = force(1) / (dx1u(ip-1)*dx2p(jv  )*dx3p(kw  ))
           !force(2) = force(2) / (dx2p(iu  )*dx2v(jp-1)*dx3p(kw  ))
           !force(3) = force(3) / (dx3p(iu  )*dx2p(jv  )*dx3w(kp-1))
           !-------------------
           !force(1) = force(1) / (x1p(ip)-x1p(ip-1)) / (dx2v(jv)-dx2v(jv-1)) / (dx3w(kw)-dx3w(kw-1))
           !force(2) = force(2) / (x1u(iu)-x1u(iu-1)) / (dx2p(jp)-dx2p(jp-1)) / (dx3w(kw)-dx3w(kw-1))
           !force(3) = force(3) / (x1u(iu)-x1u(iu-1)) / (dx2v(jv)-dx2v(jv-1)) / (dx3p(kp)-dx3p(kp-1))
           !--------------------------------------------------------------------------------------------------
           if (force(1) /= 0.) then
              pp(iu-1,jp-1,kp-1) = pp(iu-1,jp-1,kp-1) + cu1L*cp2L*cp3L*force(1) / (dx1u(iu-1)*dx2p(jp-1)*dx3p(kp-1))
              pp(iu-1,jp-1,kp  ) = pp(iu-1,jp-1,kp  ) + cu1L*cp2L*cp3U*force(1) / (dx1u(iu-1)*dx2p(jp-1)*dx3p(kp  ))
              pp(iu-1,jp  ,kp-1) = pp(iu-1,jp  ,kp-1) + cu1L*cp2U*cp3L*force(1) / (dx1u(iu-1)*dx2p(jp  )*dx3p(kp-1))
              pp(iu-1,jp  ,kp  ) = pp(iu-1,jp  ,kp  ) + cu1L*cp2U*cp3U*force(1) / (dx1u(iu-1)*dx2p(jp  )*dx3p(kp  ))
              pp(iu  ,jp-1,kp-1) = pp(iu  ,jp-1,kp-1) + cu1U*cp2L*cp3L*force(1) / (dx1u(iu  )*dx2p(jp-1)*dx3p(kp-1))
              pp(iu  ,jp-1,kp  ) = pp(iu  ,jp-1,kp  ) + cu1U*cp2L*cp3U*force(1) / (dx1u(iu  )*dx2p(jp-1)*dx3p(kp  ))
              pp(iu  ,jp  ,kp-1) = pp(iu  ,jp  ,kp-1) + cu1U*cp2U*cp3L*force(1) / (dx1u(iu  )*dx2p(jp  )*dx3p(kp-1))
              pp(iu  ,jp  ,kp  ) = pp(iu  ,jp  ,kp  ) + cu1U*cp2U*cp3U*force(1) / (dx1u(iu  )*dx2p(jp  )*dx3p(kp  ))
           end if
           !--------------------------------------------------------------------------------------------------
           if (force(2) /= 0.) then
              rr(ip-1,jv-1,kp-1) = rr(ip-1,jv-1,kp-1) + cp1L*cv2L*cp3L*force(2) / (dx1p(ip-1)*dx2v(jv-1)*dx3p(kp-1))
              rr(ip-1,jv-1,kp  ) = rr(ip-1,jv-1,kp  ) + cp1L*cv2L*cp3U*force(2) / (dx1p(ip-1)*dx2v(jv-1)*dx3p(kp  ))
              rr(ip  ,jv-1,kp-1) = rr(ip  ,jv-1,kp-1) + cp1U*cv2L*cp3L*force(2) / (dx1p(ip  )*dx2v(jv-1)*dx3p(kp-1))
              rr(ip  ,jv-1,kp  ) = rr(ip  ,jv-1,kp  ) + cp1U*cv2L*cp3U*force(2) / (dx1p(ip  )*dx2v(jv-1)*dx3p(kp  ))
              rr(ip-1,jv  ,kp-1) = rr(ip-1,jv  ,kp-1) + cp1L*cv2U*cp3L*force(2) / (dx1p(ip-1)*dx2v(jv  )*dx3p(kp-1))
              rr(ip-1,jv  ,kp  ) = rr(ip-1,jv  ,kp  ) + cp1L*cv2U*cp3U*force(2) / (dx1p(ip-1)*dx2v(jv  )*dx3p(kp  ))
              rr(ip  ,jv  ,kp-1) = rr(ip  ,jv  ,kp-1) + cp1U*cv2U*cp3L*force(2) / (dx1p(ip  )*dx2v(jv  )*dx3p(kp-1))
              rr(ip  ,jv  ,kp  ) = rr(ip  ,jv  ,kp  ) + cp1U*cv2U*cp3U*force(2) / (dx1p(ip  )*dx2v(jv  )*dx3p(kp  ))
           end if
           !--------------------------------------------------------------------------------------------------
           if (force(3) /= 0.) then
              Ap(ip-1,jp-1,kw-1) = Ap(ip-1,jp-1,kw-1) + cp1L*cp2L*cw3L*force(3) / (dx1p(ip-1)*dx2p(jp-1)*dx3w(kw-1))
              Ap(ip-1,jp  ,kw-1) = Ap(ip-1,jp  ,kw-1) + cp1L*cp2U*cw3L*force(3) / (dx1p(ip-1)*dx2p(jp  )*dx3w(kw-1))
              Ap(ip  ,jp-1,kw-1) = Ap(ip  ,jp-1,kw-1) + cp1U*cp2L*cw3L*force(3) / (dx1p(ip  )*dx2p(jp-1)*dx3w(kw-1))
              Ap(ip  ,jp  ,kw-1) = Ap(ip  ,jp  ,kw-1) + cp1U*cp2U*cw3L*force(3) / (dx1p(ip  )*dx2p(jp  )*dx3w(kw-1))
              Ap(ip-1,jp-1,kw  ) = Ap(ip-1,jp-1,kw  ) + cp1L*cp2L*cw3U*force(3) / (dx1p(ip-1)*dx2p(jp-1)*dx3w(kw  ))
              Ap(ip-1,jp  ,kw  ) = Ap(ip-1,jp  ,kw  ) + cp1L*cp2U*cw3U*force(3) / (dx1p(ip-1)*dx2p(jp  )*dx3w(kw  ))
              Ap(ip  ,jp-1,kw  ) = Ap(ip  ,jp-1,kw  ) + cp1U*cp2L*cw3U*force(3) / (dx1p(ip  )*dx2p(jp-1)*dx3w(kw  ))
              Ap(ip  ,jp  ,kw  ) = Ap(ip  ,jp  ,kw  ) + cp1U*cp2U*cw3U*force(3) / (dx1p(ip  )*dx2p(jp  )*dx3w(kw  ))
           end if
           !--------------------------------------------------------------------------------------------------
        else
           kp = 2+ls3 ! == S3p == N3p
           !--------------------------------------------------------------------------------------------------
           ! Drei Alternativen:
           !-------------------
           !cp1L = cp1L/dx1p(ip-1)
           !cp1U = cp1U/dx1p(ip  )
           !cu1L = cu1L/dx1u(iu-1)
           !cu1U = cu1U/dx1u(iu  )
           !
           !cp2L = cp2L/dx2p(jp-1)
           !cp2U = cp2U/dx2p(jp  )
           !cv2L = cv2L/dx2v(jv-1)
           !cv2U = cv2U/dx2v(jv  )
           !-------------------
           !force(1) = force(1) / (dx1u(ip-1)*dx2p(jv  ))
           !force(2) = force(2) / (dx2p(iu  )*dx2v(jp-1))
           !-------------------
           !force(1) = force(1) / (x1p(ip)-x1p(ip-1)) / (dx2v(jv)-dx2v(jv-1))
           !force(2) = force(2) / (x1u(iu)-x1u(iu-1)) / (dx2p(jp)-dx2p(jp-1))
           !--------------------------------------------------------------------------------------------------
           if (force(1) /= 0.) then
              pp(iu-1,jp-1,kp  ) = pp(iu-1,jp-1,kp  ) + cu1L*cp2L*force(1) / (dx1u(iu-1)*dx2p(jp-1))
              pp(iu-1,jp  ,kp  ) = pp(iu-1,jp  ,kp  ) + cu1L*cp2U*force(1) / (dx1u(iu-1)*dx2p(jp  ))
              pp(iu  ,jp-1,kp  ) = pp(iu  ,jp-1,kp  ) + cu1U*cp2L*force(1) / (dx1u(iu  )*dx2p(jp-1))
              pp(iu  ,jp  ,kp  ) = pp(iu  ,jp  ,kp  ) + cu1U*cp2U*force(1) / (dx1u(iu  )*dx2p(jp  ))
           end if
           !--------------------------------------------------------------------------------------------------
           if (force(2) /= 0.) then
              rr(ip-1,jv-1,kp  ) = rr(ip-1,jv-1,kp  ) + cp1L*cv2L*force(2) / (dx1p(ip-1)*dx2v(jv-1))
              rr(ip  ,jv-1,kp  ) = rr(ip  ,jv-1,kp  ) + cp1U*cv2L*force(2) / (dx1p(ip  )*dx2v(jv-1))
              rr(ip-1,jv  ,kp  ) = rr(ip-1,jv  ,kp  ) + cp1L*cv2U*force(2) / (dx1p(ip-1)*dx2v(jv  ))
              rr(ip  ,jv  ,kp  ) = rr(ip  ,jv  ,kp  ) + cp1U*cv2U*force(2) / (dx1p(ip  )*dx2v(jv  ))
           end if
           !--------------------------------------------------------------------------------------------------
        end if
        !-----------------------------------------------------------------------------------------------------
     end if
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Periodizitaet ======================================================================================
     !========================================================================================================
     if (BC_1L == -1) then
        if (particles(1,m) < x1p(S1p  )) particles(1,m) = particles(1,m) + (x1p(N1p)-x1p(S1p-1))
        if (particles(1,m) >= x1p(N1p+1)) particles(1,m) = particles(1,m) - (x1p(N1p)-x1p(S1p-1))
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == -1) then
        if (particles(2,m) < x2p(S2p  )) particles(2,m) = particles(2,m) + (x2p(N2p)-x2p(S2p-1))
        if (particles(2,m) >= x2p(N2p+1)) particles(2,m) = particles(2,m) - (x2p(N2p)-x2p(S2p-1))
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == -1 .and. dimens == 3) then
        if (particles(3,m) < x3p(S3p  )) particles(3,m) = particles(3,m) + (x3p(N3p)-x3p(S3p-1))
        if (particles(3,m) >= x3p(N3p+1)) particles(3,m) = particles(3,m) - (x3p(N3p)-x3p(S3p-1))
     end if
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Raender ============================================================================================
     !========================================================================================================
     ! - Partikel werden an Raendern zurueck ins Feld reflektiert, wenn entweder
     !   + Symmetrie vorliegt oder
     !   + die Advektionsgeschwindigkeit auf der Wand + die stationaere Sinkgeschwindigkeit ins Feld zeigt
     !     (gilt hier vorlaeufig fuer leichte & schwere Partikel)
     ! - Bei leichten Partikeln (Stp == 0.) sollte dieser Fall idealerweise nur in Ausnahmefaellen auftreten,
     !   z.B. bei grober Gitteraufloesung.
     ! - Bei schweren Partikeln (Stp /= 0.) kann dieser Fall generell auftreten. Auch weitere/andere
     !   Steuerungsoptionen koennten hierfuer sinnvoll sein.
     ! - Die korrekte Handhabung von Partikelgeschwindigkeit und -Beschleunigung bei Reflexion ist eigentlich
     !   generell nur schwer moeglich, weil allein die Partikelposition bei Reflexion nicht glatt in der Zeit
     !   ist und daher die Zeitableitungen entsprechend unstetig sind. Runge-Kutta Zeitintegration setzt aber
     !   Glattheit voraus. Daher wird hier kein Versuch unternommen, die Partikelgeschwindigkeit und
     !   Partikelbeschleunigung ebenfalls zu behandeln. Die Groesse des Fehlers haengt in jedem Fall von der
     !   Massentraegheit bzw. der Partikel-Stokes-Zahl ab.
     !--------------------------------------------------------------------------------------------------------
     i = 11
     !--------------------------------------------------------------------------------------------------------
     if (particles(1,m) < x1p(S1p) .and. (BC_1L > 0 .or. BC_1L == -2)) then
        if (BC_1L > 0) then
           !--- Konvektionsgeschwindigkeit auf dem Rand ---
           if (dimens == 3) then
              dd1 = cp2L*cp3L * work1(1 ,jp-1,kp-1)  & ! TEST!!! cp2L, jp, etc. beziehen sich auf ALTE Partikelposition (O1 genau?) ==> genauere Bestimmung wuenschenswert ...
               &  + cp2L*cp3U * work1(1 ,jp-1,kp  )  &
               &  + cp2U*cp3L * work1(1 ,jp  ,kp-1)  &
               &  + cp2U*cp3U * work1(1 ,jp  ,kp  )
           else
              dd1 = cp2L      * work1(1 ,jp-1,S3p )  &
               &  + cp2U      * work1(1 ,jp  ,S3p )
           end if
           dd1 = dd1 + particles(i,m)*gravity(1)
        end if
        
        !--- Spiegelung an Randflaeche ---
        if ((BC_1L > 0 .and. dd1 > 0.) .or. BC_1L == -2) particles(1,m) = 2.*x1p(S1p) - particles(1,m)
        if ((BC_1L > 0 .and. dd1 > 0.) .or. BC_1L == -2) particles(4,m) = -particles(4,m) ! TEST!!!
     end if
     if (particles(1,m) > x1p(N1p) .and. (BC_1U > 0 .or. BC_1U == -2)) then
        if (BC_1U > 0) then
           !--- Konvektionsgeschwindigkeit auf dem Rand ---
           if (dimens == 3) then
              dd1 = cp2L*cp3L * work1(N1,jp-1,kp-1)  &
               &  + cp2L*cp3U * work1(N1,jp-1,kp  )  &
               &  + cp2U*cp3L * work1(N1,jp  ,kp-1)  &
               &  + cp2U*cp3U * work1(N1,jp  ,kp  )
           else
              dd1 = cp2L      * work1(N1,jp-1,S3p )  &
               &  + cp2U      * work1(N1,jp  ,S3p )
           end if
           dd1 = dd1 + particles(i,m)*gravity(1)
        end if
        
        !--- Spiegelung an Randflaeche ---
        if ((BC_1U > 0 .and. dd1 < 0.) .or. BC_1U == -2) particles(1,m) = 2.*x1p(N1p) - particles(1,m)
        if ((BC_1U > 0 .and. dd1 < 0.) .or. BC_1U == -2) particles(4,m) = -particles(4,m) ! TEST!!!
     end if
     !--------------------------------------------------------------------------------------------------------
     if (particles(2,m) < x2p(S2p) .and. (BC_2L > 0 .or. BC_2L == -2)) then
        if (BC_2L > 0) then
           !--- Konvektionsgeschwindigkeit auf dem Rand ---
           if (dimens == 3) then
              dd1 = cp1L*cp3L * work2(ip-1,1 ,kp-1)  &
               &  + cp1L*cp3U * work2(ip-1,1 ,kp  )  &
               &  + cp1U*cp3L * work2(ip  ,1 ,kp-1)  &
               &  + cp1U*cp3U * work2(ip  ,1 ,kp  )
           else
              dd1 = cp1L      * work2(ip-1,1 ,S3p )  &
               &  + cp1U      * work2(ip  ,1 ,S3p )
           end if
           dd1 = dd1 + particles(i,m)*gravity(2)
        end if
        
        !--- Spiegelung an Randflaeche ---
        if ((BC_2L > 0 .and. dd1 > 0.) .or. BC_2L == -2) particles(2,m) = 2.*x2p(S2p) - particles(2,m)
        if ((BC_2L > 0 .and. dd1 > 0.) .or. BC_2L == -2) particles(5,m) = -particles(5,m) ! TEST!!!
     end if
     if (particles(2,m) > x2p(N2p) .and. (BC_2U > 0 .or. BC_2U == -2)) then
        if (BC_2U > 0) then
           !--- Konvektionsgeschwindigkeit auf dem Rand ---
           if (dimens == 3) then
              dd1 = cp1L*cp3L * work2(ip-1,N2,kp-1)  &
               &  + cp1L*cp3U * work2(ip-1,N2,kp  )  &
               &  + cp1U*cp3L * work2(ip  ,N2,kp-1)  &
               &  + cp1U*cp3U * work2(ip  ,N2,kp  )
           else
              dd1 = cp1L      * work2(ip-1,N2,S3p )  &
               &  + cp1U      * work2(ip  ,N2,S3p )
           end if
           dd1 = dd1 + particles(i,m)*gravity(2)
        end if
        
        !--- Spiegelung an Randflaeche ---
        if ((BC_2U > 0 .and. dd1 < 0.) .or. BC_2U == -2) particles(2,m) = 2.*x2p(N2p) - particles(2,m)
        if ((BC_2U > 0 .and. dd1 < 0.) .or. BC_2U == -2) particles(5,m) = -particles(5,m) ! TEST!!!
     end if
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
     if (particles(3,m) < x3p(S3p) .and. (BC_3L > 0 .or. BC_3L == -2)) then
        if (BC_3L > 0) then
           !--- Konvektionsgeschwindigkeit auf dem Rand ---
           dd1 = cp1L*cp2L * work3(ip-1,jp-1,1 )  &
            &  + cp1L*cp2U * work3(ip-1,jp  ,1 )  &
            &  + cp1U*cp2L * work3(ip  ,jp-1,1 )  &
            &  + cp1U*cp2U * work3(ip  ,jp  ,1 )
           
           dd1 = dd1 + particles(i,m)*gravity(3)
        end if
        
        !--- Spiegelung an Randflaeche ---
        if ((BC_3L > 0 .and. dd1 > 0.) .or. BC_3L == -2) particles(3,m) = 2.*x3p(S3p) - particles(3,m)
        if ((BC_3L > 0 .and. dd1 > 0.) .or. BC_3L == -2) particles(6,m) = -particles(6,m) ! TEST!!!
     end if
     if (particles(3,m) > x3p(N3p) .and. (BC_3U > 0 .or. BC_3U == -2)) then
        if (BC_3U > 0) then
           !--- Konvektionsgeschwindigkeit auf dem Rand ---
           dd1 = cp1L*cp2L * work3(ip-1,jp-1,N3)  &
            &  + cp1L*cp2U * work3(ip-1,jp  ,N3)  &
            &  + cp1U*cp2L * work3(ip  ,jp-1,N3)  &
            &  + cp1U*cp2U * work3(ip  ,jp  ,N3)
           
           dd1 = dd1 + particles(i,m)*gravity(3)
        end if
        
        !--- Spiegelung an Randflaeche ---
        if ((BC_3U > 0 .and. dd1 < 0.) .or. BC_3U == -2) particles(3,m) = 2.*x3p(N3p) - particles(3,m)
        if ((BC_3U > 0 .and. dd1 < 0.) .or. BC_3U == -2) particles(6,m) = -particles(6,m) ! TEST!!!
     end if
     end if
     !========================================================================================================
     
     m = m + 1
     
  end do
  
  
  !===========================================================================================================
  !=== symmetry (original) ===================================================================================
  !===========================================================================================================
  ! (sicherheitshalber, da oben exchange_part(:,:,:,vel) gerufen wird)
  if (BC_1L  == -2) vel(0 ,0:(N2+1),0:(N3+1),1) = 0.
  if (BC_1U  == -2) vel(N1,0:(N2+1),0:(N3+1),1) = 0.
  
  if (BC_2L  == -2) vel(0:(N1+1),0 ,0:(N3+1),2) = 0.
  if (BC_2U  == -2) vel(0:(N1+1),N2,0:(N3+1),2) = 0.
  if (dimens == 3) then
  if (BC_3L  == -2) vel(0:(N1+1),0:(N2+1),0 ,3) = 0.
  if (BC_3U  == -2) vel(0:(N1+1),0:(N2+1),N3,3) = 0.
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== ghost-cell-update =====================================================================================
  !===========================================================================================================
  if (1 == 1) then ! TEST!!! Steuerung fuer passive Partikel ueberlegen (siehe auch oben) ...
     !--------------------------------------------------------------------------------------------------------
     call exchange_part(1,1,1,pp) ! TEST!!! Es wird momentan immer in beide Richtungen kommuniziert, es waere aber jeweils nur eine Richtung notwendig (siehe auch oben)!
     call exchange_part(2,1,1,pp)
     call exchange_part(3,1,1,pp)
     
     nl(S11:N11,S21:N21,S31:N31,1) = nl(S11:N11,S21:N21,S31:N31,1) - pp(S11:N11,S21:N21,S31:N31)
     !--------------------------------------------------------------------------------------------------------
     call exchange_part(1,2,1,rr)
     call exchange_part(2,2,1,rr)
     call exchange_part(3,2,1,rr)
     
     nl(S12:N12,S22:N22,S32:N32,2) = nl(S12:N12,S22:N22,S32:N32,2) - rr(S12:N12,S22:N22,S32:N32)
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
     call exchange_part(1,3,1,Ap)
     call exchange_part(2,3,1,Ap)
     call exchange_part(3,3,1,Ap)
     
     nl(S13:N13,S23:N23,S33:N33,3) = nl(S13:N13,S23:N23,S33:N33,3) - Ap(S13:N13,S23:N23,S33:N33)
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  
  ! TEST!!! m-Schleifen evtl. zusammenfassen?
  !===========================================================================================================
  !=== Partikelaustausch zwischen Prozessoren ================================================================
  !===========================================================================================================
  !--- Austausch der Partikel zwischen den Prozessoren (rueckwaerts) ---
  if (BC_1L == 0) then
     m = 1
     countS = 0
     do while (m <= n_part)
        if (particles(1,m) < x1p(S1p)) then
           if (iB(1,1) == 1  ) particles(1,m) = particles(1,m) + (y1p(M1)-y1p(1))
           countS = countS + 1
           bufferS  (1:n_args,countS) = particles(1:n_args,m     )
           particles(1:n_args,m     ) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     end do
  end if
  
  if (BC_1U == 0) call MPI_IRECV(countR,1,MPI_INTEGER,rank1U,1,COMM_CART,req1U,merror)
  if (BC_1L == 0) call MPI_SEND (countS,1,MPI_INTEGER,rank1L,1,COMM_CART      ,merror)
  if (BC_1U == 0) call MPI_WAIT (req1U,status,merror)
  if (BC_1U == 0) call pseudocall_int(countR)
  
  too_many = 0
  if (n_part+countR > n_part_max) too_many = 1
  call MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (too_many_global > 0) then
     if (rank == 0) write(*,*) 'ERROR! Too many particles on at least one processor!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  if (BC_1U == 0) call MPI_IRECV(bufferR,n_args*countR,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
  if (BC_1L == 0) call MPI_SEND (bufferS,n_args*countS,MPI_REAL8,rank1L,2,COMM_CART      ,merror)
  if (BC_1U == 0) call MPI_WAIT (req1U,status,merror)
  if (BC_1U == 0) call pseudocall(bufferR)
  
  if (BC_1U == 0) then
     do m = 1, countR
        particles(1:n_args,n_part+m) = bufferR(1:n_args,m)
     end do
     n_part = n_part+countR
  end if
  
  !--- Austausch der Partikel zwischen den Prozessoren (vorwaerts) ---
  if (BC_1U == 0) then
     m = 1
     countS = 0
     do while (m <= n_part)
        if (particles(1,m) >= x1p(N1p+1)) then
           if (iB(1,1) == NB1) particles(1,m) = particles(1,m) - (y1p(M1)-y1p(1))
           countS = countS + 1
           bufferS  (1:n_args,countS) = particles(1:n_args,m     )
           particles(1:n_args,m     ) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     end do
  end if
  
  if (BC_1L == 0) call MPI_IRECV(countR,1,MPI_INTEGER,rank1L,1,COMM_CART,req1L,merror)
  if (BC_1U == 0) call MPI_SEND (countS,1,MPI_INTEGER,rank1U,1,COMM_CART      ,merror)
  if (BC_1L == 0) call MPI_WAIT (req1L,status,merror)
  if (BC_1L == 0) call pseudocall_int(countR)
  
  too_many = 0
  if (n_part+countR > n_part_max) too_many = 1
  call MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (too_many_global > 0) then
     if (rank == 0) write(*,*) 'ERROR! Too many particles on at least one processor!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  if (BC_1L == 0) call MPI_IRECV(bufferR,n_args*countR,MPI_REAL8,rank1L,2,COMM_CART,req1L,merror)
  if (BC_1U == 0) call MPI_SEND (bufferS,n_args*countS,MPI_REAL8,rank1U,2,COMM_CART      ,merror)
  if (BC_1L == 0) call MPI_WAIT (req1L,status,merror)
  if (BC_1L == 0) call pseudocall(bufferR)
  
  if (BC_1L == 0) then
     do m = 1, countR
        particles(1:n_args,n_part+m) = bufferR(1:n_args,m)
     end do
     n_part = n_part+countR
  end if
  !===========================================================================================================
  !--- Austausch der Partikel zwischen den Prozessoren (rueckwaerts) ---
  if (BC_2L == 0) then
     m = 1
     countS = 0
     do while (m <= n_part)
        if (particles(2,m) < x2p(S2p)) then
           if (iB(2,1) == 1  ) particles(2,m) = particles(2,m) + (y2p(M2)-y2p(1))
           countS = countS + 1
           bufferS  (1:n_args,countS) = particles(1:n_args,m     )
           particles(1:n_args,m     ) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     end do
  end if
  
  if (BC_2U == 0) call MPI_IRECV(countR,1,MPI_INTEGER,rank2U,1,COMM_CART,req2U,merror)
  if (BC_2L == 0) call MPI_SEND (countS,1,MPI_INTEGER,rank2L,1,COMM_CART      ,merror)
  if (BC_2U == 0) call MPI_WAIT (req2U,status,merror)
  if (BC_2U == 0) call pseudocall_int(countR)
  
  too_many = 0
  if (n_part+countR > n_part_max) too_many = 1
  call MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (too_many_global > 0) then
     if (rank == 0) write(*,*) 'ERROR! Too many particles on at least one processor!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  if (BC_2U == 0) call MPI_IRECV(bufferR,n_args*countR,MPI_REAL8,rank2U,2,COMM_CART,req2U,merror)
  if (BC_2L == 0) call MPI_SEND (bufferS,n_args*countS,MPI_REAL8,rank2L,2,COMM_CART      ,merror)
  if (BC_2U == 0) call MPI_WAIT (req2U,status,merror)
  if (BC_2U == 0) call pseudocall(bufferR)
  
  if (BC_2U == 0) then
     do m = 1, countR
        particles(1:n_args,n_part+m) = bufferR(1:n_args,m)
     end do
     n_part = n_part+countR
  end if
  
  !--- Austausch der Partikel zwischen den Prozessoren (vorwaerts) ---
  if (BC_2U == 0) then
     m = 1
     countS = 0
     do while (m <= n_part)
        if (particles(2,m) >= x2p(N2p+1)) then
           if (iB(2,1) == NB2) particles(2,m) = particles(2,m) - (y2p(M2)-y2p(1))
           countS = countS + 1
           bufferS  (1:n_args,countS) = particles(1:n_args,m     )
           particles(1:n_args,m     ) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     end do
  end if
  
  if (BC_2L == 0) call MPI_IRECV(countR,1,MPI_INTEGER,rank2L,1,COMM_CART,req2L,merror)
  if (BC_2U == 0) call MPI_SEND (countS,1,MPI_INTEGER,rank2U,1,COMM_CART      ,merror)
  if (BC_2L == 0) call MPI_WAIT (req2L,status,merror)
  if (BC_2L == 0) call pseudocall_int(countR)
  
  too_many = 0
  if (n_part+countR > n_part_max) too_many = 1
  call MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (too_many_global > 0) then
     if (rank == 0) write(*,*) 'ERROR! Too many particles on at least one processor!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  if (BC_2L == 0) call MPI_IRECV(bufferR,n_args*countR,MPI_REAL8,rank2L,2,COMM_CART,req2L,merror)
  if (BC_2U == 0) call MPI_SEND (bufferS,n_args*countS,MPI_REAL8,rank2U,2,COMM_CART      ,merror)
  if (BC_2L == 0) call MPI_WAIT (req2L,status,merror)
  if (BC_2L == 0) call pseudocall(bufferR)
  
  if (BC_2L == 0) then
     do m = 1, countR
        particles(1:n_args,n_part+m) = bufferR(1:n_args,m)
     end do
     n_part = n_part+countR
  end if
  !===========================================================================================================
  if (dimens == 3) then
  !--- Austausch der Partikel zwischen den Prozessoren (rueckwaerts) ---
  if (BC_3L == 0) then
     m = 1
     countS = 0
     do while (m <= n_part)
        if (particles(3,m) < x3p(S3p)) then
           if (iB(3,1) == 1  ) particles(3,m) = particles(3,m) + (y3p(M3)-y3p(1))
           countS = countS + 1
           bufferS  (1:n_args,countS) = particles(1:n_args,m     )
           particles(1:n_args,m     ) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     end do
  end if
  
  if (BC_3U == 0) call MPI_IRECV(countR,1,MPI_INTEGER,rank3U,1,COMM_CART,req3U,merror)
  if (BC_3L == 0) call MPI_SEND (countS,1,MPI_INTEGER,rank3L,1,COMM_CART      ,merror)
  if (BC_3U == 0) call MPI_WAIT (req3U,status,merror)
  if (BC_3U == 0) call pseudocall_int(countR)
  
  too_many = 0
  if (n_part+countR > n_part_max) too_many = 1
  call MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (too_many_global > 0) then
     if (rank == 0) write(*,*) 'ERROR! Too many particles on at least one processor!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  if (BC_3U == 0) call MPI_IRECV(bufferR,n_args*countR,MPI_REAL8,rank3U,2,COMM_CART,req3U,merror)
  if (BC_3L == 0) call MPI_SEND (bufferS,n_args*countS,MPI_REAL8,rank3L,2,COMM_CART      ,merror)
  if (BC_3U == 0) call MPI_WAIT (req3U,status,merror)
  if (BC_3U == 0) call pseudocall(bufferR)
  
  if (BC_3U == 0) then
     do m = 1, countR
        particles(1:n_args,n_part+m) = bufferR(1:n_args,m)
     end do
     n_part = n_part+countR
  end if
  
  !--- Austausch der Partikel zwischen den Prozessoren (vorwaerts) ---
  if (BC_3U == 0) then
     m = 1
     countS = 0
     do while (m <= n_part)
        if (particles(3,m) >= x3p(N3p+1)) then
           if (iB(3,1) == NB3) particles(3,m) = particles(3,m) - (y3p(M3)-y3p(1))
           countS = countS + 1
           bufferS  (1:n_args,countS) = particles(1:n_args,m     )
           particles(1:n_args,m     ) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     end do
  end if
  
  if (BC_3L == 0) call MPI_IRECV(countR,1,MPI_INTEGER,rank3L,1,COMM_CART,req3L,merror)
  if (BC_3U == 0) call MPI_SEND (countS,1,MPI_INTEGER,rank3U,1,COMM_CART      ,merror)
  if (BC_3L == 0) call MPI_WAIT (req3L,status,merror)
  if (BC_3L == 0) call pseudocall_int(countR)
  
  too_many = 0
  if (n_part+countR > n_part_max) too_many = 1
  call MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (too_many_global > 0) then
     if (rank == 0) write(*,*) 'ERROR! Too many particles on at least one processor!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  if (BC_3L == 0) call MPI_IRECV(bufferR,n_args*countR,MPI_REAL8,rank3L,3,COMM_CART,req3L,merror)
  if (BC_3U == 0) call MPI_SEND (bufferS,n_args*countS,MPI_REAL8,rank3U,3,COMM_CART      ,merror)
  if (BC_3L == 0) call MPI_WAIT (req3L,status,merror)
  if (BC_3L == 0) call pseudocall(bufferR)
  
  if (BC_3L == 0) then
     do m = 1, countR
        particles(1:n_args,n_part+m) = bufferR(1:n_args,m)
     end do
     n_part = n_part+countR
  end if
  end if
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Transport aus der Domain / Loeschen ===================================================================
  !===========================================================================================================
  ! - Wird nach dem Prozessoraustausch durchgefuehrt, damit die Sub-Domain Zuordnung fuer die Deposit-
  !   Berechnung aktuell ist.
  ! - Feldgroessen von dep1L_part etc. sollten gross genug sein (siehe Box).
  !-----------------------------------------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !Beachte (Partikelzuordnung zwischen den Prozessoren):
  !-----------------------------------------------------
  !
  !FIND1: DO i = S1p, N1p+1
  !   IF (particles(1,m) .LT. x1p(i)) THEN
  !      ip = i
  !      EXIT
  !
  !  ==>  ...   x1p(S1p) <= particles(1,m) < x1p(N1p+1) == x1p(S1p) <= x1p(S1p+1)  ...
  !       ...       S1p  <= ip             <     N1p+1  ==     S1p  <=     S1p+1   ...
  !
  ! mit  S1p = 2  + ls1 = 1
  !      N1p = N1 + ls1 = N1-1
  !
  !  ==>  ...   x1p(2) <= particles(1,m) < x1p(N1) == x1p(2) <= x1p(3)  ...
  !       ...       2  <= ip             <     N1  ==     2  <=     3   ...
  !
  !                       --------------------
  !                  ==>  ! 2 <= ip   < N1   !
  !                       ! 1 <= ip-1 < N1-1 !  etc. ...
  !                       --------------------
  !                     (gilt nur fuer ls1 = -1)
  !-----------------------------------------------------------------------------
  
  if (ls1 /= -1 .or. ls2 /= -1 .or. ls3 /= -1) then
     if (rank == 0) write(*,*) 'ERROR! Computation of particle deposit is only possible for ls1 = ls2 = ls3 = -1'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  i   = 10
  j   = 11
  k   = 12
  m   = 1
  dd1 = 0.
  do while (m <= n_part)
     
     ! TEST!!! Deposit-Integration wird momentan nur fuer Richtung 1 ausgefuehrt (analog zu conc)
     ! Anmerkung: - Peroidizitaet wird korrekt behandelt (siehe oben).
     !            - streng genommen ist das keine Zeitintegration:
     if (dimens == 3) then
        if ((particles(1,m) < x1p(S1p) .and. BC_1L > 0) .or.  &
          & (particles(1,m) > x1p(N1p) .and. BC_1U > 0) .or.  &
          & (particles(2,m) < x2p(S2p) .and. BC_2L > 0) .or.  &
          & (particles(2,m) > x2p(N2p) .and. BC_2U > 0) .or.  &
          & (particles(3,m) < x3p(S3p) .and. BC_3L > 0) .or.  &
          & (particles(3,m) > x3p(N3p) .and. BC_3U > 0)) then
           
           !--- aktuelle Gewichte ---
           call interpol_coeffs(m,cp1L,cp1U,cu1L,cu1U,cp2L,cp2U,cv2L,cv2U,cp3L,cp3U,cw3L,cw3U,ip,iu,jp,jv,kp,kw)
           
           if (particles(1,m) < y1p(1) .and.                                   &
            &  particles(2,m) >= y2p(1) .and. particles(2,m) <= y2p(M2) .and. &
            &  particles(3,m) >= y3p(1) .and. particles(3,m) <= y3p(M3)) then
              dep1L_part(jp-1,kp-1) = dep1L_part(jp-1,kp-1) + cp2L*cp3L
              dep1L_part(jp-1,kp  ) = dep1L_part(jp-1,kp  ) + cp2L*cp3U
              dep1L_part(jp  ,kp-1) = dep1L_part(jp  ,kp-1) + cp2U*cp3L
              dep1L_part(jp  ,kp  ) = dep1L_part(jp  ,kp  ) + cp2U*cp3U
           end if
           
           !--- Verlust kinetischer Energie ---
           dd1 = dd1 + particles(i,m)/particles(j,m)*particles(k,m)*(particles(4,m)**2 + particles(5,m)**2 + particles(6,m)**2)
           
           !--- L�schen des hinaustransprortierten Partikels bzw. Umsortieren ---
           particles(1:n_args,m) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     else
        if ((particles(1,m) < x1p(S1p) .and. BC_1L > 0) .or.  &
          & (particles(1,m) > x1p(N1p) .and. BC_1U > 0) .or.  &
          & (particles(2,m) < x2p(S2p) .and. BC_2L > 0) .or.  &
          & (particles(2,m) > x2p(N2p) .and. BC_2U > 0)) then
           
           !--- aktuelle Gewichte ---
           call interpol_coeffs(m,cp1L,cp1U,cu1L,cu1U,cp2L,cp2U,cv2L,cv2U,cp3L,cp3U,cw3L,cw3U,ip,iu,jp,jv,kp,kw)
           
           if (particles(1,m) < y1p(1) .and.                                   &
            &  particles(2,m) >= y2p(1) .and. particles(2,m) <= y2p(M2)) then
              dep1L_part(jp-1,S3p ) = dep1L_part(jp-1,S3p ) + cp2L
              dep1L_part(jp  ,S3p ) = dep1L_part(jp  ,S3p ) + cp2U
           end if
           
           !--- Verlust kinetischer Energie ---
           dd1 = dd1 + particles(i,m)/particles(j,m)*particles(k,m)*(particles(4,m)**2 + particles(5,m)**2)
           
           !--- L�schen des hinaustransprortierten Partikels bzw. Umsortieren ---
           particles(1:n_args,m) = particles(1:n_args,n_part)
           n_part = n_part - 1
        else
           m = m + 1
        end if
     end if
  end do
  
  
  ! TEST!!! das ist noch keine so tolle L�sung ...
  ! TEST!!! nur eine Spezies (i.e. n_spec = 1)
  !--- Verlust kinetischer Energie ---
  call MPI_REDUCE(dd1,mult,1,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
  
  if (rank == 0) then
     Ekin_part = Ekin_part - mult / 2.
     !IF (substep /= 1) Ekin_part = Ekin_part - bRK(substep)*mult / 2. ! TEST!!! "mult" m�sste global gespeichert werden ...
     !                  Ekin_part = Ekin_part - aRK(substep)*mult / 2.
  end if
  !===========================================================================================================
  
  
  end subroutine move_particles
  
  
  
  
  
  
  
  
  
  
  ! Anmerkung: Fuer gute Performance koennte Inlining dieser Routine entscheidend sein.
  subroutine interpol_coeffs(m,cp1L,cp1U,cu1L,cu1U,cp2L,cp2U,cv2L,cv2U,cp3L,cp3U,cw3L,cw3U,ip,iu,jp,jv,kp,kw)
  
  implicit none
  
  integer, intent(in   ) ::  m
  
  real   , intent(out  ) ::  cp1L, cp1U, cu1L, cu1U
  real   , intent(out  ) ::  cp2L, cp2U, cv2L, cv2U
  real   , intent(out  ) ::  cp3L, cp3U, cw3L, cw3U
  
  integer, intent(out  ) ::  ip, iu
  integer, intent(out  ) ::  jp, jv
  integer, intent(out  ) ::  kp, kw
  integer                ::  i, j, k
 
 
  ! Anmerkung: Evtl. lohnt es sich Geschwindigkeitsm�ssig nur um die vorherigen Koordinaten herum zu suchen ...
  !-----------------------------------------------------------------------------------------------------------
  FIND1: do i = S1p, N1p+1
     if (particles(1,m) < x1p(i)) then
        ip = i
        
        if (particles(1,m) < x1u(ip-1)) then
           iu = ip-1
        else
           iu = ip
        end if
        
        cu1L = (x1u(iu)-particles(1,m)) / (x1u(iu)-x1u(iu-1))
        cp1L = (x1p(ip)-particles(1,m)) / (x1p(ip)-x1p(ip-1))
        
        cu1U = 1.-cu1L
        cp1U = 1.-cp1L
        
        exit FIND1
     end if
  end do FIND1
  !-----------------------------------------------------------------------------------------------------------
  FIND2: do j = S2p, N2p+1
     if (particles(2,m) < x2p(j)) then
        jp = j
        
        if (particles(2,m) < x2v(jp-1)) then
           jv = jp-1
        else
           jv = jp
        end if
        
        cv2L = (x2v(jv)-particles(2,m)) / (x2v(jv)-x2v(jv-1))
        cp2L = (x2p(jp)-particles(2,m)) / (x2p(jp)-x2p(jp-1))
        
        cv2U = 1.-cv2L
        cp2U = 1.-cp2L
        
        exit FIND2
     end if
  end do FIND2
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
  FIND3: do k = S3p, N3p+1
     if (particles(3,m) < x3p(k)) then
        kp = k
        
        if (particles(3,m) < x3w(kp-1)) then
           kw = kp-1
        else
           kw = kp
        end if
        
        cw3L = (x3w(kw)-particles(3,m)) / (x3w(kw)-x3w(kw-1))
        cp3L = (x3p(kp)-particles(3,m)) / (x3p(kp)-x3p(kp-1))
        
        cw3U = 1.-cw3L
        cp3U = 1.-cp3L
        
        exit FIND3
     end if
  end do FIND3
  !ELSE ! Sollte eigentlich nicht notwendig sein, solange sauber nach "dimens" unterschieden wird!
  !  kp = 1
  !  kw = 1
  !  
  !  cw3L = 0.
  !  cp3L = 0.
  !      
  !  cw3U = 0.
  !  cp3U = 0.
  end if
  !-----------------------------------------------------------------------------------------------------------
  
  
  end subroutine interpol_coeffs
  
  
  
end module mod_particles
