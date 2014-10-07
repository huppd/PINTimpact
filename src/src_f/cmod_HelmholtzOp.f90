!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing many routines to compute derivatives ...
module cmod_HelmholtzOp
  
  
    use mod_vars!, only :cu11,cv22,cw33,cp11,cp22,cp33

    use iso_c_binding
  
!
  
  
contains
      !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
    !!
    !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
    !! \param[in] m dimension from one to three
    !! \param[in] mulI factor which coresponds to the factor of the identity part
    !! \param[in] mulL factor which coresponds to the factor of the laplace part
    !! \param[inout] phi
    !! \param[out] Lap
    !! \todo parameters: bl,bu,N, SS,NN,dimens, cu11,cp11,cv22,cp22, cw33, cp33 ...
    subroutine OP_helmholtz(    &
        dimens,                 &
        N,                      &
        bL,bU,                  &
        SS,NN,                  &
        c11,c22,c33,            &
        mulI, mulL,             &
        phi, Lap ) bind (c,name='OP_helmholtz')
  
        implicit none
  
        integer(c_int)  , intent(in  ) ::  dimens

        integer(c_int), intent(in)     :: N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)

        real(c_double),  intent(in)  :: c11(bL(1):bU(1),0:N(1))
        real(c_double),  intent(in)  :: c22(bL(2):bU(2),0:N(2))
        real(c_double),  intent(in)  :: c33(bL(3):bU(3),0:N(3))

        real(c_double)  , intent(in   ) ::  mulI
        real(c_double)  , intent(in   ) ::  mulL

        real(c_double),  intent(in )   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(out)   :: lap (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
        real                   ::  dd1
  

        !===========================================================================================================
        if( 3==dimens) then
            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = bL(3), bU(3)
                            dd1 = dd1 + c33(kk,k)*phi(i,j,k+kk)
                        end do
                        Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                    end do
                end do
            end do
        else
            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do
                        Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                    end do
                end do
            end do
        end if
    !===========================================================================================================

  
    end subroutine OP_helmholtz
  
!    !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
!    !!
!    !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
!    !! \param[in] m dimension from one to three
!    !! \param[in] mulI factor which coresponds to the factor of the identity part
!    !! \param[in] mulL factor which coresponds to the factor of the laplace part
!    !! \param[inout] phi
!    !! \param[out] Lap
!    !! \todo parameters: bl,bu,N, SS,NN,dimens, cu11,cp11,cv22,cp22, cw33, cp33 ...
!    subroutine OP_helmholtz( m, mulI, mulL, phi, Lap ) bind (c,name='OP_helmholtz')
!
!        implicit none
!
!        integer(c_int)  , intent(in   ) ::  m
!        real(c_double)  , intent(in   ) ::  mulI
!        real(c_double)  , intent(in   ) ::  mulL
!
!        real(c_double)  , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!        real(c_double)  , intent(  out) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!
!        integer                ::  i, ii
!        integer                ::  j, jj
!        integer                ::  k, kk
!
!        real                   ::  dd1
!
!
!
!        !===========================================================================================================
!        if (m == 1) then
!            !--------------------------------------------------------------------------------------------------------
!            if (dimens == 3) then
!                do k = S31, N31
!                    do j = S21, N21
!                        do i = S11, N11
!                            dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = S31, N31
!                    do j = S21, N21
!                        do i = S11, N11
!                            dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!           !--------------------------------------------------------------------------------------------------------
!        !===========================================================================================================
!        if (m == 2) then
!            !--------------------------------------------------------------------------------------------------------
!            if (dimens == 3) then
!                do k = S32, N32
!                    do j = S22, N22
!                        do i = S12, N12
!                            dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = S32, N32
!                    do j = S22, N22
!                        do i = S12, N12
!                            dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!           !--------------------------------------------------------------------------------------------------------
!        !===========================================================================================================
!        if (m == 3 .and. dimens == 3) then
!            !--------------------------------------------------------------------------------------------------------
!            do k = S33, N33
!                do j = S23, N23
!                    do i = S13, N13
!                        dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!                        !pgi$ unroll = n:8
!                        do ii = b1L+1, b1U
!                            dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do jj = b2L, b2U
!                            dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do kk = b3L, b3U
!                            dd1 = dd1 + cw33(kk,k)*phi(i,j,k+kk)
!                        end do
!                        Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                    end do
!                end do
!            end do
!        end if
!           !--------------------------------------------------------------------------------------------------------
!    !===========================================================================================================
!
!
!    end subroutine OP_helmholtz
  
!
!
!     !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
!    !!
!    !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
!    !! \param[in] m dimension from one to three
!    !! \param[in] mulI factor which coresponds to the factor of the identity part
!    !! \param[in] mulL factor which coresponds to the factor of the laplace part
!    !! \param[inout] phi
!    !! \param[out] Lap
!    subroutine OP_innerhelmholtz(    &
!        m,    &
!        SS,   &
!        NN,   &
!        mulI, &
!        mulL, &
!        phi,  &
!        Lap ) bind (c,name='OP_innerhelmholtz')
!
!        implicit none
!
!        integer(c_int), intent(in ) ::  m
!
!        integer(c_int), intent(in)    ::  SS(3)
!        integer(c_int), intent(in)    ::  NN(3)
!
!        real(c_double), intent(in   ) ::  mulI
!        real(c_double), intent(in   ) ::  mulL
!
!        real(c_double), intent(inout) ::  phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))
!        real(c_double), intent(inout) ::  Lap(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))
!
!        integer                ::  i, ii
!        integer                ::  j, jj
!        integer                ::  k, kk
!
!        real                   ::  dd1
!
!
!        !===========================================================================================================
!        if (m == 1) then
!            if (dimens == 3) then
!                do k = SS(3), NN(3)
!                    do j = SS(2), NN(2)
!                        do i = SS(1), NN(1)
!                            dd1 = 0.
!                            !pgi$ unroll = n:8
!                            do ii = b1L, b1U
!                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = SS(3), NN(3)
!                    do j = SS(2), NN(2)
!                        do i = SS(1), NN(1)
!                            dd1 = 0.
!                            !pgi$ unroll = n:8
!                            do ii = b1L, b1U
!                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!        !===========================================================================================================
!        if (m == 2) then
!            if (dimens == 3) then
!                do k = SS(3), NN(3)
!                    do j = SS(2), NN(2)
!                        do i = SS(1), NN(1)
!                            dd1 = 0.
!                            !pgi$ unroll = n:8
!                            do ii = b1L, b1U
!                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = SS(3), NN(3)
!                    do j = SS(2), NN(2)
!                        do i = SS(1), NN(1)
!                            dd1 = 0.
!                            !pgi$ unroll = n:8
!                            do ii = b1L, b1U
!                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!        !===========================================================================================================
!        if (m == 3 .and. dimens == 3) then
!            do k = SS(3), NN(3)
!                do j = SS(2), NN(2)
!                    do i = SS(1), NN(1)
!                        dd1 = 0.
!                        !pgi$ unroll = n:8
!                        do ii = b1L, b1U
!                            if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do jj = b2L, b2U
!                            if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do kk = b3L, b3U
!                            if( SS(3) <= k+kk .and. k+kk <= NN(3) ) dd1 = dd1 + cw33(kk,k)*phi(i,j,k+kk)
!                        end do
!                        Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                    end do
!                end do
!            end do
!        end if
!        !===========================================================================================================
!
!
!    end subroutine OP_innerhelmholtz
!
!
!
!
!    !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
!    !!
!    !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
!    !! \param[in] m dimension from one to three
!    !! \param[in] mulI factor which coresponds to the factor of the identity part
!    !! \param[in] mulL factor which coresponds to the factor of the laplace part
!    !! \param[inout] phi
!    !! \param[out] Lap
!    subroutine OP_HelmholtzGetRowEntries(    &
!        m,             &
!        SS,            &
!        NN,            &
!        mulI,          &
!        mulL,          &
!        i,             &
!        j,             &
!        k,             &
!        maxRowEntries, &
!        ic,            &
!        jc,            &
!        kc,            &
!        val,           &
!        rowEntries     &
!        ) bind (c,name='OP_HelmholtzGetRowEntries')
!
!        implicit none
!
!        integer(c_int), intent(in ) ::  m
!
!        integer(c_int), intent(in ) ::  i
!        integer(c_int), intent(in ) ::  j
!        integer(c_int), intent(in ) ::  k
!
!        integer(c_int), intent(in ) ::  SS(3)
!        integer(c_int), intent(in ) ::  NN(3)
!
!        integer(c_int), intent(in ) ::  maxRowEntries
!
!        integer(c_int), intent(out) ::  ic(maxRowEntries)
!        integer(c_int), intent(out) ::  jc(maxRowEntries)
!        integer(c_int), intent(out) ::  kc(maxRowEntries)
!
!        real(c_double), intent(out) ::  val(maxRowEntries)
!
!        integer(c_int), intent(out) ::  rowEntries
!
!        real(c_double), intent(in   ) ::  mulI
!        real(c_double), intent(in   ) ::  mulL
!
!
!        integer                ::  ii
!        integer                ::  jj
!        integer                ::  kk
!        !        integer                ::  counter
!
!
!        rowEntries = 1;
!        !===========================================================================================================
!        if (m == 1) then
!            if (dimens == 3) then
!                val(rowEntries) = mulI - mulL*( cu11(0,i)+cp22(0,j)+cp33(0,k) )
!                ic(rowEntries) = i
!                jc(rowEntries) = j
!                kc(rowEntries) = k
!                rowEntries = rowEntries + 1
!                !pgi$ unroll = n:8
!                do ii = b1L, -1
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cu11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do ii = 1, b1U
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cu11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = b2L, -1
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cp22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = 1, b2U
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cp22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do kk = b3L, -1
!                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
!                        val(rowEntries) = - mulL*( cp33(kk,k) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k+kk
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do kk = 1, b3U
!                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
!                        val(rowEntries) = - mulL*( cp33(kk,k) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k+kk
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!            else
!                val(rowEntries) = mulI - mulL*( cu11(0,i)+cp22(0,j) )
!                ic(rowEntries) = i
!                jc(rowEntries) = j
!                kc(rowEntries) = k
!                rowEntries = rowEntries + 1
!                !pgi$ unroll = n:8
!                do ii = b1L, -1
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cu11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do ii = 1, b1U
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cu11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = b2L, -1
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cp22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = 1, b2U
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cp22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!            end if
!        end if
!        !===========================================================================================================
!        if (m == 2) then
!            if (dimens == 3) then
!                val(rowEntries) = mulI - mulL*( cp11(0,i)+cv22(0,j)+cp33(0,k) )
!                ic(rowEntries) = i
!                jc(rowEntries) = j
!                kc(rowEntries) = k
!                rowEntries = rowEntries + 1
!                !pgi$ unroll = n:8
!                do ii = b1L, -1
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cp11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do ii = 1, b1U
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cp11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = b2L, -1
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cv22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = 1, b2U
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cv22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do kk = b3L, -1
!                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
!                        val(rowEntries) = - mulL*( cp33(kk,k) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k+kk
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do kk = 1, b3U
!                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
!                        val(rowEntries) = - mulL*( cp33(kk,k) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k+kk
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!            else
!                val(rowEntries) = mulI - mulL*( cp11(0,i)+cv22(0,j) )
!                ic(rowEntries) = i
!                jc(rowEntries) = j
!                kc(rowEntries) = k
!                rowEntries = rowEntries + 1
!                !pgi$ unroll = n:8
!                do ii = b1L, -1
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cp11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do ii = 1, b1U
!                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                        val(rowEntries) = - mulL*( cp11(ii,i) )
!                        ic(rowEntries) = i+ii
!                        jc(rowEntries) = j
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = b2L, -1
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cv22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!                !pgi$ unroll = n:8
!                do jj = 1, b2U
!                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                        val(rowEntries) = - mulL*( cv22(jj,j) )
!                        ic(rowEntries) = i
!                        jc(rowEntries) = j+jj
!                        kc(rowEntries) = k
!                        rowEntries = rowEntries + 1
!                    end if
!                end do
!            end if
!        end if
!        !===========================================================================================================
!        if (m == 3 .and. dimens == 3) then
!            val(rowEntries) = mulI - mulL*( cp11(0,i)+cp22(0,j)+cw33(0,k) )
!            ic(rowEntries) = i
!            jc(rowEntries) = j
!            kc(rowEntries) = k
!            rowEntries = rowEntries + 1
!            !pgi$ unroll = n:8
!            do ii = b1L, -1
!                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                    val(rowEntries) = - mulL*( cp11(ii,i) )
!                    ic(rowEntries) = i+ii
!                    jc(rowEntries) = j
!                    kc(rowEntries) = k
!                    rowEntries = rowEntries + 1
!                end if
!            end do
!            !pgi$ unroll = n:8
!            do ii = 1, b1U
!                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
!                    val(rowEntries) = - mulL*( cp11(ii,i) )
!                    ic(rowEntries) = i+ii
!                    jc(rowEntries) = j
!                    kc(rowEntries) = k
!                    rowEntries = rowEntries + 1
!                end if
!            end do
!            !pgi$ unroll = n:8
!            do jj = b2L, -1
!                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                    val(rowEntries) = - mulL*( cp22(jj,j) )
!                    ic(rowEntries) = i
!                    jc(rowEntries) = j+jj
!                    kc(rowEntries) = k
!                    rowEntries = rowEntries + 1
!                end if
!            end do
!            !pgi$ unroll = n:8
!            do jj = 1, b2U
!                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
!                    val(rowEntries) = - mulL*( cp22(jj,j) )
!                    ic(rowEntries) = i
!                    jc(rowEntries) = j+jj
!                    kc(rowEntries) = k
!                    rowEntries = rowEntries + 1
!                end if
!            end do
!            !pgi$ unroll = n:8
!            do kk = b3L, -1
!                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
!                    val(rowEntries) = - mulL*( cw33(kk,k) )
!                    ic(rowEntries) = i
!                    jc(rowEntries) = j
!                    kc(rowEntries) = k+kk
!                    rowEntries = rowEntries + 1
!                end if
!            end do
!            !pgi$ unroll = n:8
!            do kk = 1, b3U
!                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
!                    val(rowEntries) = - mulL*( cw33(kk,k) )
!                    ic(rowEntries) = i
!                    jc(rowEntries) = j
!                    kc(rowEntries) = k+kk
!                    rowEntries = rowEntries + 1
!                end if
!            end do
!        end if
!        !===========================================================================================================
!
!
!    end subroutine OP_HelmholtzGetRowEntries
!
!    !>  \brief computes \f$ \mathrm{lap_m = mulI phiI_m - mulL \Delta phiL_m} \f$
!    !!
!    !! \param[in] dimens dimension 2 or 3
!    !! \param[in] N local amount of
!    !! \param[in] m dimension from one to three
!    !! \param[in] mulI factor which coresponds to the factor of the identity part
!    !! \param[in] mulL factor which coresponds to the factor of the laplace part
!    !! \param[inout] phi
!    !! \param[out] Lap
!    subroutine OP_DtHelmholtz(   &
!        dimens,                 &
!        N,                      &
!        bL,bU,                  &
!        m,                      &
!        mulI,                   &
!        mulL,                   &
!        phiI,                   &
!        phiL,                   &
!        !    Lapc,                   &
!        Lap ) bind (c,name='OP_DtHelmholtz')
!
!        implicit none
!
!        integer(c_int), intent(in)    ::  dimens
!
!        integer(c_int), intent(in)    ::  N(3)
!
!        integer(c_int), intent(in)    ::  bL(3)
!        integer(c_int), intent(in)    ::  bU(3)
!
!        integer(c_int), intent(in   ) ::  m
!        real(c_double), intent(in   ) ::  mulI
!        real(c_double), intent(in   ) ::  mulL
!
!        real(c_double), intent(inout) ::  phiI(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
!        real(c_double), intent(inout) ::  phiL(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
!
!        real(c_double),  intent(out)  ::  Lap (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
!
!        integer                ::  i, ii
!        integer                ::  j, jj
!        integer                ::  k, kk
!
!        real                   ::  dd1
!
!        !        write(*,*) cu11(:,4)
!
!        !===========================================================================================================
!        if (m == 1) then
!            !--------------------------------------------------------------------------------------------------------
!            if (dimens == 3) then
!                do k = S31, N31
!                    do j = S21, N21
!                        do i = S11, N11
!                            dd1 = cu11(b1L,i)*phiL(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cu11(ii,i)*phiL(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cp22(jj,j)*phiL(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                dd1 = dd1 + cp33(kk,k)*phiL(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = S31, N31
!                    do j = S21, N21
!                        do i = S11, N11
!                            dd1 = cu11(b1L,i)*phiL(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cu11(ii,i)*phiL(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cp22(jj,j)*phiL(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!           !--------------------------------------------------------------------------------------------------------
!        !===========================================================================================================
!        if (m == 2) then
!            !--------------------------------------------------------------------------------------------------------
!            if (dimens == 3) then
!                do k = S32, N32
!                    do j = S22, N22
!                        do i = S12, N12
!                            dd1 = cp11(b1L,i)*phiL(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cp11(ii,i)*phiL(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cv22(jj,j)*phiL(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                dd1 = dd1 + cp33(kk,k)*phiL(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = S32, N32
!                    do j = S22, N22
!                        do i = S12, N12
!                            dd1 = cp11(b1L,i)*phiL(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cp11(ii,i)*phiL(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cv22(jj,j)*phiL(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!           !--------------------------------------------------------------------------------------------------------
!        !===========================================================================================================
!        if (m == 3 .and. dimens == 3) then
!            !--------------------------------------------------------------------------------------------------------
!            do k = S33, N33
!                do j = S23, N23
!                    do i = S13, N13
!                        dd1 = cp11(b1L,i)*phiL(i+b1L,j,k)
!                        !pgi$ unroll = n:8
!                        do ii = b1L+1, b1U
!                            dd1 = dd1 + cp11(ii,i)*phiL(i+ii,j,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do jj = b2L, b2U
!                            dd1 = dd1 + cp22(jj,j)*phiL(i,j+jj,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do kk = b3L, b3U
!                            dd1 = dd1 + cw33(kk,k)*phiL(i,j,k+kk)
!                        end do
!                        Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
!                    end do
!                end do
!            end do
!        end if
!           !--------------------------------------------------------------------------------------------------------
!    !===========================================================================================================
!
!
!    end subroutine OP_DtHelmholtz









  



  

  
end module cmod_HelmholtzOp
