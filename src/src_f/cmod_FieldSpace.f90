!>  Modul: cmod_FieldSpace
!!
!! transition module will become unnecessary in future
!! helps \c Pimpact::FieldSpace to extract varibales from impact
module cmod_FieldSpace


    use iso_c_binding
    use mpi

    use mod_dims
    use mod_vars

    implicit none

contains

    subroutine FS_get_dim( dim_ ) bind(c,name='FS_get_dim')
        implicit none
        INTEGER(c_int), intent(out)  ::  dim_

        dim_ = dimens

    END subroutine FS_get_dim


    subroutine SVS_get_nGlo( nGlo1_, nGlo2_, nGlo3_ ) bind(c,name='SVS_get_nGlo')
        implicit none
        INTEGER(c_int), intent(out)  :: nGlo1_,nGlo2_,nGlo3_

        nGlo1_ = M1
        nGlo2_ = M2
        nGlo3_ = M3

    end subroutine SVS_get_nGlo


    subroutine SVS_get_nLoc( nLoc1_, nLoc2_, nLoc3_ ) bind(c,name='SVS_get_nLoc')
        implicit none
        INTEGER(c_int), intent(out)  :: nLoc1_,nLoc2_,nLoc3_

        nLoc1_ = N1
        nLoc2_ = N2
        nLoc3_ = N3

    end subroutine SVS_get_nLoc


    subroutine SVS_get_sInd( sInd1_, sInd2_, sInd3_ ) bind(c,name='SVS_get_sInd')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S1p
        sInd2_ = S2p
        sInd3_ = S3p

    end subroutine SVS_get_sInd


    subroutine SVS_set_sInd( sInd_ ) bind(c,name='SVS_set_sInd')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd_(3)

        S1p = sInd_(1)
        S2p = sInd_(2)
        S3p = sInd_(3)

    end subroutine SVS_set_sInd



    subroutine SVS_get_eInd( eInd1_, eInd2_, eInd3_ ) bind(c,name='SVS_get_eInd')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N1p
        eInd2_ = N2p
        eInd3_ = N3p

    end subroutine SVS_get_eInd


    subroutine SVS_set_eInd( eInd_ ) bind(c,name='SVS_set_eInd')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd_(3)

        N1p = eInd_(1)
        N2p = eInd_(2)
        N3p = eInd_(3)

    end subroutine SVS_set_eInd


    subroutine VS_get_sIndU( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndU')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S11
        sInd2_ = S21
        sInd3_ = S31

    end subroutine VS_get_sIndU


    subroutine VS_get_eIndU( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndU')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N11
        eInd2_ = N21
        eInd3_ = N31

    end subroutine VS_get_eIndU


    subroutine VS_set_sIndU( sInd_ ) bind(c,name='VS_set_sIndU')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd_(3)

        S11 = sInd_(1)
        S21 = sInd_(2)
        S31 = sInd_(3)

    end subroutine VS_set_sIndU

    subroutine VS_set_eIndU( eInd_ ) bind(c,name='VS_set_eIndU')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd_(3)

        N11 = eInd_(1)
        N21 = eInd_(2)
        N31 = eInd_(3)

    end subroutine VS_set_eIndU

    subroutine VS_set_sIndUB( sInd_ ) bind(c,name='VS_set_sIndUB')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd_(3)

        S11B = sInd_(1)
        S21B = sInd_(2)
        S31B = sInd_(3)

    end subroutine VS_set_sIndUB

    subroutine VS_set_eIndUB( eInd_ ) bind(c,name='VS_set_eIndUB')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd_(3)

        N11B = eInd_(1)
        N21B = eInd_(2)
        N31B = eInd_(3)

    end subroutine VS_set_eIndUB


    subroutine get_sIndV( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndV')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S12
        sInd2_ = S22
        sInd3_ = S32

    end subroutine get_sIndV

    subroutine VS_set_sIndV( sInd_ ) bind(c,name='VS_set_sIndV')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd_(3)

        S12 = sInd_(1)
        S22 = sInd_(2)
        S32 = sInd_(3)

    end subroutine VS_set_sIndV
    subroutine VS_set_eIndV( eInd_ ) bind(c,name='VS_set_eIndV')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd_(3)

        N12 = eInd_(1)
        N22 = eInd_(2)
        N32 = eInd_(3)

    end subroutine VS_set_eIndV


    subroutine VS_set_sIndVB( sInd_ ) bind(c,name='VS_set_sIndVB')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd_(3)

        S12B = sInd_(1)
        S22B = sInd_(2)
        S32B = sInd_(3)

    end subroutine VS_set_sIndVB
    subroutine VS_set_eIndVB( eInd_ ) bind(c,name='VS_set_eIndVB')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd_(3)

        N12B = eInd_(1)
        N22B = eInd_(2)
        N32B = eInd_(3)

    end subroutine VS_set_eIndVB


    subroutine VS_get_eIndV( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndV')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N12
        eInd2_ = N22
        eInd3_ = N32

    end subroutine VS_get_eIndV


    subroutine VS_get_sIndW( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndW')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S13
        sInd2_ = S23
        sInd3_ = S33

    end subroutine VS_get_sIndW


    subroutine VS_get_eIndW( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndW')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N13
        eInd2_ = N23
        eInd3_ = N33

    end subroutine VS_get_eIndW

    subroutine VS_set_sIndW( sInd_ ) bind(c,name='VS_set_sIndW')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd_(3)

        S13 = sInd_(1)
        S23 = sInd_(2)
        S33 = sInd_(3)

    end subroutine VS_set_sIndW
    subroutine VS_set_eIndW( eInd_ ) bind(c,name='VS_set_eIndW')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd_(3)

        N13 = eInd_(1)
        N23 = eInd_(2)
        N33 = eInd_(3)

    end subroutine VS_set_eIndW

    subroutine VS_set_sIndWB( sInd_ ) bind(c,name='VS_set_sIndWB')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd_(3)

        S13B = sInd_(1)
        S23B = sInd_(2)
        S33B = sInd_(3)

    end subroutine VS_set_sIndWB

    subroutine VS_set_eIndWB( eInd_ ) bind(c,name='VS_set_eIndWB')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd_(3)

        N13B = eInd_(1)
        N23B = eInd_(2)
        N33B = eInd_(3)

    end subroutine VS_set_eIndWB


    subroutine get_sIndUB( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndUB')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S11B
        sInd2_ = S21B
        sInd3_ = S31B

    end subroutine get_sIndUB


    subroutine get_eIndUB( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndUB')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N11B
        eInd2_ = N21B
        eInd3_ = N31B

    end subroutine get_eIndUB


    subroutine get_sIndVB( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndVB')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S12B
        sInd2_ = S22B
        sInd3_ = S32B

    end subroutine get_sIndVB


    subroutine get_eIndVB( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndVB')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N12B
        eInd2_ = N22B
        eInd3_ = N32B

    end subroutine get_eIndVB


    subroutine get_sIndWB( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndWB')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S13B
        sInd2_ = S23B
        sInd3_ = S33B

    end subroutine get_sIndWB


    subroutine get_eIndWB( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndWB')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N13B
        eInd2_ = N23B
        eInd3_ = N33B

    end subroutine get_eIndWB


    subroutine SVS_get_bl( bl1_, bl2_, bl3_ ) bind(c,name='SVS_get_bl')
        implicit none
        INTEGER(c_int), intent(out)  :: bl1_,bl2_,bl3_

        bl1_ = b1L
        bl2_ = b2L
        bl3_ = b3L

    END subroutine SVS_get_bl


    subroutine SVS_get_bu( bu1_, bu2_, bu3_ ) bind(c,name='SVS_get_bu')
        implicit none
        INTEGER(c_int), intent(out)  :: bu1_,bu2_,bu3_

        bu1_ = b1U
        bu2_ = b2U
        bu3_ = b3U

    END subroutine SVS_get_bu


    subroutine SG_setCommCart( comm ) bind( c, name='SG_setCommCart' )
        implicit none
        integer(c_int), intent(in) :: comm

        COMM_CART = comm

    end subroutine SG_setCommCart

    subroutine SG_getCommCart( comm_ ) bind(c,name='SG_getCommCart')
        implicit none
        INTEGER(c_int), intent(out)  ::  comm_

        comm_ = COMM_CART

    END subroutine SG_getCommCart




    subroutine SG_getRank( rank_ ) bind( c, name='SG_getRank' )
        implicit none
        integer(c_int), intent(out) :: rank_

        rank_ = rank

    end subroutine SG_getRank

    subroutine SG_setRank( rank_ ) bind( c, name='SG_setRank' )
        implicit none
        integer(c_int), intent(in) :: rank_

        rank = rank_

    end subroutine SG_setRank


    subroutine SG_getIB( iB_ ) bind( c, name='SG_getIB' )
        implicit none
        integer(c_int), intent(out) :: iB_(1:3)
        integer :: i

        do i = 1, 3
            iB_(i) = iB(i,1)
        end do

    end subroutine SG_getIB


    subroutine SG_setIB( iB_ ) bind( c, name='SG_setIB' )
        implicit none
        integer(c_int), intent(in) :: iB_(1:3)
        integer :: i

        do i = 1, 3
            iB(i,1) = iB_(i)
        end do

    end subroutine SG_setIB


    subroutine set_Shift( shift_ ) bind( c, name='SG_setShift' )
        implicit none
        integer(c_int), intent(in) :: shift_(1:3)

        iShift = shift_(1)
        jShift = shift_(2)
        kShift = shift_(3)

    end subroutine set_Shift


    subroutine SG_setRankLU( rankl, ranku ) bind( c, name='SG_setRankLU' )
        implicit none
        integer(c_int), intent(in) :: rankl(1:3)
        integer(c_int), intent(in) :: ranku(1:3)

        rank1L = rankl(1)
        rank2L = rankl(2)
        rank3L = rankl(3)

        rank1U = ranku(1)
        rank2U = ranku(2)
        rank3U = ranku(3)

    end subroutine SG_setRankLU


    subroutine SG_getRankLU( rankl, ranku ) bind( c, name='SG_getRankLU' )
        implicit none
        integer(c_int), intent(out) :: rankl(3)
        integer(c_int), intent(out) :: ranku(3)

        rankl(1) = rank1L
        rankl(2) = rank2L
        rankl(3) = rank3L

        ranku(1) = rank1U
        ranku(2) = rank2U
        ranku(3) = rank3U

    end subroutine SG_getRankLU


    subroutine set_COMM_SLICE( slice1, slice2, slice3 ) bind( c, name='SG_setCommSlice' )
        implicit none
        integer(c_int), intent(in) :: slice1
        integer(c_int), intent(in) :: slice2
        integer(c_int), intent(in) :: slice3

        COMM_SLICE1 = slice1
        COMM_SLICE2 = slice2
        COMM_SLICE3 = slice3

    end subroutine set_COMM_SLICE


    subroutine set_COMM_BAR( bar1, bar2, bar3 ) bind( c, name='SG_setCommBar' )
        implicit none
        integer(c_int), intent(in) :: bar1
        integer(c_int), intent(in) :: bar2
        integer(c_int), intent(in) :: bar3

        COMM_BAR1 = bar1
        COMM_BAR2 = bar2
        COMM_BAR3 = bar3

    end subroutine set_COMM_BAR


    subroutine set_rank_SliceBar( rankSlice, rankBar ) bind( c, name='SG_setRankSliceBar' )
        implicit none
        integer(c_int), intent(in) :: rankSlice(1:3)
        integer(c_int), intent(in) :: rankBar(1:3)

        rank_slice1 = rankSlice(1)
        rank_slice2 = rankSlice(2)
        rank_slice3 = rankSlice(3)

        rank_bar1 = rankBar(1)
        rank_bar2 = rankBar(2)
        rank_bar3 = rankBar(3)

    end subroutine set_rank_SliceBar


    subroutine set_LS( n1_, n2_, n3_ ) bind( c, name='fsetLS' )
        implicit none
        integer(c_int)  :: n1_, n2_, n3_

        N1 = n1_
        N2 = n2_
        N3 = n3_

    end subroutine set_LS


    subroutine get_BCLoc( BC_1L_, BC_1U_, BC_2L_, BC_2U_, BC_3L_, BC_3U_ ) bind( c, name='fgetBCLoc' )

        implicit none

        integer(c_int), intent(out) ::  BC_1L_, BC_1U_
        integer(c_int), intent(out) ::  BC_2L_, BC_2U_
        integer(c_int), intent(out) ::  BC_3L_, BC_3U_

        BC_1L_ =  BC_1L
        BC_1U_ =  BC_1U

        BC_2L_ =  BC_2L
        BC_2U_ =  BC_2U

        BC_3L_ = BC_3L
        BC_3U_ = BC_3U

    end subroutine get_BCLoc


    subroutine set_BCLoc( BC_1L_, BC_1U_, BC_2L_, BC_2U_, BC_3L_, BC_3U_ ) bind( c, name='fsetBCLoc' )

        implicit none

        integer(c_int), intent(in) ::  BC_1L_, BC_1U_
        integer(c_int), intent(in) ::  BC_2L_, BC_2U_
        integer(c_int), intent(in) ::  BC_3L_, BC_3U_

        BC_1L =  BC_1L_
        BC_1U =  BC_1U_

        BC_2L =  BC_2L_
        BC_2U =  BC_2U_

        BC_3L = BC_3L_
        BC_3U = BC_3U_

    end subroutine set_BCLoc


    subroutine get_BC( BC_1L_, BC_1U_, BC_2L_, BC_2U_, BC_3L_, BC_3U_ ) bind( c, name='fgetBC' )

        implicit none

        integer(c_int), intent(out) ::  BC_1L_, BC_1U_
        integer(c_int), intent(out) ::  BC_2L_, BC_2U_
        integer(c_int), intent(out) ::  BC_3L_, BC_3U_

        BC_1L_ =  BC_1L_global
        BC_1U_ =  BC_1U_global

        BC_2L_ =  BC_2L_global
        BC_2U_ =  BC_2U_global

        BC_3L_ = BC_3L_global
        BC_3U_ = BC_3U_global

    end subroutine get_BC


    subroutine set_BC( BC_1L_, BC_1U_, BC_2L_, BC_2U_, BC_3L_, BC_3U_ ) bind( c, name='fsetBC' )

        implicit none

        integer(c_int), intent(in) ::  BC_1L_, BC_1U_
        integer(c_int), intent(in) ::  BC_2L_, BC_2U_
        integer(c_int), intent(in) ::  BC_3L_, BC_3U_

        BC_1L_global =  BC_1L_
        BC_1U_global =  BC_1U_

        BC_2L_global =  BC_2L_
        BC_2U_global =  BC_2U_

        BC_3L_global = BC_3L_
        BC_3U_global = BC_3U_

    end subroutine set_BC


    subroutine get_DS( L1_, L2_, L3_ ) bind( c, name='fgetDS' )

        implicit none

        real(c_double), intent(out) :: L1_, L2_, L3_

        L1_ = L1
        L2_ = L2
        L3_ = L3

    end subroutine get_DS

    subroutine set_DS( L1_, L2_, L3_ ) bind( c, name='fsetDS' )

        implicit none

        real(c_double), intent(in) :: L1_, L2_, L3_

        L1 = L1_
        L2 = L2_
        L3 = L3_

    end subroutine set_DS



    subroutine set_GS( n1_, n2_, n3_ ) bind( c, name='fsetGS' )

        use iso_c_binding

        use mod_vars

        implicit none

        integer(c_int)  :: n1_, n2_, n3_

        M1 = n1_
        M2 = n2_
        M3 = n3_

    end subroutine set_GS


    subroutine get_PGS( np1_, np2_, np3_ ) bind( c, name='fgetPGS' )

        implicit none

        integer(c_int), intent(out) :: np1_, np2_, np3_

        np1_ = NB1
        np2_ = NB2
        np3_ = NB3

    end subroutine get_PGS


    subroutine set_PGS( np1_, np2_, np3_ ) bind( c, name='fsetPGS' )

        implicit none

        integer(c_int), intent(in) :: np1_, np2_, np3_

        NB1 = np1_
        NB2 = np2_
        NB3 = np3_

    end subroutine set_PGS


end module cmod_FieldSpace
