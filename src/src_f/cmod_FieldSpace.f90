!>  Modul: cmod_FieldSpace
!!
!! helps \c Pimpact::FieldSpace to extract varibales from impact
module cmod_FieldSpace


    use iso_c_binding
    use mpi


    use mod_dims
    use mod_vars


    implicit none

contains

    subroutine get_comm( comm_ ) bind(c,name='SVS_get_comm')
        implicit none
        INTEGER(c_int), intent(out)  ::  comm_

        comm_ = COMM_CART

    END subroutine get_comm


    subroutine get_dim( dim_ ) bind(c,name='FS_get_dim')
        implicit none
        INTEGER(c_int), intent(out)  ::  dim_

        dim_ = dimens

    END subroutine get_dim


    subroutine get_nGlo( nGlo1_, nGlo2_, nGlo3_ ) bind(c,name='SVS_get_nGlo')
        implicit none
        INTEGER(c_int), intent(out)  :: nGlo1_,nGlo2_,nGlo3_

        nGlo1_ = M1
        nGlo2_ = M2
        nGlo3_ = M3

    end subroutine get_nGlo


    subroutine get_nLoc( nLoc1_, nLoc2_, nLoc3_ ) bind(c,name='SVS_get_nLoc')
        implicit none
        INTEGER(c_int), intent(out)  :: nLoc1_,nLoc2_,nLoc3_

        nLoc1_ = N1
        nLoc2_ = N2
        nLoc3_ = N3

    end subroutine get_nLoc


    subroutine get_sInd( sInd1_, sInd2_, sInd3_ ) bind(c,name='SVS_get_sInd')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S1p
        sInd2_ = S2p
        sInd3_ = S3p

    end subroutine get_sInd


    subroutine get_eInd( eInd1_, eInd2_, eInd3_ ) bind(c,name='SVS_get_eInd')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N1p
        eInd2_ = N2p
        eInd3_ = N3p

    end subroutine get_eInd


    subroutine get_sIndU( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndU')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S11
        sInd2_ = S21
        sInd3_ = S31

    end subroutine get_sIndU


    subroutine get_eIndU( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndU')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N11
        eInd2_ = N21
        eInd3_ = N31

    end subroutine get_eIndU


    subroutine get_sIndV( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndV')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S12
        sInd2_ = S22
        sInd3_ = S32

    end subroutine get_sIndV


    subroutine get_eIndV( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndV')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N12
        eInd2_ = N22
        eInd3_ = N32

    end subroutine get_eIndV


    subroutine get_sIndW( sInd1_, sInd2_, sInd3_ ) bind(c,name='VS_get_sIndW')
        implicit none
        INTEGER(c_int), intent(out)  :: sInd1_, sInd2_, sInd3_

        sInd1_ = S13
        sInd2_ = S23
        sInd3_ = S33

    end subroutine get_sIndW


    subroutine get_eIndW( eInd1_, eInd2_, eInd3_ ) bind(c,name='VS_get_eIndW')
        implicit none
        INTEGER(c_int), intent(out)  :: eInd1_, eInd2_, eInd3_

        eInd1_ = N13
        eInd2_ = N23
        eInd3_ = N33

    end subroutine get_eIndW


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


    subroutine get_bl( bl1_, bl2_, bl3_ ) bind(c,name='SVS_get_bl')
        implicit none
        INTEGER(c_int), intent(out)  :: bl1_,bl2_,bl3_

        bl1_ = b1L
        bl2_ = b2L
        bl3_ = b3L

    END subroutine get_bl


    subroutine get_bu( bu1_, bu2_, bu3_ ) bind(c,name='SVS_get_bu')
        implicit none
        INTEGER(c_int), intent(out)  :: bu1_,bu2_,bu3_

        bu1_ = b1U
        bu2_ = b2U
        bu3_ = b3U

    END subroutine get_bu


    subroutine set_COMM_CART( comm ) bind( c, name='SG_setCommCart' )
        implicit none
        integer(c_int), intent(in) :: comm

        COMM_CART = comm

    end subroutine set_COMM_CART


    subroutine set_rank( rank_ ) bind( c, name='SG_setRank' )
        implicit none
        integer(c_int), intent(in) :: rank_

        rank = rank_

    end subroutine set_rank


    subroutine set_iB( iB_ ) bind( c, name='SG_setIB' )
        implicit none
        integer(c_int), intent(in) :: iB_(1:3)
        integer :: i

        do i = 1, 3
            iB(i,1) = iB_(i)
        end do

    end subroutine set_iB


    subroutine set_Shift( shift_ ) bind( c, name='SG_setShift' )
        implicit none
        integer(c_int), intent(in) :: shift_(1:3)

        iShift = shift_(1)
        jShift = shift_(2)
        kShift = shift_(3)

    end subroutine set_Shift


    subroutine set_rankLU( rankl, ranku ) bind( c, name='SG_setRankLU' )
        implicit none
        integer(c_int), intent(in) :: rankl(1:3)
        integer(c_int), intent(in) :: ranku(1:3)

        rank1L = rankl(1)
        rank2L = rankl(2)
        rank3L = rankl(3)

        rank1U = ranku(1)
        rank2U = ranku(2)
        rank3U = ranku(3)

    end subroutine set_rankLU


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

end module cmod_FieldSpace
