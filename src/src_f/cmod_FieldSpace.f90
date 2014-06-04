!>  Modul: cmod_FieldSpace
!!
!! helps \c Pimpact::FieldSpace to extract varibales from impact
module cmod_FieldSpace


    use iso_c_binding
    !  use mpi


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


end module cmod_FieldSpace
