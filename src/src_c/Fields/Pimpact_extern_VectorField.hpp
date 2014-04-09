#pragma once
#ifndef PIMPACT_EXTERN_VECTORFIELD_HPP
#define PIMPACT_EXTERN_VECTORFIELD_HPP

//#include <iostream>
//#include "mpi.h"
//
//#include "Teuchos_RCP.hpp"
//#include "BelosTypes.hpp"
//#include "Teuchos_ScalarTraitsDecl.hpp"
//#include "Teuchos_SerialDenseMatrix.hpp"
//
//#include "Pimpact_Types.hpp"
//#include "Pimpact_FieldSpace.hpp"
//#include "Pimpact_IndexSpace.hpp"



namespace Pimpact {

extern "C" {


void VF_compNorm(
    const MPI_Fint& comm, const int& dimens,
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    const double* const phi1, const double* const phi2, const double* const phi3,
    const bool& inf_yes, const bool& two_yes,
    double& normInf, double& normTwo );


void VF_weightedNorm(
    const MPI_Fint& comm, const int& dimens,
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    const double* const phi1, const double* const phi2, const double* const phi3,
    const double* const wei1, const double* const wei2, const double* const wei3,
    double& norm );


void VF_dot(
      const MPI_Fint& comm,
      const int& dimens,
      const int& N1,  const int& N2,  const int& N3,
      const int& S1U, const int& S2U, const int& S3U,
      const int& N1U, const int& N2U, const int& N3U,
      const int& S1V, const int& S2V, const int& S3V,
      const int& N1V, const int& N2V, const int& N3V,
      const int& S1W, const int& S2W, const int& S3W,
      const int& N1W, const int& N2W, const int& N3W,
      const int& b1L, const int& b2L, const int& b3L,
      const int& b1U, const int& b2U, const int& b3U,
      const double* const phi1U, const double* const phi1V, const double* const phi1W,
      const double* const phi2U, const double* const phi2V, const double* const phi2W,
      double& scalar);


void F_exchange(
    const int& dir, const int& vel_dir,
    const int& SS1, const int& SS2, const int& SS3,
    const int& NN1, const int& NN2, const int& NN3,
    double* phi );


void VF_write( double* phiU, double* phiV, double* phiW, const int& count );


void VF_init_Zero(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_2DPoiseuilleX(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_2DPoiseuilleY(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_2DPulsatileXC(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_2DPulsatileYC(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_2DPulsatileXS(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_2DPulsatileYS(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW, const double& re, const double& om, const double& px );


void VF_init_Streaming(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_Circle(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );


void VF_init_RankineVortex(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );

void VF_init_GaussianForcing1D(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );

void VF_init_BoundaryFilter1D(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );

void VF_init_GaussianForcing2D(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );

void VF_init_BoundaryFilter2D(
    const int& N1,  const int& N2,  const int& N3,
    const int& S1U, const int& S2U, const int& S3U,
    const int& N1U, const int& N2U, const int& N3U,
    const int& S1V, const int& S2V, const int& S3V,
    const int& N1V, const int& N2V, const int& N3V,
    const int& S1W, const int& S2W, const int& S3W,
    const int& N1W, const int& N2W, const int& N3W,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phiU, double* phiV, double* phiW );

} // end of extern 'C'



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EXTERN_VECTORFIELD_HPP