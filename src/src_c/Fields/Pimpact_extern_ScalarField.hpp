#pragma once
#ifndef PIMPACT_EXTERN_SCALARFIELD_HPP
#define PIMPACT_EXTERN_SCALARFIELD_HPP



namespace Pimpact {


extern "C" {


void SF_add(
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  double* phi, const double* const  phi1, const double* const  phi2,
  const double& scalar1, const double& scalar2);


void SF_abs(
    const int& N1,  const int& N2,  const int& N3,
    const int& SS1, const int& SS2, const int& SS3,
    const int& NN1, const int& NN2, const int& NN3,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phi, const double* const  phi1 );


void SF_reciprocal(
    const int& N1,  const int& N2,  const int& N3,
    const int& SS1, const int& SS2, const int& SS3,
    const int& NN1, const int& NN2, const int& NN3,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phi, const double* const  phi1 );


void SF_compNorm(const MPI_Fint& comm,
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  double* phi,
  const bool& inf_yes, const bool& two_yes,
  double& normInf, double& normTwo);


void SF_weightedNorm(
  const MPI_Fint& comm,
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  double* phi,
  const double* const weights,
  double& norm);


void SF_dot( const MPI_Fint& comm,
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  const double* const phi1, const double* const phi2, double& scalar);


void SF_scale(
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  double* phi, const double& scalar );


void SF_scale2(
    const int& N1,  const int& N2,  const int& N3,
    const int& SS1, const int& SS2, const int& SS3,
    const int& NN1, const int& NN2, const int& NN3,
    const int& b1L, const int& b2L, const int& b3L,
    const int& b1U, const int& b2U, const int& b3U,
    double* phi, const double* const  phi1 );


void SF_random(
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  double* phi );


void SF_init(
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  double* phi, const double& scalar );


void SF_print(
  const int& N1,  const int& N2,  const int& N3,
  const int& SS1, const int& SS2, const int& SS3,
  const int& NN1, const int& NN2, const int& NN3,
  const int& b1L, const int& b2L, const int& b3L,
  const int& b1U, const int& b2U, const int& b3U,
  const double* phi );


void SF_write( double* phi, const int& count );


} // end of extern 'C'



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EXTERN_SCALARFIELD_HPP
