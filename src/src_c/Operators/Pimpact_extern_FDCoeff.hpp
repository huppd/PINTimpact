#pragma once
#ifndef PIMPACT_EXTERN_FDCOEFF_HPP
#define PIMPACT_EXTERN_FDCOEFF_HPP



namespace Pimpact {


extern "C" {


void FD_getDiffCoeff(
    const int& rank,
    const int& Nmax,
    const int& bL,
    const int& bU,
    const int& cL,
    const int& cU,
    const int& BCL,
    const int& BCU,
    const int& SShift,
    const int& grid_type,
    const int& dir,
    const int& abl,
    const int& upwind,
    const bool& mapping_yes,
    const int& dim_ncb,
    const int* const n_coeff_bound,
    const double* const xC,
    const double* const xE,
          double* const cc );



} // end of extern 'C'



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EXTERN_FDCOEFF_HPP
