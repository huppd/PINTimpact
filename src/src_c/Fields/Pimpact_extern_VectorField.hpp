#pragma once
#ifndef PIMPACT_EXTERN_VECTORFIELD_HPP
#define PIMPACT_EXTERN_VECTORFIELD_HPP




namespace Pimpact {


extern "C" {


void VF_write( double* phiU, double* phiV, double* phiW, const int& count );


void VF_init_2DPulsatileXC(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L2,
    const double* const x2,
    const double& re, const double& om, const double& px,
    double* const phiU, double* const phiV, double* const phiW );


void VF_init_2DPulsatileYC(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L1,
    const double* const x1,
    const double& re, const double& om, const double& px,
    double* const phiU, double* const phiV, double* const phiW );


void VF_init_2DPulsatileXS(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L2,
    const double* const x2,
    const double& re, const double& om, const double& px,
    double* const phiU, double* const phiV, double* const phiW );


void VF_init_2DPulsatileYS(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L1,
    const double* const x1,
    const double& re, const double& om, const double& px,
    double* const phiU, double* const phiV, double* const phiW );


void VF_init_StreamingC(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L1,
    const double* const x1,
    const double& amp,
    double* const phiU, double* const phiV, double* const phiW );


void VF_init_StreamingS(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L1,
    const double* const x1,
    const double& amp,
    double* const phiU, double* const phiV, double* const phiW );


void VF_init_Circle(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const L,
    const double* const x1,
    const double* const x2,
    double* const phiU,
    double* const phiV,
    double* const phiW );


void VF_init_RankineVortex(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const L,
    const double* const x1p,
    const double* const x2p,
    const double* const x1u,
    const double* const x2v,
    double* const phiU, double* const phiV, double* const phiW );

void VF_init_GaussianForcing1D(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L1,
    const double* const x1u,
    double* const phiU, double* const phiV, double* const phiW );

void VF_init_BoundaryFilter1D(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L1,
    const double* const x1u,
    double* const phiU, double* const phiV, double* const phiW );

void VF_init_GaussianForcing2D(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const L,
    const double* const x1p,
    const double* const x2p,
    const double* const x1u,
    const double* const x2v,
    double* const phiU, double* const phiV, double* const phiW );

void VF_init_BoundaryFilter2D(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const L,
    const double* const x1p,
    const double* const x2p,
    const double* const x1u,
    const double* const x2v,
    double* const phiU, double* const phiV, double* const phiW );

void VF_init_Vpoint(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const L,
    const double* const x1u,
    const double* const x2p,
    const double& sig,
    double* const phiU, double* const phiV, double* const phiW );

void VF_init_Disc(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const x1p,
    const double* const x2p,
    const double* const x3p,
    const double* const x1u,
    const double* const x2v,
    const double& xm, const double& ym, const double& rad,
    double* const phiU, double* const phiV, double* const phiW );

void VF_init_RotatingDisc(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const x1p,
    const double* const x2p,
    const double& xm, const double& ym, const double& omega,
    double* const phiU, double* const phiV, double* const phiW );

} // end of extern 'C'



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EXTERN_VECTORFIELD_HPP
