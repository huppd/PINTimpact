#pragma once
#ifndef PIMPACT_EXTERN_VECTORFIELD_HPP
#define PIMPACT_EXTERN_VECTORFIELD_HPP




namespace Pimpact {


extern "C" {

  extern "C" void OP_extrapolateBC2(
    const int& m,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int& dL,
    const int& dU,
    const int& BC_L,
    const int& BC_U,
    const int* const SB,
    const int* const NB,
    const double* const c,
    const double*       phi );

  void VF_write( double* phiU, double* phiV, double* phiW, const int& count );

  void VF_init_2DPulsatileX(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double& L,
    const double& t,
    const double* const x2,
    const double& re, const double& alpha, const double& px,
    double* const phiU, double* const phiV, double* const phiW );

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
    const double& xm, const double& ym, const double& rad,const double& sca,
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

  void VF_init_SHBF(
    const int& rank,
    const int& iShift,
    const int& IB1,
    const int* const M,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const dL,
    const int* const dU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const double* const y1p,
    const double* const y1u,
    const double* const x3w,
    const double* const cIup,
    const double& Re,
    const int& nonDim,
    const double& kappa,
    const double& sweep_angle_degrees,
    double* const velU,
    double* const velV,
    double* const velW );


  void VF_init_Dist(
    const int& rank,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SU,
    const int* const NU,
    const int* const SV,
    const int* const NV,
    const int* const SW,
    const int* const NW,
    const int& BC_3L_global,
    const double* const x1u,
    const double* const x1p,
    const double* const x2p,
    const double* const x3w,
    const double* const x3p,
    const int& dist_type,
    const double& vortex_ampli_prim,
    const double& vortex_x1pos,
    const double& vortex_x3pos,
    const double& vortex_radius,
    const double& vortex_band,
    double* const velU,
    double* const velV,
    double* const velW );


} // end of extern 'C'


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EXTERN_VECTORFIELD_HPP
