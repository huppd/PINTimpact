#pragma once
#ifndef PIMPACT_EXTERN_SCALARFIELD_HPP
#define PIMPACT_EXTERN_SCALARFIELD_HPP




namespace Pimpact {


extern "C" {



void write_hdf5_2D(
        const int& rank,
        const int& COMM_CART,
        const int* const M,
        const int* const BC_L_global,
        const int* const BC_U_global,
        const int* const N,
        const int* const bL,
        const int* const bU,
        const int* const SS,
        const int* const NN,
        const int* const ls,
        const int* const NB,
        const int* const iB,
        const int* const iShift,
        const int& ftype,
        const int& filecount,
        const int& filenamelen,
        const double* const phi,
        const double* const y1p,
        const double* const y2p,
//        const double* const y3p,
        const double& Re,
        const double& alpha2 );


void write_hdf_3D(
    const int& rank,
    const int& COMM_CART,
    const int* const M,
    const int* const BC_L_global,
    const int* const BC_U_global,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const int* const ls,         
    const int* const NB,
    const int* const iB,
    const int* const iShift,
    const int& dir_name,
    const int& vel_dir,
    const int& filecount,
    const int& namelen,
    const int* const stride,
    const double* const phi,
    const double* const y1p,
    const double* const y2p,
    const double* const y3p,
    const double* const y1u,
    const double* const y2v,
    const double* const y3w,
    const double& Re,
    const double& alpha2 );



void F_exchange(
    const int& dim,
    const int& comm,
    const int* const rankL,
    const int* const rankU,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const bcL,
    const int* const bcU,
    const int* const Sp,
    const int* const Np,
    const int* const SS,
    const int* const NN,
    const int& dir, const int& vel_dir,
    double* const phi );


void SF_init_2DPoiseuilleX(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const double& L2,
    const double* const x2,
          double* const phi );


void SF_init_2DPoiseuilleY(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const double& L1,
    const double* const x1,
          double* const phi );


void SF_init_2DPoiseuilleZ(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const double& L1,
    const double* const x1,
          double* const phi );


void SF_init_2DGradX(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const double& L,
    const double* const x1,
          double* const phi,
	 const double& alpha );


void SF_init_2DGradY(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const double& L,
    const double* const x2,
          double* const phi,
	 const double& alpha );


void SF_init_2DGradZ(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const double& L,
    const double* const x2,
          double* const phi,
	 const double& alpha );


void SF_init_Vpoint(
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const SS,
    const int* const NN,
    const double* const x1,
    const double* const x2,
    const double* const x3,
    const double* const xc,
		const double& amp,
    const double* const sig,
		double* const phi );


} // end of extern 'C'


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EXTERN_SCALARFIELD_HPP
