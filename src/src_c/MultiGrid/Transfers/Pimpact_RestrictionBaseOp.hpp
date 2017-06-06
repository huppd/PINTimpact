#pragma once
#ifndef PIMPACT_RESTRICTIONBASEOP_HPP
#define PIMPACT_RESTRICTIONBASEOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Space.hpp"




namespace Pimpact {


extern "C" {

  void MG_getCRVS(
    const int& N,
    const int& BC_L,
    const int& BC_U,
    const int& dd,
    const int& Nf,
    const int& bL,
    const int& bU,
    const double* const xf,
    double* const cR );

  void MG_getCRS(
    const int& N,
    const int& BC_L,
    const int& BC_U,
    const int& dd,
    const int& Nf,
    const int& bL,
    const int& bU,
    const double* const xf,
    double* const cR );

  void MG_restrictHW(
    const int dimens,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const iimax,
    const int* const dd,
    const double* const cR1,
    const double* const cR2,
    const double* const cR3,
    const double* const phif,
    double* const phic );

  void MG_restrictFW(
    const int& dimens,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const iimax,
    const int* const dd,
    const double* const cR1,
    const double* const cR2,
    const double* const cR3,
    const double* const phif,
    double* const phic );

  void MG_restrictHWV(
    const int dimens,
    const int& dir,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const SSf,
    const int* const NNf,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const SSc,
    const int* const NNc,
    const int* const iimax,
    const int* const dd,
    const double* const cRV,
    const double* const phif,
    double* const phic );

  void MG_restrictFWV(
    const int& dimens,
    const int& dir,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const SSf,
    const int* const NNf,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const SSc,
    const int* const NNc,
    const int* const iimax,
    const int* const dd,
    const double* const cRV,
    const double* const cR1,
    const double* const cR2,
    const double* const cR3,
    const double* const phif,
    double* const phic );

  void MG_RestrictGather(
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const iimax,
    const int* const n_gather,
    const bool& participate_yes,
    const int& rankc2,
    const int& comm2,
    const int* const recvR,
    const int* const dispR,
    const int* const sizsR,
    const int* const offsR,
    double* const phic );
}



/// \brief Operator that restricts from a fine space to a coarse space
///
/// \todo c++fy
/// \tparam ST type of the \c Space
template<class ST>
class RestrictionBaseOp {

  static const int dimension = ST::dimension;

  using Scalar = typename ST::Scalar;
  using Ordinal = typename ST::Ordinal;

public:

  using SpaceT = ST;

  using FSpaceT = SpaceT;
  using CSpaceT = SpaceT;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  Teuchos::RCP<const SpaceT> spaceF_; 		  ///< fine space
  Teuchos::RCP<const SpaceT> spaceC_; 		  ///< coarse space

  int rankc2_; 														  ///< rank on coarse grid that gathers

  MPI_Comm comm2_; 												  ///< gather communicator

  Teuchos::Tuple<int,dimension> nGather_;   ///< number of processor that are gathered (npf/npc)

  Teuchos::Tuple<Ordinal,dimension> iimax_; ///< coarsen N + neighbor stuff (NC-1)/nGather + 1
  Teuchos::Tuple<Ordinal,dimension> dd_;    ///< coarsening factor

  Ordinal* offsR_;
  Ordinal* sizsR_;
  Ordinal* recvR_;
  Ordinal* dispR_;


  void init( const Teuchos::Tuple<int,dimension>& np ) {

    // ------------- nGather_, iimax_
    Teuchos::Tuple<int,dimension> periodic = spaceF_->getBCGlobal()->periodic();

    const Teuchos::Tuple<int,dimension>& npF = spaceF_->getProcGrid()->getNP();
    const Teuchos::Tuple<int,dimension>& npC = spaceC_->getProcGrid()->getNP();

    Teuchos::Tuple<Ordinal,dimension> iiShift;

    //		bool gather_yes = false;
    int nGatherTotal = 1;
    for( int i=0; i<dimension; ++i ) {
      nGather_[i] = npF[i] / npC[i];
      nGatherTotal *= nGather_[i];

      if( i<3 ) {
        iimax_[i] = (spaceC_->nLoc(i) - 1)/nGather_[i] + 1;
        dd_[i] = std::max( (spaceF_->nLoc(i) - 1)/( iimax_[i] -1 ), static_cast<Ordinal>(1) ); // check

        if( spaceF_->getStencilWidths()->getLS(i)==0 )
          if( spaceF_->getBCGlobal()->getBCL(i)==BC::Neighbor || spaceF_->getBCGlobal()->getBCL(i)==BC::Periodic )
            iimax_[i] = iimax_[i]-1;
        if( spaceF_->getStencilWidths()->getLS(i)==-1 )
          if( spaceF_->getBCGlobal()->getBCU(i)==BC::Neighbor || spaceF_->getBCGlobal()->getBCU(i)==BC::Periodic )
            iimax_[i] = iimax_[i]-1;
      } else {
        iimax_[i] = spaceC_->nLoc(i);
        dd_[i]    = 1;
      }

      iiShift[i] = ( iimax_[i] - 1 )*( ( spaceF_->getProcGrid()->getIB(i) -1 )%nGather_[i] );
      if( i<3 )
        if( spaceF_->getStencilWidths()->getLS(i)==0 )
          if( spaceF_->getBCGlobal()->getBCL(i)==BC::Neighbor || spaceF_->getBCGlobal()->getBCL(i)==BC::Periodic )
            iiShift[i] = iiShift[i] - 1;
    }

    offsR_ = new Ordinal[3*nGatherTotal];
    sizsR_ = new Ordinal[3*nGatherTotal];
    recvR_ = new Ordinal[  nGatherTotal];
    dispR_ = new Ordinal[  nGatherTotal];

    // ------------- rank2_, comm2_
    if( nGatherTotal>1 ) {
      int * newRanks = new int[nGatherTotal];

      MPI_Comm commWorld = spaceF_->getProcGrid()->getCommWorld();
      MPI_Comm commTemp;
      MPI_Group baseGroup, newGroup;
      MPI_Comm_group( commWorld, &baseGroup );

      Teuchos::Tuple<int,dimension> coord;

      if( 3==dimension ) {
        for( int kk=0; kk<npC[2]; ++kk )
          for( int jj=0; jj<npC[1]; ++jj )
            for( int ii=0; ii<npC[0]; ++ii ) {
              bool member_yes = false;

              for( int k=0; k<nGather_[2]; ++k ) {
                coord[2] = ((kk*nGather_[2]+k)*np[2]/npF[2])%np[2];
                for( int j=0; j<nGather_[1]; ++j ) {
                  coord[1] = ((jj*nGather_[1]+j)*np[1]/npF[1])%np[1];
                  for( int i=0; i<nGather_[0]; ++i ) {
                    coord[0] = ((ii*nGather_[0]+i)*np[0]/npF[0])%np[0];

                    MPI_Cart_rank(
                      commWorld,
                      coord.getRawPtr(),
                      &newRanks[i+nGather_[0]*j+nGather_[0]*nGather_[1]*k] );

                    if( newRanks[i+nGather_[0]*j+nGather_[0]*nGather_[1]*k]==spaceF_->rankST() )
                      member_yes = true;

                  }
                }
              }

              MPI_Group_incl( baseGroup, nGatherTotal, newRanks, &newGroup );
              MPI_Comm_create( commWorld, newGroup, &commTemp );
              MPI_Group_free( &newGroup );

              if( member_yes ) {
                MPI_Cart_create(
                  commTemp,							// input communicator
                  3,										// number of dimensions of Cartesian grid
                  nGather_.getRawPtr(), // integer array of size ndims specifying the number of processes in each dimension
                  periodic.getRawPtr(), // logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
                  false,                // ranking may be reordered (true) or not (false) (logical)
                  &comm2_ );            // communicator with new Cartesian topology (handle)
                MPI_Comm_free( &commTemp );
                int rank_comm2 = 0;
                if( spaceC_->getProcGrid()->participating() )
                  MPI_Comm_rank( comm2_, &rank_comm2 );
                MPI_Allreduce(
                  &rank_comm2, // starting address of send buffer (choice) starting
                  &rankc2_,    // starting address of receive buffer (choice)
                  1,           // number of elements in send buffer (non-negative inte- ger)
                  MPI_INTEGER, // data type of elements of send buffer (handle)
                  MPI_SUM,     // operation (handle)
                  comm2_ );    // communicator (handle)

              }
            }
      } else {
        for( int ll=0; ll<npC[3]; ++ll )
          for( int kk=0; kk<npC[2]; ++kk )
            for( int jj=0; jj<npC[1]; ++jj )
              for( int ii=0; ii<npC[0]; ++ii ) {

                bool member_yes = false;

                for( int l=0; l<nGather_[3]; ++l ) {
                  coord[3] = ((ll*nGather_[3]+l)*np[3]/npF[3])%np[3];
                  for( int k=0; k<nGather_[2]; ++k ) {
                    coord[2] = ((kk*nGather_[2]+k)*np[2]/npF[2])%np[2];
                    for( int j=0; j<nGather_[1]; ++j ) {
                      coord[1] = ((jj*nGather_[1]+j)*np[1]/npF[1])%np[1];
                      for( int i=0; i<nGather_[0]; ++i ) {
                        coord[0] = ((ii*nGather_[0]+i)*np[0]/npF[0])%np[0];

                        MPI_Cart_rank(
                          commWorld,
                          coord.getRawPtr(),
                          &newRanks[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] + l*nGather_[0]*nGather_[1]*nGather_[2] ]  );

                        if(  newRanks[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] + l*nGather_[0]*nGather_[1]*nGather_[2] ]==spaceF_->rankST() )
                          member_yes = true;

                      }
                    }
                  }
                }

                MPI_Group_incl( baseGroup, nGatherTotal, newRanks, &newGroup );
                MPI_Comm_create( commWorld, newGroup, &commTemp );
                MPI_Group_free( &newGroup );

                if( member_yes ) {
                  MPI_Cart_create(
                    commTemp,							// input communicator
                    4,										// number of dimensions of Cartesian grid
                    nGather_.getRawPtr(), // integer array of size ndims specifying the number of processes in each dimension
                    periodic.getRawPtr(), // logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
                    false,                // ranking may be reordered (true) or not (false) (logical)
                    &comm2_ );            // communicator with new Cartesian topology (handle)
                  MPI_Comm_free( &commTemp );

                  int temp[] = {1,1,1,0};
                  MPI_Cart_sub( comm2_, temp, &commTemp );
                  MPI_Comm_free( &comm2_ );
                  comm2_ = commTemp;

                  int rank_comm2 = 0;
                  if( spaceC_->getProcGrid()->participating() )
                    MPI_Comm_rank( comm2_, &rank_comm2 );
                  MPI_Allreduce(
                    &rank_comm2, // starting address of send buffer (choice) starting
                    &rankc2_,    // starting address of receive buffer (choice)
                    1,           // number of elements in send buffer (non-negative inte- ger)
                    MPI_INTEGER, // data type of elements of send buffer (handle)
                    MPI_SUM,     // operation (handle)
                    comm2_ );    // communicator (handle)

                }
              }
      }

      MPI_Group_free( &baseGroup );
      delete[] newRanks;
      // ------------------------- offsR_, sizsR_

      if( spaceF_->getProcGrid()->participating() )  {
        int rank_comm2;
        MPI_Comm_rank( comm2_, &rank_comm2 );

        std::vector<Ordinal> offs_global(3*nGatherTotal);
        std::vector<Ordinal> sizs_global(3*nGatherTotal);

        for( Ordinal i=0; i<3*nGatherTotal; ++i ) {
          offs_global[i] = 0;
          sizs_global[i] = 0;
        }


        for( Ordinal i=0; i<3; ++i ) {
          offs_global[ i + rank_comm2*3 ] = iiShift[i];
          sizs_global[ i + rank_comm2*3 ] = iimax_[i];
        }


        MPI_Allreduce(
          offs_global.data(),
          offsR_,
          3*nGatherTotal,
          MPI_INTEGER,
          MPI_SUM,
          comm2_ );
        MPI_Allreduce(
          sizs_global.data(),
          sizsR_,
          3*nGatherTotal,
          MPI_INTEGER,
          MPI_SUM,
          comm2_ );

        Ordinal counter = 0;
        for( int k=0; k<nGather_[2]; ++k )
          for( int j=0; j<nGather_[1]; ++j )
            for( int i=0; i<nGather_[0]; ++i ) {
              recvR_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ]
                = sizsR_[ 0 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
                  * sizsR_[ 1 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
                  * sizsR_[ 2 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ];

              dispR_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ] = counter;
              counter += recvR_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ];
            }
      }
    }
  }

public:

  RestrictionBaseOp(
    const Teuchos::RCP<const SpaceT>& spaceF,
    const Teuchos::RCP<const SpaceT>& spaceC ):
    spaceF_(spaceF),
    spaceC_(spaceC),
    comm2_(MPI_COMM_NULL) {

    init( spaceF_->getProcGrid()->getNP() );
  }


  RestrictionBaseOp(
    const Teuchos::RCP<const SpaceT>& spaceF,
    const Teuchos::RCP<const SpaceT>& spaceC,
    const Teuchos::Tuple<int,dimension>& np ):
    spaceF_(spaceF),
    spaceC_(spaceC),
    comm2_(MPI_COMM_NULL) {

    init( np );
  }


  ~RestrictionBaseOp() {
    delete[] offsR_;
    delete[] sizsR_;
    delete[] recvR_;
    delete[] dispR_;
  }


  void gather( Scalar* y ) const {

    if( nGather_[0]*nGather_[1]*nGather_[2]>1 )
      MG_RestrictGather(
        spaceC_->nLoc(),
        spaceC_->bl(),
        spaceC_->bu(),
        iimax_.getRawPtr(),
        nGather_.getRawPtr(),
        spaceC_->getProcGrid()->participating(),
        rankc2_,
        MPI_Comm_c2f(comm2_),
        recvR_,
        dispR_,
        sizsR_,
        offsR_,
        y );
  }

}; // end of class RestrictionBaseOp



/// \todo colect all create methods in one file
template<template<class> class OpT, class SpaceT>
Teuchos::RCP<const OpT<SpaceT> > create(
  const Teuchos::RCP<const SpaceT>& spaceF,
  const Teuchos::RCP<const SpaceT>& spaceC ) {

  return( Teuchos::rcp( new OpT<SpaceT>(spaceF,spaceC) ) );
}


} // end of namespace Pimpact




#endif // end of #ifndef PIMPACT_RESTRICTIONBASEOP_HPP
