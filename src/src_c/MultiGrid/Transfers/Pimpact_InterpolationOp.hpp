#pragma once
#ifndef PIMPACT_INTERPOLATIONOP_HPP
#define PIMPACT_INTERPOLATIONOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_Stencil.hpp"




namespace Pimpact {



extern "C" {

  void MG_getCIS(
    const int& Nc,
    const int& Nf,
    const int& bL,
    const int& bU,
    const double* const xs,
    const int& dd,
    double* const cI );

  void MG_getCIV(
    const int& Nc,
    const int& bLc,
    const int& bUc,
    const int& SSc,
    const int& NNc,
    const int& BC_L,
    const int& BC_U,
    const int& Nf,
    const int& bLf,
    const int& bUf,
    const int& SSf,
    const double* const xc,
    const double* const xf,
    const int& dd,
    double* const cIV );

  void MG_InterpolateScatter(
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const np,
    const int* const iimax,
    const int* const n_gather,
    const bool& participating,
    const int& rank2,
    const int& comm2,
    const int* const dispI,
    const int* const offsI,
    const double* const phic );

  void MG_interpolate(
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const iimax,
    const int* const dd,
    const double* const cI1,
    const double* const cI2,
    const double* const cI3,
    const double* const phic,
    double* const phif );

  void MG_interpolateV(
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
    const double* const cIV,
    const double* const cI1,
    const double* const cI2,
    const double* const cI3,
    const double* const phif,
    double* const phic );

}



/// \brief Opetartor that interpolates from a coarse space to a fine space
///
/// \todo c++fy
/// \tparam ST type of the \c Space
template<class ST>
class InterpolationOp {

  static const int dimension = ST::dimension;

  using Scalar = typename ST::Scalar;
  using Ordinal = typename ST::Ordinal;

public:

  using SpaceT = ST;

  using FSpaceT = SpaceT;
  using CSpaceT = SpaceT;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

  using StencS = Stencil< Scalar, Ordinal, 1, 1, 2 >;
  using StencV = Stencil< Scalar, Ordinal, 0, 1, 2 >;

protected:

  Teuchos::RCP<const SpaceT> spaceC_; 		///< coarse space
  Teuchos::RCP<const SpaceT> spaceF_; 		///< fine space

  int rankc2_; 														///< rank on coarse grid that gathers

  MPI_Comm comm2_; 												///< gather communicator

  Teuchos::Tuple<int,dimension> nGather_; ///< number of proceesr that are gathered

  Teuchos::Tuple<Ordinal,dimension> iimax_;
  Teuchos::Tuple<Ordinal,dimension> dd_;

  Teuchos::ArrayRCP<Ordinal> offs_;
  Teuchos::ArrayRCP<Ordinal> sizs_;
  Teuchos::ArrayRCP<Ordinal> recv_;
  Teuchos::ArrayRCP<Ordinal> disp_;

  Teuchos::Tuple<StencS,3> cIS_;
  Teuchos::Tuple<StencV,3> cIV_;


  void init( const Teuchos::Tuple<int,dimension>& np ) {

    // ------------- nGather_, iimax_
    Teuchos::Tuple<int,dimension> periodic = spaceF_->getBCGlobal()->periodic();

    const Teuchos::Tuple<int,dimension>& npF = spaceF_->getProcGrid()->getNP();
    const Teuchos::Tuple<int,dimension>& npC = spaceC_->getProcGrid()->getNP();

    Teuchos::Tuple<Ordinal,dimension> iiShift;

    int nGatherTotal = 1;
    for( int i=0; i<dimension; ++i ) {
      nGather_[i] = npF[i] / npC[i]; // check
      nGatherTotal *= nGather_[i]; // check

      if( i<3 ) {
        iimax_[i] = (spaceC_->nLoc(i) - 1)/nGather_[i] + 1; // check
        dd_[i] = std::max( (spaceF_->nLoc(i) - 1)/( iimax_[i] -1 ), static_cast<Ordinal>(1) ); // check
        iiShift[i] = ( iimax_[i] - 1 )*( ( spaceF_->getProcGrid()->getIB(i) -1 )%nGather_[i] ); // check
      } else {
        iimax_[i] = spaceC_->nLoc(i); // check
        dd_[i] = 1; // check
        iiShift[i] = ( iimax_[i] - 1 )*( ( spaceF_->getProcGrid()->getIB(i) -1 )%nGather_[i] ); // check
      }
    }

    offs_ = Teuchos::arcp<Ordinal>( 3*nGatherTotal );
    sizs_ = Teuchos::arcp<Ordinal>( 3*nGatherTotal );
    recv_ = Teuchos::arcp<Ordinal>(   nGatherTotal );
    disp_ = Teuchos::arcp<Ordinal>(   nGatherTotal );

    MPI_Request req_o, req_s;  

    // ------------- rank2_, comm2_
    if( nGatherTotal>1 ) {
      Teuchos::ArrayRCP<int> newRanks = Teuchos::arcp<int>(nGatherTotal);

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

              MPI_Group_incl( baseGroup, nGatherTotal, newRanks.getRawPtr(), &newGroup );
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
                          &newRanks[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] + l*nGather_[0]*nGather_[1]*nGather_[2] ] );

                        if( newRanks[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] + l*nGather_[0]*nGather_[1]*nGather_[2] ]==spaceF_->rankST() )
                          member_yes = true;

                      }
                    }
                  }
                }

                MPI_Group_incl( baseGroup, nGatherTotal, newRanks.getRawPtr(), &newGroup );
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
      // ------------------------- offs_, sizs_

      if( spaceF_->getProcGrid()->participating() )  {
        int rank_comm2;
        MPI_Comm_rank( comm2_, &rank_comm2 );

        MPI_Iallgather(
            iiShift.getRawPtr(), //void* send_data,
            3,                   //int send_count,
            MPI_INTEGER,         //MPI_Datatype send_datatype,
            offs_.getRawPtr(),  //void* recv_data,
            3,                   //int recv_count,
            MPI_INTEGER,         //MPI_Datatype recv_datatype,
            comm2_,              //MPI_Comm communicator
            &req_o );              //MPI_Request 

        Teuchos::Tuple<Ordinal, 3> sizs_local;//  = Teuchos::arcp<Ordinal>(3);

        for( Ordinal i=0; i<3; ++i )
          sizs_local[i] = iimax_[i]+1;

        MPI_Iallgather(
            sizs_local.getRawPtr(), //void* send_data,
            3,                   //int send_count,
            MPI_INTEGER,         //MPI_Datatype send_datatype,
            sizs_.getRawPtr(),  //void* recv_data,
            3,                   //int recv_count,
            MPI_INTEGER,         //MPI_Datatype recv_datatype,
            comm2_,              //MPI_Comm communicator
            &req_s );              //MPI_Request 

        MPI_Wait(&req_s, MPI_STATUS_IGNORE); 

        Ordinal counter = 0;
        for( int k=0; k<nGather_[2]; ++k )
          for( int j=0; j<nGather_[1]; ++j )
            for( int i=0; i<nGather_[0]; ++i ) {
              recv_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ]
                = sizs_[ 0 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
                  * sizs_[ 1 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
                  * sizs_[ 2 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ];

              disp_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ] = counter;
              counter += recv_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ];
            }
        MPI_Wait(&req_o, MPI_STATUS_IGNORE); 
      }
    }
    // ------------------ cIS, cIV
    for( int dir=0; dir<3; ++dir ) {

      //if dd>1
      cIS_[dir] = StencS( spaceC_->nLoc(dir) );

      MG_getCIS(
        spaceC_->nLoc(dir),
        spaceF_->nLoc(dir),
        spaceC_->bl(dir),
        spaceC_->bu(dir),
        spaceF_->getCoordinatesLocal()->getX( F::S, dir  ),
        dd_[dir],
        cIS_[dir].get() );

      cIV_[dir] = StencV( spaceF_->nLoc(dir) );

      Ordinal offset = 0;
      if( 1!=nGather_[dir] )
        offset =
          ( iimax_[dir]-1 )*( spaceF_->getProcGrid()->getIB(dir)-1 - nGather_[dir]*(spaceC_->getProcGrid()->getIB(dir)-1) );

      F fdir = static_cast<F>( dir );
      MG_getCIV(
        spaceC_->nLoc(dir),
        spaceC_->bl(dir),
        spaceC_->bu(dir),
        spaceC_->sInd(fdir)[dir],
        spaceC_->eInd(fdir)[dir],
        spaceF_->getBCLocal()->getBCL(dir),
        spaceF_->getBCLocal()->getBCU(dir),
        spaceF_->nLoc(dir),
        spaceF_->bl(dir),
        spaceF_->bu(dir),
        offset,
        spaceC_->getCoordinatesLocal()->getX( fdir, dir ),
        spaceF_->getCoordinatesLocal()->getX( fdir, dir ),
        dd_[dir],
        cIV_[dir].get() );
    }

  } // end of void init( const Teuchos::Tuple<int,dimension>& np )


public:

  InterpolationOp( const Teuchos::RCP<const SpaceT>& spaceC, const Teuchos::RCP<const SpaceT>& spaceF ):
    spaceC_(spaceC),
    spaceF_(spaceF),
    comm2_(MPI_COMM_NULL) {

    init( spaceF_->getProcGrid()->getNP() );
  }

  InterpolationOp(
    const Teuchos::RCP<const SpaceT>& spaceC,
    const Teuchos::RCP<const SpaceT>& spaceF,
    const Teuchos::Tuple<int,dimension>& np ):
    spaceC_(spaceC),
    spaceF_(spaceF),
    comm2_(MPI_COMM_NULL) {

    init( np );
  }


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    y.init();

    F fType = x.getType();

    assert( x.getType()==y.getType() );

    if( F::S==fType ) {

      if( spaceC_->getProcGrid()->participating() )
        x.exchange();

      if( nGather_[0]*nGather_[1]*nGather_[2]>1 ) {
        MG_InterpolateScatter(
          spaceC_->nLoc(),
          spaceC_->bl(),
          spaceC_->bu(),
          spaceF_->np(),
          iimax_.getRawPtr(),
          nGather_.getRawPtr(),
          spaceC_->getProcGrid()->participating(),
          rankc2_,
          MPI_Comm_c2f(comm2_),
          disp_.getRawPtr(),
          offs_.getRawPtr(),
          x.getConstRawPtr() );
      }

      MG_interpolate(
        spaceC_->nLoc(),
        spaceC_->bl(),
        spaceC_->bu(),
        spaceF_->nLoc(),
        spaceF_->bl(),
        spaceF_->bu(),
        iimax_.getRawPtr(),
        dd_.getRawPtr(),
        cIS_[0].get(),
        cIS_[1].get(),
        cIS_[2].get(),
        x.getConstRawPtr(),
        y.getRawPtr() );

    } else {

      //		y.init();
      int dir = static_cast<int>( fType );

      if( spaceC_->getProcGrid()->participating() ) {
        switch( fType ) {
        case F::U:
          x.exchange(1);
          x.exchange(2);
          x.exchange(0);
          break;
        case F::V:
          x.exchange(2);
          x.exchange(0);
          x.exchange(1);
          break;
        case F::W:
          x.exchange(0);
          x.exchange(1);
          x.exchange(2);
        case F::S:
          break;
        }
      }

      if( nGather_[0]*nGather_[1]*nGather_[2]>1 ) {
        MG_InterpolateScatter(
          spaceC_->nLoc(),
          spaceC_->bl(),
          spaceC_->bu(),
          spaceF_->np(),
          iimax_.getRawPtr(),
          nGather_.getRawPtr(),
          spaceC_->getProcGrid()->participating(),
          rankc2_,
          MPI_Comm_c2f(comm2_),
          disp_.getRawPtr(),
          offs_.getRawPtr(),
          x.getConstRawPtr() );
      }

      MG_interpolateV(
        dir+1,
        spaceC_->nLoc(),
        spaceC_->bl(),
        spaceC_->bu(),
        spaceC_->sIndB(fType),
        spaceC_->eIndB(fType),
        spaceF_->nLoc(),
        spaceF_->bl(),
        spaceF_->bu(),
        spaceF_->sIndB(fType),
        spaceF_->eIndB(fType),
        iimax_.getRawPtr(),
        dd_.getRawPtr(),
        cIV_[dir].get(),
        cIS_[0].get(),
        cIS_[1].get(),
        cIS_[2].get(),
        x.getConstRawPtr(),
        y.getRawPtr() );

    }

    y.changed();
  }


  void print(  std::ostream& out=std::cout ) const {

    out << "=== Interpolation OP ===\n";
    out << "nGather:\t" << nGather_ << "\n";
    out << "dd:\t" << dd_ << "\n";
    out << "iimax:\t" << iimax_ << "\n";
    out << "rankc2:\t" << rankc2_ << "\n";
    out << "comm2:\t" << comm2_ << "\n";
    //		out << "offs:\t" << offs_ << "\n";
    //		out << "dispI:\t" << disp_ << "\n";
    out << "\n";
    for( int j=0; j<3; ++j ) {
      out << "\n Scalar dir: " << j << ":\n";
      cIS_[j].print( out );
    }
    out << "\n";

    for( int j=0; j<3; ++j ) {
      out << "\n Vector dir: " << j << ":\n";
      cIV_[j].print( out );
    }
    out << "\n";
  }

  Teuchos::RCP<const SpaceT> spaceC() const {
    return( spaceC_ );
  };
  Teuchos::RCP<const SpaceT> spaceF() const {
    return( spaceF_ );
  };

}; // end of class InterpolationOp


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_INTERPOLATIONOP_HPP
