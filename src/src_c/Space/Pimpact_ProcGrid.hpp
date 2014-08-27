#pragma once
#ifndef PIMPACT_PROCGRID_HPP
#define PIMPACT_PROCGRID_HPP

#include <ostream>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_ProcGridSize.hpp"


extern "C" {
void SVS_get_comm(MPI_Fint&);

void SG_setCommCart( const int& comm );
void SG_setRank( const int& rank_ );
void SG_setIB( const int* const iB_);
void SG_setShift( const int* const shift);
void SG_setRankLU( const int* const rankl, const int* const ranku );
void SG_getRankLU(       int* const rankl,       int* const ranku );
void SG_setCommSlice( const int& slice1, const int& slice2, const int& slice3 );
void SG_setCommBar( const int& bar1, const int& bar2, const int& bar3 );
void SG_setRankSliceBar( const int* const rankSlice, const int* const rankBar );
}


namespace Pimpact{

/// \brief ProcGrid, needs ProcGridSize, globalBoundaryConditions, GridSizeLocal
/// and sets Impact accordingly
/// \tparam dim default 3 or 4
template< class Ordinal=int, int dim=3 >
class ProcGrid {

  template< class OT, int dT >
  friend Teuchos::RCP< ProcGrid<OT,dT> > createProcGrid();

  template< class OT, int dT >
  friend Teuchos::RCP< ProcGrid<OT,dT> > createProcGrid(
      const Teuchos::RCP< GridSizeLocal<OT,dT> >& gsl,
      const Teuchos::RCP< BoundaryConditionsGlobal>& bcg,
      const Teuchos::RCP< ProcGridSize<OT,dT> >& procGridSize );

public:

  MPI_Fint commSpaceTimef_;
  MPI_Comm commSpaceTime_;

  MPI_Fint commSpacef_;
  MPI_Comm commSpace_;

  int rankS_;
  int rankST_;

  /// processor coordinates going from 1..np
  Teuchos::Tuple<int,dim> iB_;
  /// index offset going for dim 1:3 strange fortran style dim==4 c style
  Teuchos::Tuple<Ordinal,dim> shift_;

protected:

  /// rank of lower neighbour
  Teuchos::Tuple<int,dim> rankL_;

  /// rank of upper neighbour
  Teuchos::Tuple<int,dim> rankU_;

public:

  MPI_Comm commSlice_[3] ;
  MPI_Comm commBar_[3];

  Teuchos::Tuple<int,3> rankSlice_;
  Teuchos::Tuple<int,3> rankBar_;


protected:

  ProcGrid(
      MPI_Fint commf,
      MPI_Comm comm,
      int rank,
      Teuchos::Tuple<int,3> rankL,
      Teuchos::Tuple<int,3> rankU ):
        commSpaceTimef_(commf),
        commSpaceTime_(comm),
        commSpacef_(commf),
        commSpace_(comm),
        rankS_(rank),
        rankST_(rank),
        iB_(),
        rankL_(rankL),
        rankU_(rankU),
        commSlice_(),
        commBar_(),
        rankSlice_(),
        rankBar_() {}

  ProcGrid(
      const Teuchos::RCP< GridSizeLocal<Ordinal,dim> >& gridSizeLocal,
      const Teuchos::RCP< BoundaryConditionsGlobal >& bcg,
      const Teuchos::RCP< ProcGridSize<Ordinal,dim> >& procGridSize ):
        iB_(),rankL_(),rankU_(),commSlice_(),commBar_(),rankSlice_(),rankBar_() {

    TEUCHOS_TEST_FOR_EXCEPTION(
        3!=dim && 4!=dim ,
        std::logic_error,
        "Pimpact::ProcGrid::cotor: dim:"<<dim<<" not implemented!!!\n");

    Teuchos::Tuple<int,dim> ijkB;      // mpi grid coordinates
    Teuchos::Tuple<int,dim> periodic;  // array for mpir to signal where periodic grid is, from manual should be bool

    // init periodic
    for( int i=0; i<3; ++i ) {
      if( bcg->getBCL(i)==PeriodicBC )
        periodic[i] = 1;
      else
        periodic[i] = 0;
    }
    if( 4==dim ) periodic[3] = 1;


    // true means ranking may be reorderd
    // comm_cart comm with cartesian grid informations
    //    int* bla = procGridSize->getRawPtr();
    //    int mpierror =
    MPI_Cart_create(
        MPI_COMM_WORLD,
        dim,
        procGridSize->getRawPtr(),
        periodic.getRawPtr(),
        true,
        &commSpaceTime_ );

    commSpaceTimef_ = MPI_Comm_c2f( commSpaceTime_ );

    if( 3==dim ) {

      commSpace_ = commSpaceTime_;

      commSpacef_ = MPI_Comm_c2f( commSpace_ );

      // gets rank from COMM_CART
      MPI_Comm_rank( commSpace_, &rankS_ );
      rankST_ = rankS_;

    }
    else if( 4==dim ) {

      {
        int temp[] = {1,1,1,0};
        MPI_Cart_sub( commSpaceTime_, temp, &commSpace_ );
      }

      commSpacef_ = MPI_Comm_c2f( commSpace_ );

      // gets rank from COMM_CART
      MPI_Comm_rank( commSpace_, &rankS_ );
      MPI_Comm_rank( commSpaceTime_, &rankST_ );
      //      std::cout << "rank in commSpace: " << rankS_ << "\trank in commSpaceTime: " << rankST_ << "\n";
    }

    // gets coordinates in xyz direction from rank and comm_cart
    MPI_Cart_coords(
        commSpaceTime_,
        rankST_,
        dim,
        ijkB.getRawPtr() );

    //    std::cout << "rankST: " << rankST_ << " coord: " << ijkB << "\n";

    //    MPI_Cart_coords(
    //        commSpace_,
    //        rankS_,
    //        dim,
    //        ijkB.getRawPtr() );
    //    std::cout << "rankS: " << rankS_ << " coord: " << ijkB << "\n";

    // stores coordinates in a fortran fasion?
    for( int i=0; i<dim; ++i )
      iB_[i] = ijkB[i] + 1;
    // computes index offset
    for( int i=0; i<3; ++i )
      shift_[i] = (iB_[i]-1)*( gridSizeLocal->get(i)-1 );
    if( 4==dim )
      shift_[3] = (iB_[3]-1)*( gridSizeLocal->get(3) );


    for( int i = 0; i<3; ++i )
      MPI_Cart_shift( commSpace_, i, 1, &rankL_[i], &rankU_[i] );
    if( 4==dim )
      MPI_Cart_shift( commSpaceTime_, 3, 1, &rankL_[3], &rankU_[3] );
    //        !                             ^ ^   ^      ^
    //        !                             | |   |      |
    //        !                             d d   r      r
    //        !                             i i   a      a
    //        !                             r s   n      n
    //        !                             e p   k      k
    //        !                             c l   s      d
    //        !                             t a   o      e
    //        !                             i c   u      s
    //        !                             o m   r      t
    //        !                             n e   c
    //        !                               n   e
    //        !                               t


    {
      int temp[] = {0,1,1};
      MPI_Cart_sub( commSpace_, temp, &commSlice_[0] );
    }{
      int temp[] = {1,0,1};
      MPI_Cart_sub( commSpace_, temp, &commSlice_[1] );
    }{
      int temp[] = {1,1,0};
      MPI_Cart_sub( commSpace_, temp, &commSlice_[2] );
    }

    {
      int temp[] = {1,0,0};
      MPI_Cart_sub( commSpace_, temp, &commBar_[0] );
    }{
      int temp[] = {0,1,0};
      MPI_Cart_sub( commSpace_, temp, &commBar_[1] );
    }{
      int temp[] = {0,0,1};
      MPI_Cart_sub( commSpace_, temp, &commBar_[2] );
    }

    for( int i=0; i<3; ++i ) {
      MPI_Comm_rank( commSlice_[i], &rankSlice_[i] );

      MPI_Comm_rank( commBar_[i], &rankBar_[i] );
    }


    // maybe dispensable...
    MPI_Errhandler_set(commSpaceTime_, MPI_ERRORS_ARE_FATAL );
    MPI_Errhandler_set(commSpace_, MPI_ERRORS_ARE_FATAL );
    for( int i=0; i<3; ++i ) {
      MPI_Errhandler_set( commSlice_[i], MPI_ERRORS_ARE_FATAL );
      MPI_Errhandler_set( commBar_[i],   MPI_ERRORS_ARE_FATAL );
    }

  }


public:

  void set_Impact() {
    SG_setCommCart( commSpacef_ );
    SG_setRank( rankS_ );
    SG_setIB( iB_.getRawPtr() );
    SG_setShift( shift_.getRawPtr() );
    SG_setRankLU( rankL_.getRawPtr(), rankU_.getRawPtr() );
    SG_setCommSlice(
        MPI_Comm_c2f(commSlice_[0]),
        MPI_Comm_c2f(commSlice_[1]),
        MPI_Comm_c2f(commSlice_[2]) );
    SG_setCommBar(
        MPI_Comm_c2f(commBar_[0]),
        MPI_Comm_c2f(commBar_[1]),
        MPI_Comm_c2f(commBar_[2]) );

    SG_setRankSliceBar( rankSlice_.getRawPtr(), rankBar_.getRawPtr() );
  }

  void print( std::ostream& out=std::cout ) const {
    out << "\t---ProcessorGrid: ---\n";
    out << "\trankS: " <<rankS_<<"\n";
    out << "\trankST: " <<rankST_<<"\n";
    out << "\trankL: " <<rankL_<<"\n";
    out << "\trankU: " <<rankU_<<"\n";
    out << "\tproc coordinate: " << iB_ << "\n";
    out << "\toffset: " << shift_ << "\n";
  }

  const int& getRankL( int i ) const { return( rankL_[i] ); }
  const int& getRankU( int i ) const { return( rankU_[i] ); }

  const int* getRankL() const { return( rankL_.getRawPtr() ); }
  const int* getRankU() const { return( rankU_.getRawPtr() ); }

};

/// \relates ProcGrid
template< class O=int, int d=3 >
Teuchos::RCP< ProcGrid<O,d> > createProcGrid() {


  MPI_Fint comm;
  SVS_get_comm( comm );

  int rank;
  SG_setRank( rank );

  Teuchos::Tuple<int,3> rankL;
  Teuchos::Tuple<int,3> rankU;
  SG_getRankLU( rankL.getRawPtr(), rankU.getRawPtr() );

  return(
      Teuchos::rcp( new ProcGrid<O,d>(
          comm,
          MPI_Comm_f2c(comm),
          rank,
          rankL,
          rankU ) ) );

}

/// \relates ProcGrid
template< class O=int, int d=3 >
Teuchos::RCP< ProcGrid<O,d> > createProcGrid(
    const Teuchos::RCP< GridSizeLocal<O,d> >& gsl,
    const Teuchos::RCP< BoundaryConditionsGlobal>& bcg,
    const Teuchos::RCP< ProcGridSize<O,d> >& procGridSize ) {

  return(
      Teuchos::rcp( new ProcGrid<O,d>( gsl, bcg, procGridSize ) ) );

}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PROCGRID_HPP
