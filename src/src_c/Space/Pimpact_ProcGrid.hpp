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
void SG_setCommCart( const int& comm );
void SG_setRank( const int& rank_ );
void SG_setIB( const int* const iB_);
void SG_setShift( const int* const shift);
void SG_setRankLU( const int* const rankl, const int* const ranku );
void SG_setCommSlice( const int& slice1, const int& slice2, const int& slice3 );
void SG_setCommBar( const int& bar1, const int& bar2, const int& bar3 );
void SG_setRankSliceBar( const int* const rankSlice, const int* const rankBar );


}


namespace Pimpact{

/// \brief ProcGrid
/// \tparam dim default 3 or 4
/// for creation BCGlobal/ ProcGridSize
template< class Ordinal, int dim=3 >
class ProcGrid {

  MPI_Fint commSpaceTimef_;
  MPI_Comm commSpaceTime_;

  MPI_Fint commSpacef_;
  MPI_Comm commSpace_;

  int rank_;

  Teuchos::Tuple<int,dim> iB_;
  Teuchos::Tuple<int,dim> shift_;

  Teuchos::Tuple<int,dim> rankL_;
  Teuchos::Tuple<int,dim> rankU_;

  MPI_Comm commSlice_[3] ;
  MPI_Comm commBar_[3];

  Teuchos::Tuple<int,3> rankSlice_;
  Teuchos::Tuple<int,3> rankBar_;

  //friend createProcGrid<Ordinal,dim>;

public:

  ProcGrid(
      const Teuchos::RCP< GridSizeLocal<Ordinal,dim> >& gridSizeLocal,
      const Teuchos::RCP< BoundaryConditionsGlobal >& bcg,
      const Teuchos::RCP< ProcGridSize<Ordinal,dim> >& procGridSize ):
        iB_(),rankL_(),rankU_(),commSlice_(),commBar_(),rankSlice_(),rankBar_() {

    if( 3!=dim && 4!=dim ) {
      std::cout<< "Pimpact::ProcGrid::cotor: dim:"<<dim<<" not implemented!!!\n";
      return;
    }

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

    //    // Macht keinen messbaren Unterschied (falls doch irgendwann, dann sollte der CALL auch auf die Grobgitter-Kommunikatoren auch angewandt werden!):
    //    MPI_Cart_create(
    //        MPI_COMM_WORLD,
    //        dim,
    //        procGridSize->getRawPtr(),
    //        periodic->getRawPtr(),
    //        false,
    //        commSpaceTime_)
    //    MPI_Cart_map(
    //        commSpaceTime_,
    //        dim,
    //        procGridSize->getRawPtr(),
    //        periodic->getRawPtr(),
    //        rank )

    // true means ranking may be reorderd
    // comm_cart comm with cartesian grid informations
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
      MPI_Comm_rank( commSpace_, &rank_ );

    }
    else if( 4==dim ) {

      {
        int temp[] = {1,1,1,0};
        MPI_Cart_sub( commSpaceTime_, temp, &commSpace_ );
      }

      commSpacef_ = MPI_Comm_c2f( commSpace_ );

      // gets rank from COMM_CART
      MPI_Comm_rank( commSpace_, &rank_ );
    }

    // gets coordinates in xyz direction from rank and comm_cart
    MPI_Cart_coords(
        commSpaceTime_,
        rank_,
        dim,
        ijkB.getRawPtr() );

    // stores coordinates in a fortran fasion?
    for( int i=0; i<dim; ++i ) {
      iB_[i] = ijkB[i] + 1;

      // computes index ofset
      shift_[i] = (iB_[i]-1)*( gridSizeLocal->get(i)-1 );

      MPI_Cart_shift( commSpaceTime_, i, 1, &rankL_[i], &rankU_[i] );
      //      !                             ^ ^   ^      ^
      //      !                             | |   |      |
      //      !                             d d   r      r
      //      !                             i i   a      a
      //      !                             r s   n      n
      //      !                             e p   k      k
      //      !                             c l   s      d
      //      !                             t a   o      e
      //      !                             i c   u      s
      //      !                             o m   r      t
      //      !                             n e   c
      //      !                               n   e
      //      !                               t

    }

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

    set_Impact();

  }

  void set_Impact() {
    SG_setCommCart( commSpacef_ );
    SG_setRank( rank_ );
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
    out << "\trank: " <<rank_<<"\n";
    out << "\trankL: " <<rankL_<<"\n";
    out << "\trankU: " <<rankU_<<"\n";
    out << "\tiB: " << iB_ << "\n";
  }


};


/// \relates ProcGrid
template< class O=int, int d=3>
Teuchos::RCP< ProcGrid<O,d> > createProcGrid(
    const Teuchos::RCP< GridSizeLocal<O,d> >& gsl,
    const Teuchos::RCP< BoundaryConditionsGlobal>& bcg,
    const Teuchos::RCP< ProcGridSize<O,d> >& procGridSize ) {

  return(
      Teuchos::rcp( new ProcGrid<O,d>( gsl, bcg, procGridSize ) ) );

}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PROCGRID_HPP
