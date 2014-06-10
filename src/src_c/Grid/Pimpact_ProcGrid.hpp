#pragma once
#ifndef PIMPACT_PROCGRID_HPP
#define PIMPACT_PROCGRID_HPP

#include <ostream>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_FieldSpace.hpp"


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

//  Teuchos::Tuple<MPI_Comm,3> commSlice_;
//  Teuchos::Tuple<MPI_Comm,3> commBar_;
  MPI_Comm commSlice_[3] ;
  MPI_Comm commBar_[3];

  Teuchos::Tuple<int,3> rankSlice_;
  Teuchos::Tuple<int,3> rankBar_;

  //friend createProcGrid<Ordinal,dim>;

public:

  ProcGrid(
      const Teuchos::RCP<const FieldSpace<Ordinal> >& fieldSpace,
      const Teuchos::RCP< BoundaryConditionsGlobal >& bcg,
      const Teuchos::RCP< ProcGridSize<Ordinal,dim> >& procGridSize ):
        iB_(),rankL_(),rankU_(),commSlice_(),commBar_(),rankSlice_(),rankBar_() {

    if( 3!=dim && 4!=dim ) {
      std::cout<< "Pimpact::ProcGrid::cotor: dim:"<<dim<<" not implemented!!!\n";
    }
    Teuchos::Tuple<int,dim> ijkB;      // mpi grid coordinates
    Teuchos::Tuple<int,dim> periodic;  // array for mpir to signal where periodic grid is, from manual should be bool

    // init periodic
    for( int i=0; i<dim; ++i ) {
      if( bcg->getBCL(i)==PeriodicBC )
        periodic[i] = 1;
      else
        periodic[i] = 0;
    }

    //    ! Macht keinen messbaren Unterschied (falls doch irgendwann, dann sollte der CALL auch auf die Grobgitter-Kommunikatoren auch angewandt werden!):
    //    !CALL MPI_CART_CREATE(MPI_COMM_WORLD,3,(/NB1,NB2,NB3/),periodic,.FALSE.,COMM_CART,merror)
    //    !CALL MPI_CART_MAP(COMM_CART,3,(/NB1,NB2,NB3/),periodic,rank,merror)

    if( 3==dim ) {

      commSpaceTime_ = MPI_COMM_NULL;
      commSpaceTimef_ = MPI_Comm_c2f( commSpaceTime_ );

      // true means ranking may be reorderd
      // comm_cart comm with cartesian grid informations
      MPI_Cart_create(
          MPI_COMM_WORLD,
          dim,
          procGridSize->getRawPtr(),
          periodic.getRawPtr(),
          true,
          &commSpace_ );

      commSpacef_ = MPI_Comm_c2f( commSpace_ );

      // gets rank from COMM_CART
      MPI_Comm_rank( commSpace_, &rank_ ); //  (COMM_CART,rank,merror);
    }

    // gets coordinates in xyz direction from rank and comm_cart
    MPI_Cart_coords(
        commSpace_,
        rank_,
        dim,
        ijkB.getRawPtr() );

    // stores coordinates in a fortran fasion?
    for( int i=0; i<dim; ++i ) {
      iB_[i] =ijkB[i]+1;

      // computes index ofset
      shift_[i] = (iB_[i]-1)*( fieldSpace->nLoc_[i]-1 );

      MPI_Cart_shift( commSpace_, i, 1, &rankL_[i], &rankU_[i] );
      //      MPI_CART_SHIFT( COMM_CART, 1, 1, rank2L, rank2U, merror );
      //      MPI_CART_SHIFT( COMM_CART, 2, 1, rank3L, rank3U, merror );
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

  void print() const {
    std::cout<<"\t---ProcessorGrid: ---\n";
    std::cout <<"\trank: " <<rank_<<"\n";
    std::cout <<"\trankL: " <<rankL_<<"\n";
    std::cout <<"\trankU: " <<rankU_<<"\n";
    std::cout <<"\tiB: " << iB_ << "\n";
  }


};


/// \relates ProcGrid
template< class O, int d>
Teuchos::RCP< ProcGrid<O,d> > createProcGrid(
    const Teuchos::RCP<const FieldSpace<O> >& fieldSpace,
    const Teuchos::RCP< BoundaryConditionsGlobal>& bcg,
    const Teuchos::RCP< ProcGridSize<O,d> >& procGridSize ) {

  return(
      Teuchos::rcp( new ProcGrid< O, d >(fieldSpace,bcg,procGridSize) ) );

}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PROCGRID_HPP
