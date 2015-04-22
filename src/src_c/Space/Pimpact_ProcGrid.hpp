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




namespace Pimpact{



/// \brief ProcGrid, needs ProcGridSize, globalBoundaryConditions, GridSizeLocal
/// provides information about neighboring mpi processes
/// \ingroup Space
/// \tparam dim  3 or 4
template< class Ordinal, int dim>
class ProcGrid {

  template< class OT, int dT >
  friend Teuchos::RCP<const ProcGrid<OT,dT> > createProcGrid();

  template< class OT, int dT >
  friend Teuchos::RCP<const ProcGrid<OT,dT> > createProcGrid(
      const Teuchos::RCP<const GridSizeLocal<OT,dT> >& gsl,
      const Teuchos::RCP<const BoundaryConditionsGlobal<dT> >& bcg,
      const Teuchos::RCP<const ProcGridSize<OT,dT> >& procGridSize );

protected:

	/// used for multigrid if proc is involved in coarser grid
	bool participating_;

	/// communicator for 4 dimensional procGrid (used for exchange, reductions(for dim==4))
  MPI_Comm commWorld_;

	/// communicator for 3 dimensional sub procGrid of 4dim  (used for reduction, SF_write )
  MPI_Comm commSub_;

  int rankWorld_;
  int rankSub_;

  /// processor coordinates(index Block) fortranstyle going from 1..np
  Teuchos::Tuple<int,dim> iB_;

  /// index offset going for dim 1:3 strange fortran style dim==4 c style
  Teuchos::Tuple<Ordinal,dim> shift_;

  /// rank of lower neighbour
  Teuchos::Tuple<int,dim> rankL_;

  /// rank of upper neighbour
  Teuchos::Tuple<int,dim> rankU_;

//public:
//
//  MPI_Comm commSlice_[3] ;
//  MPI_Comm commBar_[3];
//
//  Teuchos::Tuple<int,3> rankSlice_;
//  Teuchos::Tuple<int,3> rankBar_;
//

  ProcGrid(
      const Teuchos::RCP<const GridSizeLocal<Ordinal,dim> >& gridSizeLocal,
      const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcg,
      const Teuchos::RCP<const ProcGridSize<Ordinal,dim> >& procGridSize ):
        participating_(true),
				commWorld_(),
				commSub_(),
				rankWorld_(0),
				rankSub_(0),
				iB_(),
				shift_(),
				rankL_(),
				rankU_()//,commSlice_(),commBar_(),rankSlice_(),rankBar_()
	{

    TEUCHOS_TEST_FOR_EXCEPTION(
        3!=dim && 4!=dim ,
        std::logic_error,
        "Pimpact::ProcGrid::cotor: dim:"<<dim<<" not implemented!!!\n");

    Teuchos::Tuple<int,dim> ijkB;      // mpi grid coordinates
    Teuchos::Tuple<int,dim> periodic = bcg->periodic();  // array for mpir to signal where periodic grid is, from manual should be bool

    // init periodic
//    for( int i=0; i<3; ++i ) {
//      if( bcg->getBCL(i)==PeriodicBC )
//        periodic[i] = 1;
//      else
//        periodic[i] = 0;
//    }
//    if( 4==dim ) periodic[3] = 1;


    // true means ranking may be reorderd
    // comm_cart comm with cartesian grid informations
    //    int* bla = procGridSize->getRawPtr();
    //    int mpierror =

    MPI_Cart_create(
        MPI_COMM_WORLD,
        dim,
        procGridSize->get(),
        periodic.getRawPtr(),
        true,
        &commWorld_ );


    if( 3==dim ) {

      commSub_ = commWorld_;

      // gets rank from COMM_CART
      MPI_Comm_rank( commSub_, &rankSub_ );
      rankWorld_ = rankSub_;

    }
    else if( 4==dim ) {

			int temp[] = {1,1,1,0};
			MPI_Cart_sub( commWorld_, temp, &commSub_ );

      // gets rank from COMM_CART
      MPI_Comm_rank( commWorld_, &rankWorld_ );
      MPI_Comm_rank( commSub_, &rankSub_ );
      //      std::cout << "rank in commSub: " << rankSub_ << "\trank in commWorld: " << rankWorld_ << "\n";
    }

    // gets coordinates in xyz direction from rankWorld and commWorld
    MPI_Cart_coords(
        commWorld_,
        rankWorld_,
        dim,
        ijkB.getRawPtr() );

		//    std::cout << "rankWorld: " << rankWorld_ << " coord: " << ijkB << "\n";


    // stores coordinates in a fortran fasion
    for( int i=0; i<dim; ++i )
      iB_[i] = ijkB[i] + 1;

    // computes index offset
    for( int i=0; i<3; ++i )
      shift_[i] = (iB_[i]-1)*( gridSizeLocal->get(i)-1 );
    if( 4==dim )
      shift_[3] = (iB_[3]-1)*( gridSizeLocal->get(3) );


    for( int i = 0; i<dim; ++i )
      MPI_Cart_shift( commWorld_, i, 1, &rankL_[i], &rankU_[i] );
    //                            ^  ^          ^      ^
    //                            |  |          |      |
    //                            d  d          r      r
    //                            i  i          a      a
    //                            r  s          n      n
    //                            e  p          k      k
    //                            c  l          s      d
    //                            t  a          o      e
    //                            i  c          u      s
    //                            o  m          r      t
    //                            n  e          c
    //                               n          e
    //                               t       


//    {
//      int temp[] = {0,1,1};
//      MPI_Cart_sub( commSub_, temp, &commSlice_[0] );
//    }{
//      int temp[] = {1,0,1};
//      MPI_Cart_sub( commSub_, temp, &commSlice_[1] );
//    }{
//      int temp[] = {1,1,0};
//      MPI_Cart_sub( commSub_, temp, &commSlice_[2] );
//    }

//    {
//      int temp[] = {1,0,0};
//      MPI_Cart_sub( commSub_, temp, &commBar_[0] );
//    }{
//      int temp[] = {0,1,0};
//      MPI_Cart_sub( commSub_, temp, &commBar_[1] );
//    }{
//      int temp[] = {0,0,1};
//      MPI_Cart_sub( commSub_, temp, &commBar_[2] );
//    }

//    for( int i=0; i<3; ++i ) {
//      MPI_Comm_rank( commSlice_[i], &rankSlice_[i] );
//
//      MPI_Comm_rank( commBar_[i], &rankBar_[i] );
//    }


    // maybe dispensable...
    MPI_Errhandler_set(commWorld_, MPI_ERRORS_ARE_FATAL );
    MPI_Errhandler_set(commSub_, MPI_ERRORS_ARE_FATAL );
//    for( int i=0; i<3; ++i ) {
//      MPI_Errhandler_set( commSlice_[i], MPI_ERRORS_ARE_FATAL );
//      MPI_Errhandler_set( commBar_[i],   MPI_ERRORS_ARE_FATAL );
//    }

  }

public:

	/// \note should only be used if knowing what to do(e.g. MG)
	/// \warning no safety checks
	ProcGrid(
			bool participating,
			MPI_Comm commWorld,
			MPI_Comm commSub,
			int rankWorld,
			int rankSub,
			Teuchos::Tuple<int,dim> ib,
			Teuchos::Tuple<Ordinal,dim> shift,
			Teuchos::Tuple<int,dim> rankL,
			Teuchos::Tuple<int,dim> rankU ):
		participating_(participating),
		commWorld_(commWorld),
		commSub_(commSub),
		rankWorld_(rankWorld),
		rankSub_(rankSub),
		iB_(ib),
		shift_(shift),
		rankL_(rankL),
		rankU_(rankU)	{}

  void print( std::ostream& out=std::cout ) const {
    out << "\t---ProcessorGrid: ---\n";
    out << "\tparticipating: " <<participating_<<"\n";
    out << "\trankSub: " <<rankSub_<<"\n";
    out << "\trankWorld: " <<rankWorld_<<"\n";
    out << "\trankL: " <<rankL_<<"\n";
    out << "\trankU: " <<rankU_<<"\n";
    out << "\tproc coordinate: " << iB_ << "\n";
    out << "\toffset: " << shift_ << "\n";
  }

  const bool& participating() const { return( participating_ ); }

  const MPI_Comm& getCommWorld() const { return( commWorld_ ); }
  const MPI_Comm& getCommS() const { return( commSub_ ); }

  const int& getRank() const { return( rankWorld_ ); }

  const int& getRankS() const { return( rankSub_ ); }

  const int* getRankL() const { return( rankL_.getRawPtr() ); }
  const int* getRankU() const { return( rankU_.getRawPtr() ); }

  const int& getIB( int i ) const { return( iB_[i] ); }
  const int* getIB() const { return( iB_.getRawPtr() ); }

  const int& getShift( int i ) const { return( shift_[i] ); }
  const int* getShift() const { return( shift_.getRawPtr() ); }

  const int& getRankL( int i ) const { return( rankL_[i] ); }
  const int& getRankU( int i ) const { return( rankU_[i] ); }

}; // end of class ProcGrid



/// \relates ProcGrid
template<class O, int d>
Teuchos::RCP<const ProcGrid<O,d> > createProcGrid(
    const Teuchos::RCP<const GridSizeLocal<O,d> >& gsl,
    const Teuchos::RCP<const BoundaryConditionsGlobal<d> >& bcg,
    const Teuchos::RCP<const ProcGridSize<O,d> >& procGridSize ) {

  return(
      Teuchos::rcp( new ProcGrid<O,d>( gsl, bcg, procGridSize ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_PROCGRID_HPP
