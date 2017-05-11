#pragma once
#ifndef PIMPACT_PROCGRID_HPP
#define PIMPACT_PROCGRID_HPP


#include <ostream>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact{



/// \brief ProcGrid, needs ProcGridSize, globalBoundaryConditions
/// provides information about neighboring mpi processes
/// \ingroup SpaceObject
/// \tparam dim  3 or 4 as soon as time is own class ->sd
template< class OrdinalT, int dim>
class ProcGrid {

public:

	using TO = const Teuchos::Tuple<OrdinalT,dim>;

protected:

	template< class OT, int dT >
  friend Teuchos::RCP<const ProcGrid<OT,dT> > createProcGrid(
      const Teuchos::Tuple<OT,dT>& procGridSize,
      const Teuchos::RCP<const BoundaryConditionsGlobal<dT> >& bcg,
		 	bool participating	);

	/// processor grid size
  TO procGridSize_;

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

  /// rank of lower neighbour
  Teuchos::Tuple<int,dim> rankL_;

  /// rank of upper neighbour
  Teuchos::Tuple<int,dim> rankU_;


	MPI_Comm commBar_[dim];           ///< sub comm along dimension
	Teuchos::Tuple<int,dim> rankBar_; ///< sub rank should be equivalent to iB_ rankcoordinate??

	MPI_Comm commSlice_[3];           ///< sub comm containt slice( slice(time) not necessary, because this is in subcommunicator, a bit messy. I know
	Teuchos::Tuple<int,3> rankSlice_;

	ProcGrid(
			TO& procGridSize,
			const Teuchos::RCP<const BoundaryConditionsGlobal<dim> >& bcg,
			bool participating ):
		procGridSize_(procGridSize),
		participating_(participating),
		commWorld_(),
		commSub_(),
		rankWorld_(0),
		rankSub_(0),
		iB_(),
		rankL_(),
		rankU_(),
		commBar_(),
		rankBar_(),
		commSlice_(),
		rankSlice_() {

		//
		// --- tests ---
		// 
		static_assert( 3==dim || 4==dim, "spatial dimension not valid" );

		for( int i=0; i<dim; ++i )
			TEUCHOS_TEST_FOR_EXCEPT( procGridSize_[i]<1 );

		int commSize;
		MPI_Comm_size( MPI_COMM_WORLD, &commSize );

		int procSize = 1;
		for( int i=0; i<dim; ++i )
			procSize *= procGridSize_[i];

		TEUCHOS_TEST_FOR_EXCEPTION( procSize!=commSize, std::logic_error, "procSize: " << procSize << " != commSize: " << commSize );

		Teuchos::Tuple<int,dim> ijkB;                        // mpi grid coordinates

		//
		// -- commWorld_ ---
		// 
		MPI_Cart_create(
				MPI_COMM_WORLD,              // communicator without Cartesian information
				dim,                         // number of dimensions
				procGridSize_.getRawPtr(),   // number of processors in each dimension
				bcg->periodic().getRawPtr(), // array for mpi to signal where periodic grid is, from manual should be bool
				true,                        // true means ranking may be reorderd
				&commWorld_ );               // new communicator with Cartesian information


		//
		// -- commSub_(just spatial comm without time) ---
		// -- rankWorld_(just spatial comm without time) ---
		// -- rankSub_(just spatial comm without time) ---
		// 
		if( 3==dim ) {

			commSub_ = commWorld_;
			// gets rank from COMM_CART
			if( commSub_==MPI_COMM_NULL )
				rankSub_ = -1;
			else 
				MPI_Comm_rank( commSub_, &rankSub_ ); // get rank
			rankWorld_ = rankSub_;

		}
		else if( 4==dim ) {

			int temp[] = {1,1,1,0};
			MPI_Cart_sub( commWorld_, temp, &commSub_ );

			// gets rank from COMM_CART
			if( commWorld_==MPI_COMM_NULL )
				rankWorld_ = -1;
			else
				MPI_Comm_rank( commWorld_, &rankWorld_ );
			if( commSub_==MPI_COMM_NULL )
				rankSub_ = -1;
			else 
				MPI_Comm_rank( commSub_, &rankSub_ );
		}

		//
		// -- iB_ ---
    // gets coordinates in xyz direction from rankWorld and commWorld
		//
		if( commWorld_!=MPI_COMM_NULL )
			MPI_Cart_coords(
					commWorld_,
					rankWorld_,
					dim,
					ijkB.getRawPtr() );
		else
			for( int i=0; i<dim; ++i )
				ijkB[i] = 0;


    // stores coordinates in a fortran fasion
    for( int i=0; i<dim; ++i )
      iB_[i] = ijkB[i] + 1;


		//
		// --- rankL_ rankU_ ---
		//
		if( commWorld_!=MPI_COMM_NULL )
			for( int i = 0; i<dim; ++i )
				MPI_Cart_shift(
						commWorld_,
						i,             // direction
						1,             // displacement
						&rankL_[i],    // rank source
						&rankU_[i] );  // rank destination

		//
		// --- commBar, rankBar ---
		//
		Teuchos::Tuple<int,dim> temp;
		for( int i=0; i<dim; ++i ) {
			for( int j=0; j<dim; ++j )
				temp[j] = 0;
			temp[i] = 1;

			if( commWorld_==MPI_COMM_NULL )
				commBar_[i] = MPI_COMM_NULL;
			else
				MPI_Cart_sub(
						commWorld_,       // input communicator
						temp.getRawPtr(), // included sub dimensions
						&commBar_[i] ); // output sub communicator

			// gets rank from COMM_CART
			if( commBar_[i]==MPI_COMM_NULL )
				rankBar_[i] = -1;
			else
				MPI_Comm_rank( commBar_[i], &rankBar_[i] );
		}

		// --- commSlice, rankSlice ---
		//
		temp[dim-1] = 0; // setting time or Z zero
		for( int i=0; i<3; ++i ) {
			for( int j=0; j<3; ++j )
				temp[j] = 1;
			temp[i] = 0;

			if( commWorld_==MPI_COMM_NULL )
				commSlice_[i] = MPI_COMM_NULL;
			else
				MPI_Cart_sub(
						commWorld_,       // input communicator
						temp.getRawPtr(), // included sub dimensions
						&commSlice_[i] ); // output sub communicator

			// gets rank from COMM_CART
			if( commSlice_[i]==MPI_COMM_NULL )
				rankSlice_[i] = -1;
			else
				MPI_Comm_rank( commSlice_[i], &rankSlice_[i] );
		}

		// maybe dispensable...
		//if( commWorld_!=MPI_COMM_NULL )
			//MPI_Comm_create_errhandler(commWorld_, MPI_ERRORS_ARE_FATAL ); // is deprecated
		//if( commSub_!=MPI_COMM_NULL )
			//MPI_Comm_create_errhandler(commSub_, MPI_ERRORS_ARE_FATAL );
  }

public:

	/// \note should only be used if knowing what to do(e.g. MG)
	/// \warning no safety checks
	/// \note todo make privat, make coarsening here.
	ProcGrid(
			const TO& procGridSize,
			const bool& participating,
			const MPI_Comm& commWorld,
			const MPI_Comm& commSub,
			const int& rankWorld,
			const int& rankSub,
			const Teuchos::Tuple<int,dim>& ib,
			const Teuchos::Tuple<int,dim>& rankL,
			const Teuchos::Tuple<int,dim>& rankU ):
		procGridSize_(procGridSize),
		participating_(participating),
		commWorld_(commWorld),
		commSub_(commSub),
		rankWorld_(rankWorld),
		rankSub_(rankSub),
		iB_(ib),
		rankL_(rankL),
		rankU_(rankU)	{
			std::cout << "rank: " << rankWorld_ << "\tcomm: " << commWorld_ << "\n";
		}

	/// \todo figure out why freeing is not working, desirable for many Utest, destructor called twice?
	~ProcGrid() {
		if( commWorld_!=MPI_COMM_NULL )
			//std::cout << "rank: " << rankWorld_ << "\tcomm: " << commWorld_ << "\n";
			//MPI_Comm_free( &commWorld_ );
		if( commSub_!=MPI_COMM_NULL )
			MPI_Comm_free( &commSub_ );
		//for( int i=0; i<dim; ++i )
			//if( commSlice_[i]!=MPI_COMM_NULL )
				//MPI_Comm_free( &(commSlice_[i]) );
	}

  void print( std::ostream& out=std::cout ) const {
    out << "\t---ProcessorGrid ---\n";
		out << "\tProcGridSize: " << procGridSize_ << " ---\n";
    out << "\tparticipating: " <<participating_<<"\n";
    out << "\trankSub: " <<rankSub_<<"\n";
    out << "\trankWorld: " <<rankWorld_<<"\n";
    out << "\trankL: " <<rankL_<<"\n";
    out << "\trankU: " <<rankU_<<"\n";
    out << "\tproc coordinate: " << iB_ << "\n";
		out << "\trankBar: " << rankBar_ << "\n";
		out << "\trankSlice: " << rankSlice_ << "\n";
  }

  constexpr const bool& participating() { return( participating_ ); }

	constexpr const OrdinalT& getNP( int i ) { return( procGridSize_[i] ); }
	constexpr const TO&       getNP()        { return( procGridSize_ ); }

  constexpr const MPI_Comm& getCommWorld() { return( commWorld_ ); }
  constexpr const MPI_Comm& getCommS()     { return( commSub_ ); }

  constexpr const int& getRank()  { return( rankWorld_ ); }

  constexpr const int& getRankS() { return( rankSub_ ); }

  constexpr const int* getRankL() { return( rankL_.getRawPtr() ); }
  constexpr const int* getRankU() { return( rankU_.getRawPtr() ); }

  constexpr const int& getIB( const int& i ) { return( iB_[i] ); }
  constexpr const TO & getIB() { return( iB_ ); }

  constexpr const int& getRankL( const int& i ) { return( rankL_[i] ); }
  constexpr const int& getRankU( const int& i ) { return( rankU_[i] ); }

  constexpr const MPI_Comm& getCommBar( const int& i ) { return( commBar_[i] ); }
  constexpr const int&      getRankBar( const int& i ) { return( rankBar_[i] ); }

  constexpr const MPI_Comm& getCommSlice( const int& i ) { return( commSlice_[i] ); }
  constexpr const int&      getRankSlice( const int& i ) { return( rankSlice_[i] ); }

}; // end of class ProcGrid



/// \relates ProcGrid
template<class O, int d>
Teuchos::RCP<const ProcGrid<O,d> > createProcGrid(
    const Teuchos::Tuple<O,d>& procGridSize,
    const Teuchos::RCP<const BoundaryConditionsGlobal<d> >& bcg,
	 	bool participating=true	) {

  return(
      Teuchos::rcp(
				new ProcGrid<O,d>(
					procGridSize,
					bcg,
					participating ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_PROCGRID_HPP
