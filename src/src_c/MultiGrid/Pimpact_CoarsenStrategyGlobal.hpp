#pragma once
#ifndef PIMPACT_COARSENSTRATEGYGLOBAL_HPP
#define PIMPACT_COARSENSTRATEGYGLOBAL_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {



/// \brief first model implementation, where one coarsens until every processor has 3 dofs.
///
/// \tparam FSpaceT space on finest level not necessary the same for the coarse grids( difference in \c dim_nc).
/// it should be that the coardinates are taken from the fine grid and are halved.
/// \todo compute coordinates correctly (grid stretching only on finest grid, afterwards simple coarsening),
/// cleanest version would be to use same grid stretching on every level, makes
/// interpolation and restriction slightly more complicated.
///\todo the \c CSpaceT parameter allows even smaller Grids nLoc=1 but therefore
///the exception in GridSilzeLocal has to be adapted on StencilWidths
/// \ingroup MG
/// \todo add template parameter for coarses gridSize, and some procGrid stuff
template<class FSpaceT,class CSpaceT, int cgsize=9>
class CoarsenStrategyGlobal {

  typedef typename FSpaceT::Scalar  Scalar;
  typedef typename FSpaceT::Ordinal Ordinal;

  /// should be same for finest space and coarse spaces
  static const int dimension = FSpaceT::dimension;

	typedef typename Teuchos::Tuple<Ordinal,dimension> TO;

  /// can be different for finest and coarse spaces
  static const int dimNCF = FSpaceT::dimNC;

  static const int dimNCC = CSpaceT::dimNC;

public:

	// \todo make interface if spectral refinment is desired or not
  static std::vector<Teuchos::RCP<const CSpaceT> > getMultiSpace(
      const Teuchos::RCP<const FSpaceT> space,
      int maxGrids=10 ) {


		Teuchos::RCP<const CSpaceT> tempSpace = createSpace<CSpaceT,FSpaceT>( space );

    std::vector<Teuchos::RCP<const CSpaceT> > multiSpace( 1, tempSpace );

//		Teuchos::Tuple<Ordinal,4> nGlo = *space->getGridSizeGlobal();
		GridSizeGlobal<Ordinal,dimension> nGlo = *space->getGridSizeGlobal();

		bool spectralT = space->getStencilWidths()->spectralT();

		Teuchos::Tuple<bool,4> coarsen_dir;

    TO npWorld = space->getProcGrid()->getNP();
		TO np = npWorld;
		TO npNew = np;
    TO ibWorld = space->getProcGrid()->getIB();
		TO stride;

		for( Ordinal i=0; i<dimension; ++i ) 
			stride[i] = 1;

		// just to figure out amount of grids
		for( int i=1; i<maxGrids; ++i ) {
			bool coarsen_yes = false;
			bool procGridChanged=false;

			for( int dir=0; dir<dimension; ++dir ) {
				if( dir<3 ) {
					// figure out if global problem can be coarsened
					if( ( (nGlo[dir]-1)%2 )==0 && nGlo[dir]>=cgsize ) {
						nGlo[dir] = (nGlo[dir]-1)/2 + 1;
						coarsen_yes = true;

					}
					// figure out if processor grid has to be coarsened
					npNew[dir] = 1;
					for( int p=2; p<=np[dir]; ++p ) {
						if( ((nGlo[dir]-1)%p)==0 && (np[dir]%p)==0 )
							if( ((nGlo[dir]-1)/p)%2==0 && (nGlo[dir]-1)/p>=2 )
								npNew[dir] = p;
					}
					//--- enforce gathering of coarsest grid to one processor ---
					if( ( (nGlo[dir]-1)%2!=0 || nGlo[dir]<cgsize || i==maxGrids-1 ) )
						npNew[dir] = 1;

				}
				else
					// figure out if global problem can be coarsened
					if( (!spectralT) && ( nGlo[dir]%2 )==0 && nGlo[dir]>1 && nGlo[dir]/2>=space->getProcGrid()->getNP(dir) ) {
						nGlo[dir] = nGlo[dir]/2;
						coarsen_yes = true;
					}
				if( np[dir]!=npNew[dir] )
					procGridChanged=true;

			}

//			std::cout << "grid: " << i << "\n";
//			std::cout << "corasen_yes: " << coarsen_yes << "\n";
//			std::cout << "procGridChanged: " << procGridChanged << "\n";
//			std::cout << "nGlo: " << nGlo << "\n";
//			std::cout << "npWorld: " << npWorld << "\n";
//			std::cout << "np: " << np << "\n";
//			std::cout << "npNew: " << npNew << "\n";

			if( coarsen_yes ) {
				if( procGridChanged )
//					multiSpace.push_back( createCoarseSpace( multiSpace.back(), nGlo, coarsen_dir, stride, npWorld, ibWorld, i==maxGrids-1 ) );
					multiSpace.push_back( createCoarseSpace( multiSpace.back(), nGlo, npNew, stride, npWorld, ibWorld ) );
				else
					multiSpace.push_back( createSpace( multiSpace.back(), nGlo ) );
			}
			np = npNew;


		}

		// not working on brutus
		//multiSpace.shrink_to_fit();

    return( multiSpace );

  }

protected:

	/// \todo redo input GridSizeGlobal, ProcGridSize
	template<class SpaceT>
	static Teuchos::RCP< const SpaceT > createCoarseSpace(
			const Teuchos::RCP<const SpaceT>& space,
			const GridSizeGlobal<typename SpaceT::Ordinal,SpaceT::dimension>& newGridSizeGlobal,
			const TO& npNew,
			TO& stride,
			const TO& npWorld,
			const TO& ibWorld ) {

		auto stencilWidths = space->getStencilWidths();

		auto domainSize = space->getDomainSize();
		auto boundaryConditionsGlobal = space->getBCGlobal();


		// --- coarsen gridSizeGlobal ---
		auto gridSizeGlobal =
			createGridSizeGlobal<Ordinal,dimension>( newGridSizeGlobal );


		/// redo procGrid create new communicator, make procgrid from communicator
		TO np = space->getProcGrid()->getNP();
		TO nGather;
		TO ib;
		for( Ordinal dir=0; dir<dimension; ++dir ) {
			nGather[dir] = np[dir]/npNew[dir];
			stride[dir] *= nGather[dir];
			ib[dir]      = ( ibWorld[dir] - 1 )*npNew[dir]/npWorld[dir] + 1;
		}

		bool participating = false;
		MPI_Comm commWorld = space->getProcGrid()->getCommWorld();
		MPI_Comm commSub ;//= space->getProcGrid()->getCommS();

		int rankWorld = space->getProcGrid()->getRank();
		int rankSub = space->getProcGrid()->getRank(); // necessary?

		Teuchos::Tuple<int,dimension> rankL;
		Teuchos::Tuple<int,dimension> rankU;
		for( int dir=0; dir<dimension; ++dir )
			MPI_Cart_shift(
					commWorld,      // communicator with Cartesian structure
					dir,            // coordinate dimension of shift
					stride[dir],    // displacement
					&rankL[dir],    // rank of source process
					&rankU[dir] );	// rank of destination process

		Ordinal gather_yes = 1;
		for( int i=0; i<dimension; ++i ) {
			gather_yes *= nGather[i];
		}
		if( gather_yes>1 ) {
			int n = 1;
			for( int i=0; i<dimension; ++i ) 
				n *= npNew[i];
			int* newRanks = new int[n];

			TO rankCoord;

			for( int i=0; i<npNew[0]; ++i) {
				rankCoord[0] = (i*stride[0])%npWorld[0];
				for( int j=0; j<npNew[1]; ++j ) {
					rankCoord[1] = (j*stride[1])%npWorld[1];
					for( int k=0; k<npNew[2]; ++k ) {
						rankCoord[2] = (k*stride[2])%npWorld[2];
						if( 4==dimension ) {
							for( int l=0; l<npNew[3]; ++l ) {
								rankCoord[3] = (l*stride[3])%npWorld[3];
								MPI_Cart_rank(
										commWorld,									                    // comm
										rankCoord.getRawPtr(),                          // processor coordinate
										&newRanks[i
										         +j*npNew[0]
										         +k*npNew[0]*npNew[1]
										         +l*npNew[0]*npNew[1]*npNew[2] ] );   	// according rank to coordinate
								if( rankWorld==newRanks[i + j*npNew[0] + k*npNew[0]*npNew[1] + l*npNew[0]*npNew[1]*npNew[2] ] )
									participating = true;
							}
						}
						else {
							MPI_Cart_rank(
									commWorld,									                    // comm
									rankCoord.getRawPtr(),                          // processor coordinate
									&newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ] ); // according rank to coordinate
							if( rankWorld==newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ])
								participating = true;
						}
					}
				}
			}

			MPI_Comm commTemp;
			MPI_Group baseGroup, newGroup;

			MPI_Comm_group( commWorld, &baseGroup );
			MPI_Group_incl( baseGroup, n, newRanks, &newGroup );
			MPI_Comm_create( commWorld, newGroup, &commTemp );
			MPI_Group_free( &baseGroup );
			MPI_Group_free( &newGroup );

			Teuchos::Tuple<int,dimension> periodic = boundaryConditionsGlobal->periodic();

			if( participating ) {

				MPI_Cart_create(
						commTemp,		          // communicator without Cartesian information
						dimension,            // number of dimensions
						npNew.getRawPtr(),    // number of processors in each dimension
						periodic.getRawPtr(),	// array for mpi to signal which dimension is periodic
						false,                // false means ranking is not reordered
						&commSub );           // new communicator with Cartesian information
				if( 4==dimension ) {
					MPI_Comm commTemp_;
					int temp[] = {1,1,1,0};
					MPI_Cart_sub( commSub, temp, &commTemp_ );
					MPI_Comm_free( &commSub );
					commSub = commTemp_;
				}
				MPI_Comm_free( &commTemp );
			}
			else
				commSub=MPI_COMM_NULL;

			delete[] newRanks;

		}

		if( commSub==MPI_COMM_NULL )
			rankSub = -1;
		else 
			MPI_Comm_rank( commSub, &rankSub ); // get rank

		auto procGrid =
			Teuchos::rcp(
					new ProcGrid<Ordinal,dimension>(
						npNew,
						participating,
						commWorld,
						commSub,
						rankWorld,
						rankSub,
						ib,
						rankL,
						rankU ) );


	auto gridSizeLocal =
		Pimpact::createGridSizeLocal<Ordinal,dimension,dimNCC>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

	auto boundaryConditionsLocal =
		createBoudaryConditionsLocal<Ordinal,dimension>( 
				boundaryConditionsGlobal,
				procGrid );

	auto indexSpace =
		Pimpact::createIndexSpace<Ordinal,dimension>(
				stencilWidths,
				gridSizeLocal,
				boundaryConditionsLocal,
				procGrid );

	auto coordGlobal =
		Pimpact::createGridCoordinatesGlobal<Scalar,Ordinal,dimension>(
				gridSizeGlobal,
				domainSize,
				Teuchos::tuple( None, None, None) );

	auto  coordLocal =
		Pimpact::createGridCoordinatesLocal<Scalar,Ordinal,dimension>(
				stencilWidths,
				domainSize,
				gridSizeGlobal,
				gridSizeLocal,
				boundaryConditionsGlobal,
				boundaryConditionsLocal,
				procGrid,
				coordGlobal );

	auto interV2S =
		Pimpact::createInterpolateV2S<Scalar,Ordinal,dimension>(
				indexSpace,
				gridSizeLocal,
				stencilWidths,
				domainSize,
				boundaryConditionsLocal,
				coordLocal );

	return(
			 Teuchos::rcp(
					 new SpaceT(
							 stencilWidths,
							 indexSpace,
							 gridSizeGlobal,
							 gridSizeLocal,
							 procGrid,
							 coordGlobal,
							 coordLocal,
							 domainSize,
							 boundaryConditionsGlobal,
							 boundaryConditionsLocal,
							 interV2S )
			 )
	);

  }


}; // end of class CoarsenStrategyGlobal

} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::CoarsenStrategyGlobal< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::CoarsenStrategyGlobal< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >;
#endif

#endif // end of #ifndef PIMPACT_COARSENSTRATEGYGLOBAL_HPP
