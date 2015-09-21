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

  static std::vector<Teuchos::RCP<const CSpaceT> > getMultiSpace(
      const Teuchos::RCP<const FSpaceT> space,
      int maxGrids=10 ) {


		auto tempSpace = createCoarseSpaceT( space );

    std::vector<Teuchos::RCP<const CSpaceT> > multiSpace( 1, tempSpace );

		Teuchos::Tuple<Ordinal,4> nGlo = space->getGridSizeGlobal()->getTuple();

//    for( int i=0; i<dimension; ++i )
//      nGlo[i] = space->nGlo(i);

    bool coarsen_yes;
		Teuchos::Tuple<bool,4> coarsen_dir;
    TO NB;
    TO IB;

		TO stride;
		for( Ordinal i=0; i<dimension; ++i ) {
			stride[i] = 1;
			NB[i] = space->getNProc( i );
			IB[i] = space->getProcGrid()->getIB( i );
		}

		// just to figure out amount of grids
    for( int i=1; i<maxGrids; ++i ) {
      coarsen_yes = false;

//			std::cout << "rank: " << space->rankST() << "nGlo: " << nGlo << "\n";
      for( int j=0; j<dimension; ++j ) {
        coarsen_dir[j] = false;
        if( j<3 ) {
				 if( ( (nGlo[j]-1)%2 )==0 && nGlo[j]>=cgsize ) {
            nGlo[j] = (nGlo[j]-1)/2 + 1;
            coarsen_yes = true;
            coarsen_dir[j] = true;
          }
        }
        else
          if( ( (nGlo[j])%2 )==0 && nGlo[j]>1 && (nGlo[j]/2)>=space->getNProc(j) ) {
            nGlo[j] = (nGlo[j])/2;
            coarsen_yes = true;
            coarsen_dir[j] = true;
          }
      }
			if( 3==dimension ) {
				if( nGlo[3]>1 ) {
					nGlo[3] = (nGlo[3]) - 1;
					coarsen_yes = true;
					coarsen_dir[3] = true;
				}
			}

			if( coarsen_yes ) {
//				std::ofstream file;
//				std::string fname = "mgs.txt";
//				fname.insert( 3, std::to_string( (long long)space->rankST() ) );
//				file.open( fname, std::ofstream::out | std::ofstream::app );
//				file << "\n\ngrid: " <<  i << "\n\n";
//				multiSpace.back()->print(file);
//				file.close();

				multiSpace.push_back( createCoarseSpace( multiSpace.back(), coarsen_dir, nGlo, stride, NB, IB, i==maxGrids-1 ) );
			}


    }
//		std::ofstream file;
//		std::string fname = "mgs.txt";
//		fname.insert( 3, std::to_string( (long long)space->rankST() ) );
//		file.open( fname, std::ofstream::out | std::ofstream::app );
//		file << "\n\ngrid: " <<  -1 << "\n\n";
//		multiSpace.back()->print(file);
//		file.close();


		// not working on brutus
		//multiSpace.shrink_to_fit();

    return( multiSpace );

  }

protected:

  template<class SpaceT>
  static Teuchos::RCP< const SpaceT > createCoarseSpace(
      const Teuchos::RCP<const SpaceT>& space,
      const Teuchos::Tuple<bool,4>& coarsen_dir,
      const Teuchos::Tuple<Ordinal,4>& gridSizeGlobalTup,
		 	TO& stride,
			const TO& NB,
			const TO& IB,
			bool maxGrid_yes=false	) {

    auto stencilWidths = space->getStencilWidths();

    auto domain = space->getDomain();
    auto boundaryConditionsGlobal = domain->getBCGlobal();
    auto boundaryConditionsLocal = domain->getBCLocal();


		// coarsen gridSizeGlobal
//		std::cout << "rank: " << space->rankST()<< "\t gridSizGLobal: " << gridSizeGlobalTup << "\n";
    auto gridSizeGlobal = createGridSizeGlobal<Ordinal,dimension>( gridSizeGlobalTup );

		auto procGridSize = space->getProcGridSize();
		auto procGrid = space->getProcGrid();

    TO np = space->getProcGridSize()->getTuple();
		TO npNew = np;
//		std::cout << "rank: " << space->rankST()<< "\t procGridSizOld: " << npNew << "\n";
    for( int i=0; i<dimension; ++i ) {
			if( i<3 ) {
				if( coarsen_dir[i] ) {
					npNew[i] = 1;
					for( int j=2; j<np[i]; ++j ) {
						if( ((gridSizeGlobalTup[i]-1)%j)==0 && (np[i]%j)==0 )
							if( ((gridSizeGlobalTup[i]-1)/j)%2==0 && (gridSizeGlobalTup[i]-1)/j>=2 )
								npNew[i] = j;
					}
					// --- enforce gathering of coarsest grid to one processor ---
					if( ( (gridSizeGlobalTup[i]-1)%2!=0 || gridSizeGlobalTup[i]<cgsize || maxGrid_yes ) && npNew[i]>1 )
						npNew[i] = 1;
				}
			}
			else {
				if( coarsen_dir[3] ) {
					for( int j=1; j<np[3]; ++j ) {
						if( (gridSizeGlobalTup[3]%j)==0 && (gridSizeGlobalTup[3]/j)>=2 )
							npNew[3] = j;
					}
					if( (gridSizeGlobalTup[3]%2)!=0  && npNew[3]>1 )
						npNew[3] = 1;
				}
			}

		}
//		std::cout << "rank: " << space->rankST()<< "\t procGridSiz: " << npNew << "\n";

		bool procGridChanged=false;
		for( int i=0; i<dimension; ++i )
			if( np[i]!=npNew[i] )
				procGridChanged=true;

		auto gridSizeLocal = space->getGridSizeLocal();
//			Pimpact::createGridSizeLocal<Ordinal,dimension,dimNCC>(
//					gridSizeGlobal, procGridSize, stencilWidths );

		if( procGridChanged ) { 

			/// redo procGrid create new communicator, make procgrid from communicator
			TO nGather;
			TO ib;
			for( Ordinal i=0; i<dimension; ++i ) {
				nGather[i] = np[i]/npNew[i];
				stride[i] *= nGather[i];
				ib[i] = ( IB[i] - 1 )*npNew[i]/NB[i] + 1;
			}

			bool participating = false;
			MPI_Comm commWorld = space->getProcGrid()->getCommWorld();
			MPI_Comm commSub ;//= space->getProcGrid()->getCommS();

			int rankWorld = space->getProcGrid()->getRank();
			int rankSub = space->getProcGrid()->getRank();

			Teuchos::Tuple<int,dimension> rankL;
			Teuchos::Tuple<int,dimension> rankU;
			for( int i=0; i<dimension; ++i )
				MPI_Cart_shift( commWorld, i, stride[i], &rankL[i], &rankU[i] );

			Ordinal gather_yes = 1;
			for( int i=0; i<dimension; ++i ) {
				gather_yes *= nGather[i];
			}
			if( gather_yes>1 ) {
				int n = 1;
				for( int i=0; i<dimension; ++i ) 
					n *= npNew[i];
				int* newRanks = new int[n];
//				std::cout << "rankWorld: " << rankWorld << "n: " << n << "\n";

				TO rankCoord;

				for( int i=0; i<npNew[0]; ++i) 
					for( int j=0; j<npNew[1]; ++j )
						for( int k=0; k<npNew[2]; ++k ) {
							if( 4==dimension ) {
//								for( int l=0; l<npNew[3]; ++l ) {
//									rankCoord[0] = i; rankCoord[1]=j; rankCoord[2] = k; rankCoord[3] = l;
//									MPI_Cart_rank(
//											commWorld,
//											rankCoord.getRawPtr(),
//											&newRanks[i+j*npNew[0]+k*npNew[0]*npNew[2]+l*npNew[0]*npNew[2]*npNew[3] ] );
//									if( rankWorld==newRanks[i+j*npNew[0]+k*npNew[0]*npNew[2]+l*npNew[0]*npNew[2]*npNew[3] ] )
//										participating = true;
//								}
							}
							else {
								rankCoord[0] = (i*stride[0])%NB[0];
								rankCoord[1] = (j*stride[1])%NB[1];
								rankCoord[2] = (k*stride[2])%NB[2];
								MPI_Cart_rank(
										commWorld,
										rankCoord.getRawPtr(),
										&newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ] );
								if( rankWorld==newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ])
									participating = true;
//							std::cout << "rankWorld: " << rankWorld << "\t strid: "<< stride << "\n\tt" << "rankCOrd: " << rankCoord << "\n";

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
					MPI_Cart_create( commTemp, dimension, npNew.getRawPtr(), periodic.getRawPtr(), false, &commSub );
					MPI_Comm_free( &commTemp );

				}
				delete[] newRanks;

			}

			procGridSize = createProcGridSize( npNew, commSub, participating );
			
			gridSizeLocal =
				Pimpact::createGridSizeLocal<Ordinal,dimension,dimNCC>(
						gridSizeGlobal, procGridSize, stencilWidths );

			TO shift;
			// computes index offset
			for( int i=0; i<3; ++i )
				shift[i] = (ib[i]-1)*( gridSizeLocal->get(i)-1 );
			if( 4==dimension )
				shift[3] = (ib[3]-1)*( gridSizeLocal->get(3) );

			procGrid =
				Teuchos::rcp(
						new ProcGrid<Ordinal,dimension>(
							participating,
							commWorld,
							commSub,
							rankWorld,
							rankSub,
							ib,
							shift,
							rankL,
							rankU ) );


		} 
		else {

			gridSizeLocal =
				Pimpact::createGridSizeLocal<Ordinal,dimension,dimNCC>(
						gridSizeGlobal, procGridSize, stencilWidths );

			// necessary because offset has to be recomputed
			procGrid =
				Pimpact::createProcGrid<Ordinal,dimension>(
						gridSizeLocal,
						domain->getBCGlobal(),
						procGridSize,
						space()->getProcGrid()->participating()
						);

		}


		boundaryConditionsLocal =
			createBoudaryConditionsLocal( 
					boundaryConditionsGlobal,
					procGridSize,
					procGrid );

		domain =
			createDomain(
					space->getDomain()->getDomainSize(),
					boundaryConditionsGlobal,
					boundaryConditionsLocal );

    auto indexSpace = Pimpact::createIndexSpace<Ordinal,dimension>(
        stencilWidths, gridSizeLocal, boundaryConditionsLocal );

    auto  coordGlobal = Pimpact::createGridCoordinatesGlobal<Scalar,Ordinal,dimension>(
        gridSizeGlobal,
        domain->getDomainSize() );

    auto  coordLocal = Pimpact::createGridCoordinatesLocal<Scalar,Ordinal,dimension>(
        stencilWidths,
        domain->getDomainSize(),
        gridSizeGlobal,
        gridSizeLocal,
        domain->getBCGlobal(),
        domain->getBCLocal(),
        procGrid,
        coordGlobal );

    auto interV2S =
        Pimpact::createInterpolateV2S<Scalar,Ordinal,dimension>(
            procGrid,
            gridSizeLocal,
            stencilWidths,
            domain,
            coordLocal );

    return(
         Teuchos::rcp(
             new SpaceT(
                 stencilWidths,
                 indexSpace,
                 gridSizeGlobal,
                 gridSizeLocal,
                 procGridSize,
                 procGrid,
                 coordGlobal,
                 coordLocal,
                 domain,
                 interV2S )
         )
    );

  }

  static Teuchos::RCP< const CSpaceT > createCoarseSpaceT(
      const Teuchos::RCP<const FSpaceT>& space ) {

    auto stencilWidths = createStencilWidths<dimension,dimNCC>();

    auto domain = space->getDomain();

    auto boundaryConditionsLocal = domain->getBCLocal();

    auto procGridSize = space->getProcGridSize();

    auto gridSizeGlobalTup = space->getGridSizeGlobal()->getTuple();

    auto gridSizeGlobal = space->getGridSizeGlobal();

    auto gridSizeLocal = space->getGridSizeLocal();

    auto procGrid = space->getProcGrid();

    auto indexSpace = space->getIndexSpace();

    auto  coordGlobal = space->getCoordinatesGlobal();

    auto  coordLocal =
			Pimpact::createGridCoordinatesLocal(
					stencilWidths,
					domain->getDomainSize(),
					gridSizeGlobal,
					gridSizeLocal,
					domain->getBCGlobal(),
					domain->getBCLocal(),
					procGrid,
					coordGlobal );

		auto interV2S =
			Pimpact::createInterpolateV2S(
					procGrid,
					gridSizeLocal,
					stencilWidths,
					domain,
					coordLocal );

    return(
         Teuchos::rcp(
             new CSpaceT(
                 stencilWidths,
                 indexSpace,
                 gridSizeGlobal,
                 gridSizeLocal,
                 procGridSize,
                 procGrid,
                 coordGlobal,
                 coordLocal,
                 domain,
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
