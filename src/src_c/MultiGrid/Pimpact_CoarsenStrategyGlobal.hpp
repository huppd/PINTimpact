#pragma once
#ifndef PIMPACT_COARSENSTRATEGYGLOBAL_HPP
#define PIMPACT_COARSENSTRATEGYGLOBAL_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {



/// \brief first model implementation, where one coarsens until every processor has 3 dofs.
///
/// \tparam FSpaceT space on finest level not necessary the same for the coarse grids( difference in \c dim_nc).
/// it should be that the coardinates are taken from the fine grid and are halved.
/// \note compute coordinates correctly (grid stretching only on finest grid, afterwards simple coarsening),
/// cleanest version would be to use same grid stretching on every level, makes
/// interpolation and restriction slightly more complicated.
///\note the \c CSpaceT parameter allows even smaller Grids nLoc=1 but therefore
///the exception in GridSilzeLocal has to be adapted on StencilWidths
/// \ingroup MG
/// \note add template parameter for coarses gridSize, and some procGrid stuff
template<class FSpaceT,class CST, int cgsize=9>
class CoarsenStrategyGlobal {

public:

	using SpaceT = FSpaceT;
	using CSpaceT = CST;

protected:

  using Scalar = typename FSpaceT::Scalar;
  using Ordinal = typename FSpaceT::Ordinal;

  /// should be same for finest space and coarse spaces
  static const int sdim = FSpaceT::sdim;
  static const int dimension = FSpaceT::dimension;

	using TO = typename Teuchos::Tuple<Ordinal,dimension>;

  /// can be different for finest and coarse spaces
  static const int dimNCF = FSpaceT::dimNC;

  static const int dimNCC = CSpaceT::dimNC;

public:

	/// \todo make interface if spectral refinment is desired or not
	static std::vector<Teuchos::RCP<const CSpaceT> > getMultiSpace(
			const Teuchos::RCP<const FSpaceT> space,
			int maxGrids=10 ) {

		Teuchos::RCP<const CSpaceT> tempSpace = createSpace<CSpaceT,FSpaceT>( space );

		std::vector<Teuchos::RCP<const CSpaceT> > multiSpace( 1, tempSpace );

		GridSizeGlobal<Ordinal,sdim> nGlo = *space->getGridSizeGlobal();

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

					//multiSpace.push_back( createCoarseSpace( multiSpace.back(), coarsen_dir, nGlo, stride, NB, IB, i==maxGrids-1 ) );
					//}


			//}
			//std::ofstream file;
			//std::string fname = "mgs.txt";
			//fname.insert( 3, std::to_string( (long long)space->rankST() ) );
			//file.open( fname, std::ofstream::out | std::ofstream::app );
			//file << "\n\ngrid: " <<  -1 << "\n\n";
			//multiSpace.back()->print(file);
			//file.close();

		//}
		//// not working on brutus
		////multiSpace.shrink_to_fit();
		//return( multiSpace );
	//}

//protected:

  //template<class SpaceT>
  //static Teuchos::RCP< const SpaceT > createCoarseSpace(
      //const Teuchos::RCP<const SpaceT>& space,
      //const Teuchos::Tuple<bool,dimension>& coarsen_dir,
      //const Teuchos::Tuple<Ordinal,4>& gridSizeGlobalTup,
			 //TO& stride,
			//const TO& NB,
			//const TO& IB,
			//bool maxGrid_yes=false	) {

    //auto stencilWidths = space->getStencilWidths();

    //auto domain = space->getDomain();
    //auto boundaryConditionsGlobal = domain->getBCGlobal();
    //auto boundaryConditionsLocal = domain->getBCLocal();


		//// coarsen gridSizeGlobal
////		std::cout << "rank: " << space->rankST()<< "\t gridSizGLobal: " << gridSizeGlobalTup << "\n";
    //auto gridSizeGlobal = createGridSizeGlobal<Ordinal,dimension>( gridSizeGlobalTup );

		//auto procGridSize = space->getProcGridSize();
		//auto procGrid = space->getProcGrid();

    //TO np = space->getProcGridSize()->getTuple();
		//TO npNew = np;
////		std::cout << "rank: " << space->rankST()<< "\t procGridSizOld: " << npNew << "\n";
////    for( int i=0; i<dimension; ++i ) {
    //for( int i=0; i<3; ++i ) {
			//if( i<3 ) {
				//if( coarsen_dir[i] ) {
					//npNew[i] = 1;
					//for( int j=2; j<np[i]; ++j ) {
						//if( ((gridSizeGlobalTup[i]-1)%j)==0 && (np[i]%j)==0 )
							//if( ((gridSizeGlobalTup[i]-1)/j)%2==0 && (gridSizeGlobalTup[i]-1)/j>=2 )
								//npNew[i] = j;
					//}
					//// --- enforce gathering of coarsest grid to one processor ---
					//if( ( (gridSizeGlobalTup[i]-1)%2!=0 || gridSizeGlobalTup[i]<cgsize || maxGrid_yes ) && npNew[i]>1 )
						//npNew[i] = 1;
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
					multiSpace.push_back( createSpace( multiSpace.back(), nGlo, npNew, stride, npWorld, ibWorld ) );
				else
					multiSpace.push_back( createSpace( multiSpace.back(), nGlo ) );
			}
			np = npNew;


		}

		// not working on brutus
		//multiSpace.shrink_to_fit();

    return( multiSpace );

  }


}; // end of class CoarsenStrategyGlobal


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::CoarsenStrategyGlobal< Pimpact::Space<double,int,3,3,4>, Pimpact::Space<double,int,3,3,2> >;
extern template class Pimpact::CoarsenStrategyGlobal< Pimpact::Space<double,int,3,4,4>, Pimpact::Space<double,int,3,4,2> >;
#endif


#endif // end of #ifndef PIMPACT_COARSENSTRATEGYGLOBAL_HPP
