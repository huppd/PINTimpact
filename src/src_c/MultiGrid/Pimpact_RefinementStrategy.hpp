#pragma once
#ifndef PIMPACT_REFINEMENTSTRATEGY_HPP
#define PIMPACT_REFINEMENTSTRATEGY_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {



/// \brief simple refinement
template<class SpaceT>
class RefinementStrategy {


  
	typedef typename SpaceT::Scalar  Scalar;
	typedef typename SpaceT::Ordinal Ordinal;

	/// should be same for finest space and coarse spaces
	static const int dimension = SpaceT::dimension;

	static const int dimNC = SpaceT::dimNC;

	typedef typename Teuchos::Tuple<Ordinal,dimension> TO;


public:

	static Teuchos::RCP< const SpaceT > createRefinedSpace(
			const Teuchos::RCP<const SpaceT>& space,
			const Teuchos::Tuple<bool,4>& refine_dir ) {

		auto stencilWidths = space->getStencilWidths();

    auto domain = space->getDomain();
    auto boundaryConditionsGlobal = domain->getBCGlobal();
    auto boundaryConditionsLocal = domain->getBCLocal();

		// refine global gridsize
    auto gridSizeGlobalTup = space->getGridSizeGlobal()->getTuple();

		for( int i=0; i<3; ++i )
		 if( refine_dir[i] )
			 gridSizeGlobalTup[i] = (gridSizeGlobalTup[i]-1)*2 +1;
		if( refine_dir[3] )
			gridSizeGlobalTup[3] = gridSizeGlobalTup[3]+2;

    auto gridSizeGlobal = createGridSizeGlobal<Ordinal,dimension>( gridSizeGlobalTup );

    auto procGridSize = space->getProcGridSize();

    auto gridSizeLocal =
			Pimpact::createGridSizeLocal<Ordinal,dimension,dimNC>(
					gridSizeGlobal, procGridSize, stencilWidths );

    auto procGrid = Pimpact::createProcGrid<Ordinal,dimension>(
        gridSizeLocal, domain->getBCGlobal(), procGridSize );

    auto indexSpace =
			Pimpact::createIndexSpace<Ordinal,dimension>(
					stencilWidths, gridSizeLocal, boundaryConditionsLocal );

    auto coordGlobal = Pimpact::createGridCoordinatesGlobal<Scalar,Ordinal,dimension>(
        gridSizeGlobal,
        domain->getDomainSize() );

    auto coordLocal = Pimpact::createGridCoordinatesLocal<Scalar,Ordinal,dimension>(
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


}; // end of class RefinementStrategy

} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::RefinementStrategy< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::RefinementStrategy< Pimpact::Space<double,int,3,4> >;
#endif

#endif // end of #ifndef PIMPACT_REFINEMENTSTRATEGY_HPP