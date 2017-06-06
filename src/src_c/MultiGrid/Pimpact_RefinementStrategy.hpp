#pragma once
#ifndef PIMPACT_REFINEMENTSTRATEGY_HPP
#define PIMPACT_REFINEMENTSTRATEGY_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {



/// \brief simple refinement
template<class SpaceT>
class RefinementStrategy {


  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  /// should be same for finest space and coarse spaces
  static const int sdim = SpaceT::sdim;
  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  using TO = typename Teuchos::Tuple<Ordinal,dimension>;

public:

  static Teuchos::RCP< const SpaceT > createRefinedSpace(
    const Teuchos::RCP<const SpaceT>& space,
    const Teuchos::Tuple<int,4>& refine_dir ) {

    auto stencilWidths = space->getStencilWidths();

    auto domainSize = space->getDomainSize();
    auto boundaryConditionsGlobal = space->getBCGlobal();
    auto boundaryConditionsLocal = space->getBCLocal();

    // refine global gridsize
    Teuchos::Tuple<Ordinal,4> gridSizeGlobalTup = *space->getGridSizeGlobal();

    for( int i=0; i<3; ++i )
      if( refine_dir[i] )
        gridSizeGlobalTup[i] = (gridSizeGlobalTup[i]-1)*2 +1;
    if( refine_dir[3] )
      gridSizeGlobalTup[3] = gridSizeGlobalTup[3]+refine_dir[3];

    auto gridSizeGlobal = createGridSizeGlobal<Ordinal,sdim>( gridSizeGlobalTup );

    auto procGrid = space->getProcGrid();

    auto gridSizeLocal =
      Pimpact::createGridSizeLocal<Ordinal,sdim,dimension,dimNC>(
        gridSizeGlobal,
        procGrid,
        stencilWidths );

    auto indexSpace =
      Pimpact::createIndexSpace<Ordinal,sdim,dimension>(
        stencilWidths,
        gridSizeLocal,
        boundaryConditionsLocal,
        procGrid );

    auto coordGlobal =
      Pimpact::createCoordinatesGlobal<Scalar,Ordinal,sdim,dimension>(
        gridSizeGlobal,
        domainSize,
        space->getCoordinatesGlobal()->getStretchParameter() );

    auto coordLocal = Pimpact::createCoordinatesLocal<Scalar,Ordinal,sdim,dimension>(
                        stencilWidths,
                        domainSize,
                        gridSizeGlobal,
                        gridSizeLocal,
                        boundaryConditionsGlobal,
                        boundaryConditionsLocal,
                        procGrid,
                        coordGlobal );

    auto interV2S =
      Pimpact::createInterpolateV2S<Scalar,Ordinal,sdim,dimension>(
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
                interV2S ) ) );
  }


}; // end of class RefinementStrategy

} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::RefinementStrategy< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::RefinementStrategy< Pimpact::Space<double,int,3,4> >;
#endif

#endif // end of #ifndef PIMPACT_REFINEMENTSTRATEGY_HPP
