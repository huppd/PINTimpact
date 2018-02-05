#pragma once
#ifndef PIMPACT_REFINEMENTSTRATEGY_HPP
#define PIMPACT_REFINEMENTSTRATEGY_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Grid.hpp"




namespace Pimpact {



/// \brief simple refinement
template<class GridT>
class RefinementStrategy {


  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  /// should be same for finest grid and coarse grids
  static const int sdim = GridT::sdim;
  static const int dimension = GridT::dimension;

  static const int dimNC = GridT::dimNC;

  using TO = typename Teuchos::Tuple<Ordinal, dimension>;

public:

  static Teuchos::RCP<const GridT > createRefinedGrid(
    const Teuchos::RCP<const GridT>& grid,
    const Teuchos::Tuple<int, 4>& refine_dir) {

    auto stencilWidths = grid->getStencilWidths();

    auto domainSize = grid->getDomainSize();
    auto boundaryConditionsGlobal = grid->getBCGlobal();
    auto boundaryConditionsLocal = grid->getBCLocal();

    // refine global gridsize
    Teuchos::Tuple<Ordinal, 4> gridSizeGlobalTup = *grid->getGridSizeGlobal();

    for(int i=0; i<3; ++i)
      if(refine_dir[i])
        gridSizeGlobalTup[i] = (gridSizeGlobalTup[i]-1)*2 + 1;
    if(refine_dir[3])
      gridSizeGlobalTup[3] = gridSizeGlobalTup[3]+refine_dir[3];

    auto gridSizeGlobal = createGridSizeGlobal<Ordinal, sdim>(gridSizeGlobalTup);

    auto procGrid = grid->getProcGrid();

    auto gridSizeLocal =
      Pimpact::createGridSizeLocal<Ordinal, sdim, dimension, dimNC>(
        gridSizeGlobal,
        procGrid,
        stencilWidths);

    auto indexSpace =
      Pimpact::createIndexSpace<Ordinal, sdim, dimension>(
        stencilWidths,
        gridSizeLocal,
        boundaryConditionsLocal,
        procGrid);

    auto coordGlobal =
      Pimpact::createCoordinatesGlobal<Scalar, Ordinal, sdim, dimension>(
        gridSizeGlobal,
        domainSize,
        grid->getCoordinatesGlobal()->getStretchParameter());

    auto coordLocal = Pimpact::createCoordinatesLocal<Scalar, Ordinal, sdim, dimension>(
                        stencilWidths,
                        domainSize,
                        gridSizeGlobal,
                        gridSizeLocal,
                        boundaryConditionsGlobal,
                        boundaryConditionsLocal,
                        procGrid,
                        coordGlobal);

    auto interV2S =
      Pimpact::createInterpolateV2S<Scalar, Ordinal, sdim, dimension>(
        indexSpace,
        gridSizeLocal,
        stencilWidths,
        domainSize,
        boundaryConditionsLocal,
        coordLocal);

    return Teuchos::rcp(
        new GridT(
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
          interV2S));
  }


}; // end of class RefinementStrategy

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_REFINEMENTSTRATEGY_HPP
