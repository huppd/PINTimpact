#pragma once
#ifndef PIMPACT_COARSENSTRATEGY_HPP
#define PIMPACT_COARSENSTRATEGY_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Grid.hpp"




namespace Pimpact {



/// \brief first model implementation, where one coarsens until every processor has 3 dofs.
///
/// \tparam FGridT grid on finest level not necessary the same for the coarse grids(difference in \c dim_nc).
/// it should be that the coardinates are taken from the fine grid and are halved.
/// \note compute coordinates correctly (grid stretching only on finest grid, afterwards simple coarsening),
/// cleanest version would be to use same grid stretching on every level, makes
/// interpolation and restriction slightly more complicated.
///\note the \c CGridT parameter allows even smaller Grids nLoc=1 but therefore
///the exception in GridSilzeLocal has to be adapted on StencilWidths
/// \ingroup MG
template<class FGridT, class CGT>
class CoarsenStrategy {

public:

  using GridT = FGridT;
  using CGridT = CGT;

protected:

  using Scalar = typename FGridT::Scalar;
  using Ordinal = typename FGridT::Ordinal;

  /// should be same for finest grid and coarse grids
  static const int sdim = FGridT::sdim;
  static const int dimension = FGridT::dimension;

  using TO = typename Teuchos::Tuple<Ordinal, dimension>;

  /// can be different for finest and coarse grids
  static const int dimNCF = FGridT::dimNC;

  static const int dimNCC = CGridT::dimNC;

public:

  static std::vector<Teuchos::RCP<const CGridT> > getMultiGrid(
    const Teuchos::RCP<const FGridT> grid,
    int maxGrids=10) {

    // creating low order grid
    Teuchos::RCP<const CGridT> tempGrid = createGrid<CGridT, FGridT>(grid);

    std::vector<Teuchos::RCP<const CGridT> > multiGrid(1, tempGrid);

    Teuchos::Tuple<Ordinal, dimension> nLoc = *grid->getGridSizeLocal();
    GridSizeGlobal<Ordinal, sdim> nGlo = *grid->getGridSizeGlobal();

    for(Ordinal i=1; i<maxGrids; ++i) {
      bool coarsen_yes = false;
      for(Ordinal j=0; j<dimension; ++j) {
        if(j<3) {
          if(((nLoc[j]-1)%2)==0 && ((nLoc[j]-1)/2 + 1)%2!=0 && (nLoc[j]-1)/2 + 1>3) {
            nLoc[j] = (nLoc[j]-1)/2 + 1;
            nGlo[j] = (nGlo[j]-1)/2 + 1;
            coarsen_yes = true;
          }
        } else if(!grid->getStencilWidths()->spectralT() && ((nLoc[j])%2)==0 && nLoc[j]>1) {
          nLoc[j] = (nLoc[j])/2;
          nGlo[j] = (nGlo[j])/2;
          coarsen_yes = true;
        }
      }
      if(coarsen_yes)
        multiGrid.push_back(createGrid(multiGrid.back(), nGlo));
    }


    // not working on brutus
    //multiGrid.shrink_to_fit();

    return multiGrid;
  }


}; // end of class CoarsenStrategy


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COARSENSTRATEGY_HPP
